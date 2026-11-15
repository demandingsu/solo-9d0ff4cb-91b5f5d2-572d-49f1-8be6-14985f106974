"""FastAPI 入口：VCF 变异归一化服务。"""
from __future__ import annotations

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI, File, HTTPException, Request, UploadFile

from .parse import ParseError
from .schemas import (
    DatasetCreate,
    DatasetDetail,
    DatasetOut,
    JobDetail,
    PublishRequest,
    PublishResponse,
    ReportCreate,
    ReportDetail,
    ReportOut,
)
from .service import process_upload
from .store import (
    ChromMismatchError,
    JobNotReadyError,
    NotFoundError,
    NotPublishedError,
    Store,
    VersionConflictError,
)

MAX_UPLOAD_BYTES = 100 * 1024 * 1024  # 单文件 100MB 上限


@asynccontextmanager
async def lifespan(app: FastAPI):
    app.state.store = Store(os.environ.get("DATABASE_PATH", "data/variants.db"))
    yield


app = FastAPI(
    title="VCF 变异归一化服务",
    version="1.0.0",
    lifespan=lifespan,
    description=(
        "面向测序中心数据员的跨批次变异归一化 API。\n\n"
        "**流程**：创建数据集 → 上传单染色体 FASTA + 双等位 SNV/插入/缺失 VCF（生成暂存任务）→ "
        "按查询所得版本号发布（乐观锁整体替换）。\n\n"
        "**归一化**：逐条核对 REF → 消去共有前后缀 → 保留锚点 → 串联重复左对齐 → 归并等价位点（保留各原始行号）。\n\n"
        "**批次交集报告**：提交两个数据集号，按归一化后的染色体、坐标、REF/ALT 找出共同位点，"
        "保存双方原始行号与所用版本，返回报告号供复核。\n\n"
        "**一致性**：处理期间查询仍返回完整旧版；并发发布仅版本匹配者成功；坏行或发布失败不改动旧版；"
        "报告生成于单事务快照（并发替换不撕裂），保存失败不留半份，替换后旧报告仍可按号读取。"
    ),
)


def get_store(request: Request) -> Store:
    return request.app.state.store


@app.get("/", summary="服务信息")
def root() -> dict:
    return {"service": "VCF 变异归一化服务", "docs": "/docs", "health": "/healthz"}


@app.get("/healthz", summary="健康检查")
def healthz() -> dict:
    return {"status": "ok"}


@app.post("/datasets", status_code=201, response_model=DatasetOut, summary="创建数据集")
def create_dataset(body: DatasetCreate, request: Request):
    return get_store(request).create_dataset(body.name)


@app.get("/datasets", response_model=list[DatasetOut], summary="列出全部数据集")
def list_datasets(request: Request):
    return get_store(request).list_datasets()


@app.get(
    "/datasets/{dataset_id}",
    response_model=DatasetDetail,
    summary="查询数据集当前已发布版本（含归一化位点）",
)
def get_dataset(dataset_id: int, request: Request):
    detail = get_store(request).get_dataset_detail(dataset_id)
    if detail is None:
        raise HTTPException(status_code=404, detail=f"数据集 {dataset_id}  不存在")
    return detail


@app.post(
    "/datasets/{dataset_id}/uploads",
    status_code=201,
    response_model=JobDetail,
    summary="上传 FASTA + VCF 并归一化（生成暂存任务，不影响已发布数据）",
)
async def upload(
    dataset_id: int,
    request: Request,
    fasta: UploadFile = File(..., description="单染色体 FASTA"),
    vcf: UploadFile = File(..., description="仅含双等位 SNV/插入/缺失的 VCF"),
):
    fasta_bytes = await fasta.read()
    vcf_bytes = await vcf.read()
    if len(fasta_bytes) > MAX_UPLOAD_BYTES or len(vcf_bytes) > MAX_UPLOAD_BYTES:
        raise HTTPException(status_code=413, detail="文件超过 100MB 限制")
    try:
        fasta_text = fasta_bytes.decode("utf-8")
        vcf_text = vcf_bytes.decode("utf-8")
    except UnicodeDecodeError:
        raise HTTPException(status_code=400, detail="文件必须为 UTF-8 文本")

    try:
        chrom, total, variants, errors = process_upload(fasta_text, vcf_text)
    except ParseError as exc:
        raise HTTPException(status_code=400, detail=f"文件解析失败：{exc}")

    try:
        return get_store(request).create_job(
            dataset_id,
            fasta.filename or "ref.fa",
            vcf.filename or "variants.vcf",
            chrom,
            total,
            errors,
            variants,
        )
    except NotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc))


@app.get(
    "/datasets/{dataset_id}/jobs/{job_id}",
    response_model=JobDetail,
    summary="查询任务状态（坏行明细 / 暂存的归一化结果）",
)
def get_job(dataset_id: int, job_id: int, request: Request):
    job = get_store(request).get_job(dataset_id, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail=f"任务 {job_id} 不存在")
    return job


@app.post(
    "/datasets/{dataset_id}/publish",
    response_model=PublishResponse,
    summary="按期望版本发布任务结果（乐观锁整体替换）",
    responses={
        409: {"description": "版本不匹配，发布被拒绝，旧版本保持不变"},
        422: {"description": "任务含坏行，不可发布，旧版本保持不变"},
    },
)
def publish(dataset_id: int, body: PublishRequest, request: Request):
    try:
        return get_store(request).publish(dataset_id, body.job_id, body.expected_version)
    except NotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    except VersionConflictError as exc:
        raise HTTPException(
            status_code=409,
            detail={
                "message": "版本不匹配，发布被拒绝（旧版本保持不变）",
                "expected_version": body.expected_version,
                "current_version": exc.current_version,
            },
        )
    except JobNotReadyError as exc:
        raise HTTPException(
            status_code=422,
            detail={
                "message": "任务含坏行，不可发布（旧版本保持不变）",
                "status": exc.status,
                "errors": exc.errors,
            },
        )


@app.post(
    "/reports",
    status_code=201,
    response_model=ReportDetail,
    summary="生成两个数据集已发布变异的交集报告",
    responses={
        404: {"description": "数据集不存在"},
        409: {"description": "数据集尚未发布"},
        422: {"description": "两个数据集染色体不同（或缺失），无法比较"},
    },
)
def create_report(body: ReportCreate, request: Request):
    """按归一化后的染色体、坐标、REF/ALT 找出两个数据集的共同位点。

    报告在单事务快照内生成：生成期间若有并发替换，报告仍来自同一时点的
    两份完整已发布数据；保存失败整体回滚，不留半份报告。
    """
    try:
        return get_store(request).create_report(body.dataset_a_id, body.dataset_b_id)
    except NotFoundError as exc:
        raise HTTPException(status_code=404, detail=str(exc))
    except NotPublishedError as exc:
        raise HTTPException(status_code=409, detail=str(exc))
    except ChromMismatchError as exc:
        raise HTTPException(status_code=422, detail=str(exc))


@app.get("/reports", response_model=list[ReportOut], summary="列出全部交集报告（不含明细）")
def list_reports(request: Request):
    return get_store(request).list_reports()


@app.get(
    "/reports/{report_id}",
    response_model=ReportDetail,
    summary="按报告号读取交集报告（已发布数据被替换后仍可读）",
)
def get_report(report_id: int, request: Request):
    report = get_store(request).get_report(report_id)
    if report is None:
        raise HTTPException(status_code=404, detail=f"报告 {report_id} 不存在")
    return report
