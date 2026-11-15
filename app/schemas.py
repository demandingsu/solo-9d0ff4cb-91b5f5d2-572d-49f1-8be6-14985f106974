"""API 请求/响应模型。"""
from __future__ import annotations

from pydantic import BaseModel, Field


class DatasetCreate(BaseModel):
    name: str = Field(min_length=1, max_length=200, description="数据集名称")


class VariantOut(BaseModel):
    pos: int = Field(description="归一化后的 1-based 坐标")
    ref: str = Field(description="归一化后的 REF")
    alt: str = Field(description="归一化后的 ALT")
    lines: list[int] = Field(description="归并到该位点的所有原始 VCF 行号（含头部行计数，从 1 起）")


class DatasetOut(BaseModel):
    id: int
    name: str
    version: int = Field(description="当前已发布版本，0 表示尚未发布")
    published_job_id: int | None
    created_at: str
    updated_at: str


class DatasetDetail(DatasetOut):
    chrom: str | None = Field(description="已发布数据所属染色体")
    variant_count: int
    variants: list[VariantOut]


class RowErrorOut(BaseModel):
    line: int = Field(description="原始 VCF 物理行号（从 1 起，含头部行）")
    reason: str
    chrom: str | None = None
    pos: str | None = None
    ref: str | None = None
    alt: str | None = None


class JobOut(BaseModel):
    id: int
    dataset_id: int
    status: str = Field(description="ready=可发布；failed=含坏行，不可发布")
    fasta_name: str
    vcf_name: str
    chrom: str | None
    total_records: int = Field(description="VCF 数据行总数")
    bad_count: int = Field(description="坏行数")
    errors: list[RowErrorOut]
    created_at: str


class JobDetail(JobOut):
    variants: list[VariantOut] | None = Field(
        default=None, description="暂存的归一化结果（仅 ready 任务存在，发布前不影响已发布数据）"
    )


class PublishRequest(BaseModel):
    job_id: int = Field(description="待发布的任务 ID")
    expected_version: int = Field(description="查询所得的当前版本号；不匹配则拒绝发布（乐观锁）")


class PublishResponse(BaseModel):
    id: int
    version: int = Field(description="发布后的新版本号")
    published_job_id: int
    variant_count: int
