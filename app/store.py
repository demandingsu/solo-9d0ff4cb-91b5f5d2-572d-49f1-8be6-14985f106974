"""SQLite 存储层。

数据模型：
- datasets：数据集，version 从 0 起，每次成功发布 +1；
- jobs / job_variants：上传处理后的暂存结果（不影响已发布数据）；
- variants：当前已发布的归一化位点（查询只读这张表）；
- reports / report_items：批次交集报告（报告头 + 共同位点行），持久保存。

发布语义：单事务内校验 expected_version 并整体替换 variants，
版本不匹配或任务含坏行时回滚，已发布数据保持不变。

报告语义：单事务（BEGIN IMMEDIATE）内读取两个数据集的当前版本与全部
已发布位点、计算交集并写入报告——读取期间并发发布被写锁阻塞，报告
必然来自同一时点的两份完整已发布数据；任何失败整体回滚，不留半份报告；
报告提交后与已发布数据解耦，后续替换不影响按报告号读取。
"""
from __future__ import annotations

import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

SCHEMA = """
CREATE TABLE IF NOT EXISTS datasets (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    name             TEXT NOT NULL,
    version          INTEGER NOT NULL DEFAULT 0,
    published_job_id INTEGER,
    created_at       TEXT NOT NULL,
    updated_at       TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS jobs (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_id    INTEGER NOT NULL REFERENCES datasets(id),
    status        TEXT NOT NULL CHECK (status IN ('ready', 'failed')),
    fasta_name    TEXT NOT NULL,
    vcf_name      TEXT NOT NULL,
    chrom         TEXT,
    total_records INTEGER NOT NULL,
    bad_count     INTEGER NOT NULL,
    errors_json   TEXT NOT NULL DEFAULT '[]',
    created_at    TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_jobs_dataset ON jobs(dataset_id);

CREATE TABLE IF NOT EXISTS job_variants (
    job_id     INTEGER NOT NULL REFERENCES jobs(id),
    chrom      TEXT NOT NULL,
    pos        INTEGER NOT NULL,
    ref        TEXT NOT NULL,
    alt        TEXT NOT NULL,
    lines_json TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_job_variants_job ON job_variants(job_id);

CREATE TABLE IF NOT EXISTS variants (
    dataset_id INTEGER NOT NULL REFERENCES datasets(id),
    chrom      TEXT NOT NULL,
    pos        INTEGER NOT NULL,
    ref        TEXT NOT NULL,
    alt        TEXT NOT NULL,
    lines_json TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_variants_dataset ON variants(dataset_id);

CREATE TABLE IF NOT EXISTS reports (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_a_id INTEGER NOT NULL REFERENCES datasets(id),
    dataset_b_id INTEGER NOT NULL REFERENCES datasets(id),
    version_a    INTEGER NOT NULL,
    version_b    INTEGER NOT NULL,
    chrom        TEXT NOT NULL,
    common_count INTEGER NOT NULL,
    created_at   TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS report_items (
    report_id    INTEGER NOT NULL REFERENCES reports(id),
    pos          INTEGER NOT NULL,
    ref          TEXT NOT NULL,
    alt          TEXT NOT NULL,
    lines_a_json TEXT NOT NULL,
    lines_b_json TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_report_items_report ON report_items(report_id);
"""


def _utcnow() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class NotFoundError(Exception):
    """数据集或任务不存在。"""


class VersionConflictError(Exception):
    """发布时版本不匹配（乐观锁冲突）。"""

    def __init__(self, current_version: int):
        super().__init__(f"版本不匹配：当前版本为 {current_version}")
        self.current_version = current_version


class JobNotReadyError(Exception):
    """任务含坏行，不可发布。"""

    def __init__(self, status: str, errors: list[dict]):
        super().__init__("任务含坏行，不可发布")
        self.status = status
        self.errors = errors


class NotPublishedError(Exception):
    """数据集尚未发布，无法生成交集报告。"""


class ChromMismatchError(Exception):
    """两个数据集已发布数据的染色体不同（或缺失），无法比较。"""


class Store:
    def __init__(self, path: str):
        self.path = path
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        conn = self._connect()
        try:
            conn.executescript(SCHEMA)
            conn.execute("PRAGMA journal_mode=WAL")  # 写事务期间允许并发读旧版本
        finally:
            conn.close()

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.path, timeout=30, check_same_thread=False, isolation_level=None)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        conn.execute("PRAGMA busy_timeout = 30000")
        return conn

    # ---------- 数据集 ----------

    def create_dataset(self, name: str) -> dict:
        now = _utcnow()
        conn = self._connect()
        try:
            cur = conn.execute(
                "INSERT INTO datasets (name, version, created_at, updated_at) VALUES (?, 0, ?, ?)",
                (name, now, now),
            )
            row = conn.execute("SELECT * FROM datasets WHERE id = ?", (cur.lastrowid,)).fetchone()
            return self._dataset_dict(row)
        finally:
            conn.close()

    def list_datasets(self) -> list[dict]:
        conn = self._connect()
        try:
            rows = conn.execute("SELECT * FROM datasets ORDER BY id").fetchall()
            return [self._dataset_dict(r) for r in rows]
        finally:
            conn.close()

    def get_dataset_detail(self, dataset_id: int) -> dict | None:
        conn = self._connect()
        try:
            row = conn.execute("SELECT * FROM datasets WHERE id = ?", (dataset_id,)).fetchone()
            if row is None:
                return None
            vrows = conn.execute(
                "SELECT chrom, pos, ref, alt, lines_json FROM variants "
                "WHERE dataset_id = ? ORDER BY pos, ref, alt",
                (dataset_id,),
            ).fetchall()
            detail = self._dataset_dict(row)
            detail["chrom"] = vrows[0]["chrom"] if vrows else None
            detail["variant_count"] = len(vrows)
            detail["variants"] = [self._variant_dict(v) for v in vrows]
            return detail
        finally:
            conn.close()

    # ---------- 任务（上传暂存） ----------

    def create_job(
        self,
        dataset_id: int,
        fasta_name: str,
        vcf_name: str,
        chrom: str,
        total_records: int,
        errors: list,
        variants: list,
    ) -> dict:
        """把上传处理结果整体写入暂存区（单事务），不触碰已发布数据。"""
        status = "ready" if not errors else "failed"
        now = _utcnow()
        conn = self._connect()
        try:
            conn.execute("BEGIN IMMEDIATE")
            if conn.execute("SELECT 1 FROM datasets WHERE id = ?", (dataset_id,)).fetchone() is None:
                raise NotFoundError(f"数据集 {dataset_id} 不存在")
            cur = conn.execute(
                "INSERT INTO jobs (dataset_id, status, fasta_name, vcf_name, chrom, "
                "total_records, bad_count, errors_json, created_at) "
                "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                (
                    dataset_id,
                    status,
                    fasta_name,
                    vcf_name,
                    chrom,
                    total_records,
                    len(errors),
                    json.dumps([e.to_dict() for e in errors], ensure_ascii=False),
                    now,
                ),
            )
            job_id = cur.lastrowid
            if status == "ready":
                conn.executemany(
                    "INSERT INTO job_variants (job_id, chrom, pos, ref, alt, lines_json) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    [(job_id, chrom, v.pos, v.ref, v.alt, json.dumps(v.lines)) for v in variants],
                )
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()
        job = self.get_job(dataset_id, job_id)
        assert job is not None
        return job

    def get_job(self, dataset_id: int, job_id: int) -> dict | None:
        conn = self._connect()
        try:
            row = conn.execute(
                "SELECT * FROM jobs WHERE id = ? AND dataset_id = ?", (job_id, dataset_id)
            ).fetchone()
            if row is None:
                return None
            variants = None
            if row["status"] == "ready":
                vrows = conn.execute(
                    "SELECT pos, ref, alt, lines_json FROM job_variants "
                    "WHERE job_id = ? ORDER BY pos, ref, alt",
                    (job_id,),
                ).fetchall()
                variants = [self._variant_dict(v) for v in vrows]
            return self._job_dict(row, variants)
        finally:
            conn.close()

    # ---------- 发布（乐观版本锁 + 整体替换） ----------

    def publish(self, dataset_id: int, job_id: int, expected_version: int) -> dict:
        """按期望版本发布任务结果。

        在 BEGIN IMMEDIATE 写事务内完成全部校验与替换：
        - 版本不匹配 -> VersionConflictError（409）；
        - 任务含坏行 -> JobNotReadyError（422）；
        任何失败都回滚，已发布数据保持不变。
        """
        conn = self._connect()
        try:
            conn.execute("BEGIN IMMEDIATE")

            ds = conn.execute("SELECT * FROM datasets WHERE id = ?", (dataset_id,)).fetchone()
            if ds is None:
                raise NotFoundError(f"数据集 {dataset_id} 不存在")
            job = conn.execute(
                "SELECT * FROM jobs WHERE id = ? AND dataset_id = ?", (job_id, dataset_id)
            ).fetchone()
            if job is None:
                raise NotFoundError(f"任务 {job_id} 不存在")

            if ds["version"] != expected_version:
                raise VersionConflictError(ds["version"])
            if job["status"] != "ready":
                raise JobNotReadyError(job["status"], json.loads(job["errors_json"]))

            conn.execute("DELETE FROM variants WHERE dataset_id = ?", (dataset_id,))
            conn.execute(
                "INSERT INTO variants (dataset_id, chrom, pos, ref, alt, lines_json) "
                "SELECT ?, chrom, pos, ref, alt, lines_json FROM job_variants WHERE job_id = ?",
                (dataset_id, job_id),
            )
            cur = conn.execute(
                "UPDATE datasets SET version = version + 1, published_job_id = ?, updated_at = ? "
                "WHERE id = ? AND version = ?",
                (job_id, _utcnow(), dataset_id, expected_version),
            )
            if cur.rowcount != 1:  # 持写锁期间理论上不会发生，兜底
                raise VersionConflictError(expected_version)

            conn.commit()
            count = conn.execute(
                "SELECT COUNT(*) AS c FROM variants WHERE dataset_id = ?", (dataset_id,)
            ).fetchone()["c"]
            return {
                "id": dataset_id,
                "version": expected_version + 1,
                "published_job_id": job_id,
                "variant_count": count,
            }
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()

    # ---------- 批次交集报告 ----------

    def create_report(self, dataset_a_id: int, dataset_b_id: int) -> dict:
        """生成两个数据集当前已发布位点的交集报告。

        在单个 BEGIN IMMEDIATE 事务内完成全部读取与写入：
        - 读取期间并发发布（同样需要写锁）被阻塞，两份已发布数据来自
          同一时点的一致快照，版本号与位点内容不会跨版本撕裂；
        - 交集计算与报告写入在同一事务提交，任何失败整体回滚，
          不会留下半份报告；
        - 报告提交后与已发布数据解耦，后续替换不影响按报告号读取。

        拒绝情形：数据集不存在（NotFoundError）、尚未发布
        （NotPublishedError）、染色体不同或缺失（ChromMismatchError）。
        """
        conn = self._connect()
        try:
            conn.execute("BEGIN IMMEDIATE")
            snap_a = self._published_snapshot(conn, dataset_a_id)
            snap_b = self._published_snapshot(conn, dataset_b_id)

            chrom_a, chrom_b = snap_a["chrom"], snap_b["chrom"]
            if chrom_a is None or chrom_b is None:
                raise ChromMismatchError("已发布数据缺少染色体信息，无法比较")
            if chrom_a != chrom_b:
                raise ChromMismatchError(
                    f"染色体不同：数据集 {dataset_a_id} 为 {chrom_a}，"
                    f"数据集 {dataset_b_id} 为 {chrom_b}"
                )

            items = self._intersect(snap_a["variants"], snap_b["variants"])
            report_id = self._insert_report(
                conn,
                dataset_a_id,
                dataset_b_id,
                snap_a["version"],
                snap_b["version"],
                chrom_a,
                items,
            )
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            conn.close()
        report = self.get_report(report_id)
        assert report is not None
        return report

    def get_report(self, report_id: int) -> dict | None:
        """按报告号读取报告（含共同位点明细）；不存在返回 None。"""
        conn = self._connect()
        try:
            row = conn.execute("SELECT * FROM reports WHERE id = ?", (report_id,)).fetchone()
            if row is None:
                return None
            irows = conn.execute(
                "SELECT pos, ref, alt, lines_a_json, lines_b_json FROM report_items "
                "WHERE report_id = ? ORDER BY pos, ref, alt",
                (report_id,),
            ).fetchall()
            return self._report_dict(row, irows)
        finally:
            conn.close()

    def list_reports(self) -> list[dict]:
        """列出全部报告（不含共同位点明细）。"""
        conn = self._connect()
        try:
            rows = conn.execute("SELECT * FROM reports ORDER BY id").fetchall()
            return [self._report_dict(r, None) for r in rows]
        finally:
            conn.close()

    # ---------- 报告内部辅助 ----------

    def _published_snapshot(self, conn: sqlite3.Connection, dataset_id: int) -> dict:
        """在调用方事务内读取数据集当前已发布版本、染色体与全部位点。"""
        ds = conn.execute("SELECT * FROM datasets WHERE id = ?", (dataset_id,)).fetchone()
        if ds is None:
            raise NotFoundError(f"数据集 {dataset_id} 不存在")
        if ds["version"] == 0:
            raise NotPublishedError(f"数据集 {dataset_id} 尚未发布")
        vrows = conn.execute(
            "SELECT chrom, pos, ref, alt, lines_json FROM variants "
            "WHERE dataset_id = ? ORDER BY pos, ref, alt",
            (dataset_id,),
        ).fetchall()
        chrom = vrows[0]["chrom"] if vrows else None
        if chrom is None:  # 发布了空 VCF 的边缘情况：染色体取自任务记录
            job = conn.execute(
                "SELECT chrom FROM jobs WHERE id = ?", (ds["published_job_id"],)
            ).fetchone()
            chrom = job["chrom"] if job else None
        return {
            "version": ds["version"],
            "chrom": chrom,
            "variants": [self._variant_dict(v) for v in vrows],
        }

    @staticmethod
    def _intersect(variants_a: list[dict], variants_b: list[dict]) -> list[dict]:
        """按归一化后的 (pos, ref, alt) 求交集，保留双方各自的原始行号。

        variants_a 已按 (pos, ref, alt) 排序，交集结果保持同样顺序。
        """
        lines_b = {(v["pos"], v["ref"], v["alt"]): v["lines"] for v in variants_b}
        items = []
        for v in variants_a:
            key = (v["pos"], v["ref"], v["alt"])
            if key in lines_b:
                items.append(
                    {
                        "pos": v["pos"],
                        "ref": v["ref"],
                        "alt": v["alt"],
                        "lines_a": v["lines"],
                        "lines_b": lines_b[key],
                    }
                )
        return items

    @staticmethod
    def _insert_report(
        conn: sqlite3.Connection,
        dataset_a_id: int,
        dataset_b_id: int,
        version_a: int,
        version_b: int,
        chrom: str,
        items: list[dict],
    ) -> int:
        """在调用方事务内写入报告头与全部报告行（随同事务提交或回滚）。"""
        cur = conn.execute(
            "INSERT INTO reports (dataset_a_id, dataset_b_id, version_a, version_b, "
            "chrom, common_count, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (dataset_a_id, dataset_b_id, version_a, version_b, chrom, len(items), _utcnow()),
        )
        report_id = cur.lastrowid
        conn.executemany(
            "INSERT INTO report_items (report_id, pos, ref, alt, lines_a_json, lines_b_json) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            [
                (
                    report_id,
                    it["pos"],
                    it["ref"],
                    it["alt"],
                    json.dumps(it["lines_a"]),
                    json.dumps(it["lines_b"]),
                )
                for it in items
            ],
        )
        return report_id

    # ---------- 行转换 ----------

    @staticmethod
    def _dataset_dict(row: sqlite3.Row) -> dict:
        return {
            "id": row["id"],
            "name": row["name"],
            "version": row["version"],
            "published_job_id": row["published_job_id"],
            "created_at": row["created_at"],
            "updated_at": row["updated_at"],
        }

    @staticmethod
    def _variant_dict(row: sqlite3.Row) -> dict:
        return {
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
            "lines": json.loads(row["lines_json"]),
        }

    @staticmethod
    def _job_dict(row: sqlite3.Row, variants: list[dict] | None = None) -> dict:
        d = {
            "id": row["id"],
            "dataset_id": row["dataset_id"],
            "status": row["status"],
            "fasta_name": row["fasta_name"],
            "vcf_name": row["vcf_name"],
            "chrom": row["chrom"],
            "total_records": row["total_records"],
            "bad_count": row["bad_count"],
            "errors": json.loads(row["errors_json"]),
            "created_at": row["created_at"],
        }
        if variants is not None:
            d["variants"] = variants
        return d

    @staticmethod
    def _report_dict(row: sqlite3.Row, item_rows: list[sqlite3.Row] | None) -> dict:
        d = {
            "id": row["id"],
            "dataset_a_id": row["dataset_a_id"],
            "dataset_b_id": row["dataset_b_id"],
            "version_a": row["version_a"],
            "version_b": row["version_b"],
            "chrom": row["chrom"],
            "common_count": row["common_count"],
            "created_at": row["created_at"],
        }
        if item_rows is not None:
            d["items"] = [
                {
                    "chrom": row["chrom"],
                    "pos": r["pos"],
                    "ref": r["ref"],
                    "alt": r["alt"],
                    "lines_a": json.loads(r["lines_a_json"]),
                    "lines_b": json.loads(r["lines_b_json"]),
                }
                for r in item_rows
            ]
        return d
