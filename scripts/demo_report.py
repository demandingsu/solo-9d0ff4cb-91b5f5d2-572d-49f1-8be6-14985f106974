#!/usr/bin/env python3
"""批次交集报告端到端演示：共同位点、双方行号与所用版本、并发替换一致性、旧报告复核。

用法：
    python3 scripts/demo_report.py [BASE_URL]     # 默认 http://localhost:8000

前置：服务已启动（如 docker compose up --build）。
"""
import json
import sys
import threading
import urllib.error
import urllib.request
from pathlib import Path

BASE = sys.argv[1].rstrip("/") if len(sys.argv) > 1 else "http://localhost:8000"
SAMPLES = Path(__file__).resolve().parent.parent / "samples"

# batch1 / batch2 发布后（归一化）的位点 -> 原始行号，用于校验报告快照不撕裂
BATCH1_LINES = {
    (4, "G", "GTA"): [5, 6],
    (4, "GTA", "G"): [7],
    (17, "A", "G"): [8],
    (18, "T", "G"): [9],
    (20, "G", "GA"): [10],
}
BATCH2_LINES = {
    (4, "G", "GTA"): [5, 6, 7],
    (4, "GTA", "G"): [8, 9],
    (17, "A", "G"): [10],
    (18, "T", "G"): [11],
    (20, "G", "GA"): [12],
}


def req(method, path, body=None, files=None):
    """最小 HTTP 客户端（仅用标准库）。files: {字段名: (文件名, bytes)}。"""
    url = BASE + path
    if files:
        boundary = "----demo-boundary-7d1a9f3c"
        payload = b""
        for field, (filename, data) in files.items():
            payload += (
                f"--{boundary}\r\n"
                f'Content-Disposition: form-data; name="{field}"; filename="{filename}"\r\n'
                f"Content-Type: application/octet-stream\r\n\r\n"
            ).encode() + data + b"\r\n"
        payload += f"--{boundary}--\r\n".encode()
        request = urllib.request.Request(
            url, data=payload, method="POST",
            headers={"Content-Type": f"multipart/form-data; boundary={boundary}"},
        )
    else:
        data = json.dumps(body).encode() if body is not None else None
        request = urllib.request.Request(
            url, data=data, method=method, headers={"Content-Type": "application/json"}
        )
    try:
        with urllib.request.urlopen(request) as resp:
            return resp.status, json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        return exc.code, json.loads(exc.read().decode())


def show(title, status, payload):
    print(f"\n=== {title}（HTTP {status}）===")
    print(json.dumps(payload, ensure_ascii=False, indent=2))


def upload(dataset_id, vcf_name, fasta_name="chr1.fa"):
    return req(
        "POST",
        f"/datasets/{dataset_id}/uploads",
        files={
            "fasta": (fasta_name, (SAMPLES / fasta_name).read_bytes()),
            "vcf": (vcf_name, (SAMPLES / vcf_name).read_bytes()),
        },
    )


def publish_batch(dataset_id, vcf_name, expected_version, fasta_name="chr1.fa"):
    status, job = upload(dataset_id, vcf_name, fasta_name)
    assert status == 201 and job["status"] == "ready", job
    status, pub = req("POST", f"/datasets/{dataset_id}/publish",
                      {"job_id": job["id"], "expected_version": expected_version})
    assert status == 200, pub
    return pub


def make_report(a, b):
    return req("POST", "/reports", {"dataset_a_id": a, "dataset_b_id": b})


def check_snapshot(report):
    """校验单份报告内部一致：version_a 与 lines_a 必须来自同一个已发布版本。"""
    expected = BATCH1_LINES if report["version_a"] % 2 == 1 else BATCH2_LINES
    got = {(it["pos"], it["ref"], it["alt"]): it["lines_a"] for it in report["items"]}
    return got == expected


def main():
    status, _ = req("GET", "/healthz")
    if status != 200:
        sys.exit(f"服务不可达：{BASE}")

    # ---- 1. 准备两个批次的数据集并发布 ----
    status, a = req("POST", "/datasets", {"name": "qc-批次1"})
    status, b = req("POST", "/datasets", {"name": "qc-批次2"})
    a_id, b_id = a["id"], b["id"]
    publish_batch(a_id, "batch1.vcf", 0)
    publish_batch(b_id, "batch2.vcf", 0)
    print(f"\n=== 1. 两个批次已发布：数据集 {a_id} v1（batch1）、数据集 {b_id} v1（batch2）===")

    # ---- 2. 生成交集报告 ----
    status, report = make_report(a_id, b_id)
    show("2. 生成交集报告：共同位点 + 双方原始行号 + 所用版本", status, report)
    assert status == 201 and report["common_count"] == 5
    report_id = report["id"]

    # ---- 3. 并发替换演示：报告生成期间，数据集 A 被反复重新发布 ----
    print("\n=== 3. 并发替换演示：后台反复重新发布数据集 A，同时连续生成报告 ===")
    stop = threading.Event()
    errors = []

    def republish():
        i = 0
        while not stop.is_set():
            vcf = "batch2.vcf" if i % 2 == 0 else "batch1.vcf"
            _, d = req("GET", f"/datasets/{a_id}")
            try:
                publish_batch(a_id, vcf, d["version"])
                i += 1
            except AssertionError as exc:
                errors.append(exc)

    worker = threading.Thread(target=republish)
    worker.start()
    concurrent_reports = []
    try:
        for _ in range(6):
            status, rep = make_report(a_id, b_id)
            assert status == 201, rep
            concurrent_reports.append(rep)
    finally:
        stop.set()
        worker.join(timeout=30)

    assert not errors, errors
    _, final_a = req("GET", f"/datasets/{a_id}")
    print(f"    期间数据集 A 被替换到 v{final_a['version']}，生成 {len(concurrent_reports)} 份报告：")
    for rep in concurrent_reports:
        ok = check_snapshot(rep)
        print(f"    报告 {rep['id']}：version_a={rep['version_a']} version_b={rep['version_b']} "
              f"common={rep['common_count']} 快照一致={ok}")
        assert ok, f"报告 {rep['id']} 快照撕裂"
    print("    >>> 每份报告的版本号与位点内容均来自同一时点的两份完整已发布数据，无撕裂。")

    # ---- 4. 替换后旧报告仍可按号读取 ----
    status, reread = req("GET", f"/reports/{report_id}")
    same = reread == report
    show(f"4. 替换后按报告号复核旧报告 {report_id}：内容不变={same}", status,
         {"version_a": reread["version_a"], "version_b": reread["version_b"],
          "common_count": reread["common_count"], "首条明细": reread["items"][0]})
    assert status == 200 and same

    # ---- 5. 拒绝场景 ----
    status, ds_new = req("POST", "/datasets", {"name": "qc-未发布"})
    status, rej = make_report(a_id, ds_new["id"])
    show("5a. 数据集尚未发布：应 409", status, rej)
    assert status == 409

    status, c = req("POST", "/datasets", {"name": "qc-chr2批次"})
    publish_batch(c["id"], "batch3_chr2.vcf", 0, fasta_name="chr2.fa")
    status, rej = make_report(a_id, c["id"])
    show("5b. 染色体不同（chr1 vs chr2）：应 422", status, rej)
    assert status == 422

    status, rej = make_report(a_id, 999999)
    show("5c. 数据集不存在：应 404", status, rej)
    assert status == 404

    status, reports = req("GET", "/reports")
    print(f"\n=== 6. 报告列表（{len(reports)} 份，均持久保存供复核）===")
    for r in reports:
        print(f"    报告 {r['id']}：数据集 {r['dataset_a_id']}@v{r['version_a']} "
              f"x 数据集 {r['dataset_b_id']}@v{r['version_b']} -> 共同位点 {r['common_count']}")

    print("\n演示完成：交集报告快照一致、原子保存、替换后旧报告可按号复核。")


if __name__ == "__main__":
    main()
