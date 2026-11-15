#!/usr/bin/env python3
"""端到端演示：等价变异归并、版本乐观锁、坏行隔离。

用法：
    python3 scripts/demo.py [BASE_URL]     # 默认 http://localhost:8000

前置：服务已启动（如 docker compose up --build）。
"""
import json
import sys
import urllib.error
import urllib.request
from pathlib import Path

BASE = sys.argv[1].rstrip("/") if len(sys.argv) > 1 else "http://localhost:8000"
SAMPLES = Path(__file__).resolve().parent.parent / "samples"


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


def upload(vcf_name):
    return req(
        "POST",
        f"/datasets/{DATASET_ID}/uploads",
        files={
            "fasta": ("chr1.fa", (SAMPLES / "chr1.fa").read_bytes()),
            "vcf": (vcf_name, (SAMPLES / vcf_name).read_bytes()),
        },
    )


def core(variants):
    return [(v["pos"], v["ref"], v["alt"]) for v in variants]


def main():
    global DATASET_ID

    status, _ = req("GET", "/healthz")
    if status != 200:
        sys.exit(f"服务不可达：{BASE}")

    status, ds = req("POST", "/datasets", {"name": "demo-等价变异"})
    show("1. 创建数据集", status, ds)
    DATASET_ID = ds["id"]

    status, job1 = upload("batch1.vcf")
    show("2. 上传批次 1（生成暂存任务）", status, job1)

    status, d = req("GET", f"/datasets/{DATASET_ID}")
    show("3. 处理期间/发布前查询：仍是完整旧版 v0", status,
         {"version": d["version"], "variant_count": d["variant_count"]})

    status, pub1 = req("POST", f"/datasets/{DATASET_ID}/publish",
                       {"job_id": job1["id"], "expected_version": 0})
    show("4. 按查询所得版本发布：v0 -> v1", status, pub1)

    status, d1 = req("GET", f"/datasets/{DATASET_ID}")
    show("5. 查询 v1 归一化结果", status, d1)

    status, job2 = upload("batch2.vcf")
    show("6. 上传批次 2（同一批变异的不同写法）", status, job2)

    status, conflict = req("POST", f"/datasets/{DATASET_ID}/publish",
                           {"job_id": job2["id"], "expected_version": 0})
    show("7. 用过期版本发布：应 409，旧版不变", status, conflict)

    status, pub2 = req("POST", f"/datasets/{DATASET_ID}/publish",
                       {"job_id": job2["id"], "expected_version": 1})
    show("8. 用当前版本发布：v1 -> v2", status, pub2)

    status, d2 = req("GET", f"/datasets/{DATASET_ID}")
    same = core(d1["variants"]) == core(d2["variants"])
    print(f"\n>>> 两批写法不同的 VCF 归一化结果完全一致：{same}")
    print("    批次1 行号 vs 批次2 行号（TA 重复区插入）：",
          d1["variants"][0]["lines"], "vs", d2["variants"][0]["lines"])
    assert same

    status, job3 = upload("bad.vcf")
    show("9. 上传含坏行 VCF：任务 failed，逐行给出原因", status, job3)

    status, rejected = req("POST", f"/datasets/{DATASET_ID}/publish",
                           {"job_id": job3["id"], "expected_version": 2})
    show("10. 发布失败任务：应 422，旧版不变", status, rejected)

    status, d3 = req("GET", f"/datasets/{DATASET_ID}")
    show("11. 最终查询：仍为 v2，内容未受坏行影响", status,
         {"version": d3["version"], "variant_count": d3["variant_count"]})

    assert d3["version"] == 2 and core(d3["variants"]) == core(d2["variants"])
    print("\n演示完成：等价变异已归并，坏行与版本冲突均未影响已发布数据。")


if __name__ == "__main__":
    DATASET_ID = None
    main()
