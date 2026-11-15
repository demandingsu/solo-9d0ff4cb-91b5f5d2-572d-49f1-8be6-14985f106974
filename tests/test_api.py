"""API 端到端测试：上传→暂存→发布、等价归并、版本乐观锁、坏行隔离、并发发布。"""
import threading
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

SAMPLES = Path(__file__).resolve().parent.parent / "samples"


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def create_dataset(client, name="ds") -> int:
    r = client.post("/datasets", json={"name": name})
    assert r.status_code == 201, r.text
    return r.json()["id"]


def upload(client, dataset_id: int, vcf_name: str):
    with open(SAMPLES / "chr1.fa", "rb") as fa, open(SAMPLES / vcf_name, "rb") as vcf:
        return client.post(
            f"/datasets/{dataset_id}/uploads",
            files={"fasta": ("chr1.fa", fa, "text/plain"), "vcf": (vcf_name, vcf, "text/plain")},
        )


def publish(client, dataset_id: int, job_id: int, expected_version: int):
    return client.post(
        f"/datasets/{dataset_id}/publish",
        json={"job_id": job_id, "expected_version": expected_version},
    )


def core(variants):
    """忽略行号，只比较 (pos, ref, alt)。"""
    return [(v["pos"], v["ref"], v["alt"]) for v in variants]


def test_full_flow_and_equivalent_merge(client):
    did = create_dataset(client)

    # 上传批次 1：任务就绪，但已发布数据仍是空旧版 v0
    r = upload(client, did, "batch1.vcf")
    assert r.status_code == 201, r.text
    job1 = r.json()
    assert job1["status"] == "ready"
    assert job1["total_records"] == 6 and job1["bad_count"] == 0

    detail = client.get(f"/datasets/{did}").json()
    assert detail["version"] == 0 and detail["variants"] == []

    # 等价位点归并：6 行归并为 5 个位点，TA 重复区插入归并了两行（物理行号 5、6）
    staged = {tuple([v["pos"], v["ref"], v["alt"]]): v["lines"] for v in job1["variants"]}
    assert staged == {
        (4, "G", "GTA"): [5, 6],
        (4, "GTA", "G"): [7],
        (17, "A", "G"): [8],
        (18, "T", "G"): [9],
        (20, "G", "GA"): [10],
    }

    # 发布 v0 -> v1
    r = publish(client, did, job1["id"], expected_version=0)
    assert r.status_code == 200, r.text
    assert r.json()["version"] == 1

    d1 = client.get(f"/datasets/{did}").json()
    assert d1["version"] == 1 and d1["chrom"] == "chr1"
    assert core(d1["variants"]) == sorted(staged.keys())

    # 批次 2：同一批变异的不同写法
    job2 = upload(client, did, "batch2.vcf").json()
    assert job2["status"] == "ready"

    # 过期版本发布 -> 409，旧版不变
    r = publish(client, did, job2["id"], expected_version=0)
    assert r.status_code == 409
    assert r.json()["detail"]["current_version"] == 1
    assert client.get(f"/datasets/{did}").json()["version"] == 1

    # 正确版本发布 v1 -> v2
    r = publish(client, did, job2["id"], expected_version=1)
    assert r.status_code == 200, r.text
    assert r.json()["version"] == 2

    # 两批不同写法的 VCF 归一化结果完全一致（仅原始行号不同）
    d2 = client.get(f"/datasets/{did}").json()
    assert core(d2["variants"]) == core(d1["variants"])
    lines2 = {tuple([v["pos"], v["ref"], v["alt"]]): v["lines"] for v in d2["variants"]}
    assert lines2[(4, "G", "GTA")] == [5, 6, 7]
    assert lines2[(4, "GTA", "G")] == [8, 9]


def test_bad_rows_fail_job_and_never_touch_published(client):
    did = create_dataset(client)
    job1 = upload(client, did, "batch1.vcf").json()
    assert publish(client, did, job1["id"], 0).status_code == 200
    before = client.get(f"/datasets/{did}").json()

    # 含 5 条坏行：任务 failed，错误逐行给出
    r = upload(client, did, "bad.vcf")
    assert r.status_code == 201, r.text
    bad_job = r.json()
    assert bad_job["status"] == "failed"
    assert bad_job["total_records"] == 5 and bad_job["bad_count"] == 5
    reasons = {e["line"]: e["reason"] for e in bad_job["errors"]}
    assert "REF 与参考不符" in reasons[5]
    assert "双等位" in reasons[6]
    assert "MNP" in reasons[7]
    assert "超出染色体末端" in reasons[8]
    assert "不一致" in reasons[9]

    # 失败任务不可发布 -> 422，已发布数据不变
    r = publish(client, did, bad_job["id"], expected_version=1)
    assert r.status_code == 422
    after = client.get(f"/datasets/{did}").json()
    assert after["version"] == before["version"] == 1
    assert core(after["variants"]) == core(before["variants"])


def test_concurrent_publish_only_version_match_wins(client):
    did = create_dataset(client)
    j1 = upload(client, did, "batch1.vcf").json()["id"]
    j2 = upload(client, did, "batch2.vcf").json()["id"]

    barrier = threading.Barrier(3)
    statuses = []

    def pub(job_id):
        barrier.wait(timeout=10)
        statuses.append(publish(client, did, job_id, expected_version=0).status_code)

    threads = [threading.Thread(target=pub, args=(j,)) for j in (j1, j2)]
    for t in threads:
        t.start()
    barrier.wait(timeout=10)
    for t in threads:
        t.join(timeout=30)

    assert sorted(statuses) == [200, 409]
    assert client.get(f"/datasets/{did}").json()["version"] == 1


def test_not_found_and_validation(client):
    assert client.get("/datasets/999").status_code == 404
    assert upload(client, 999, "batch1.vcf").status_code == 404
    did = create_dataset(client)
    assert publish(client, did, 999, 0).status_code == 404
    # 文件级致命错误 -> 400
    r = client.post(
        f"/datasets/{did}/uploads",
        files={
            "fasta": ("ref.fa", b">chr1\nACGT\n>chr2\nACGT\n", "text/plain"),
            "vcf": ("v.vcf", (SAMPLES / "batch1.vcf").read_bytes(), "text/plain"),
        },
    )
    assert r.status_code == 400
    assert "单染色体" in r.json()["detail"]


def test_get_job_detail(client):
    did = create_dataset(client)
    job = upload(client, did, "batch1.vcf").json()
    r = client.get(f"/datasets/{did}/jobs/{job['id']}")
    assert r.status_code == 200
    assert r.json()["status"] == "ready"
    assert len(r.json()["variants"]) == 5
    assert client.get(f"/datasets/{did}/jobs/999").status_code == 404
