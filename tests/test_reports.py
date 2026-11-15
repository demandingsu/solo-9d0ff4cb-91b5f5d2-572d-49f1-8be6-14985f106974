"""批次交集报告 API 测试：生成、复核读取、拒绝场景、并发替换一致性、保存原子性。"""
import threading
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

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


def upload(client, dataset_id: int, vcf_name: str, fasta_name: str = "chr1.fa"):
    with open(SAMPLES / fasta_name, "rb") as fa, open(SAMPLES / vcf_name, "rb") as vcf:
        return client.post(
            f"/datasets/{dataset_id}/uploads",
            files={"fasta": (fasta_name, fa, "text/plain"), "vcf": (vcf_name, vcf, "text/plain")},
        )


def publish(client, dataset_id: int, job_id: int, expected_version: int):
    return client.post(
        f"/datasets/{dataset_id}/publish",
        json={"job_id": job_id, "expected_version": expected_version},
    )


def publish_batch(client, dataset_id: int, vcf_name: str, expected_version: int,
                  fasta_name: str = "chr1.fa") -> dict:
    job = upload(client, dataset_id, vcf_name, fasta_name).json()
    assert job["status"] == "ready", job
    r = publish(client, dataset_id, job["id"], expected_version)
    assert r.status_code == 200, r.text
    return r.json()


def make_report(client, a: int, b: int):
    return client.post("/reports", json={"dataset_a_id": a, "dataset_b_id": b})


def item_map(report):
    """报告明细 -> {(pos, ref, alt): (lines_a, lines_b)}。"""
    return {
        (it["pos"], it["ref"], it["alt"]): (it["lines_a"], it["lines_b"])
        for it in report["items"]
    }


def test_report_basic_and_reread_after_replacement(client):
    a = create_dataset(client, "批次1")
    b = create_dataset(client, "批次2")
    publish_batch(client, a, "batch1.vcf", 0)
    publish_batch(client, b, "batch2.vcf", 0)

    # 生成报告：共同位点 + 双方原始行号 + 所用版本
    r = make_report(client, a, b)
    assert r.status_code == 201, r.text
    report = r.json()
    assert report["dataset_a_id"] == a and report["dataset_b_id"] == b
    assert report["version_a"] == 1 and report["version_b"] == 1
    assert report["chrom"] == "chr1"
    assert report["common_count"] == 5
    assert item_map(report) == {
        (4, "G", "GTA"): ([5, 6], [5, 6, 7]),
        (4, "GTA", "G"): ([7], [8, 9]),
        (17, "A", "G"): ([8], [10]),
        (18, "T", "G"): ([9], [11]),
        (20, "G", "GA"): ([10], [12]),
    }

    # 双方已发布数据都被替换后，旧报告仍可按报告号读取且内容不变
    publish_batch(client, a, "batch2.vcf", 1)
    publish_batch(client, b, "batch1.vcf", 1)
    r = client.get(f"/reports/{report['id']}")
    assert r.status_code == 200
    assert r.json() == report

    # 新报告反映新所用版本（双方行号互换）
    r = make_report(client, a, b)
    assert r.status_code == 201, r.text
    report2 = r.json()
    assert report2["version_a"] == 2 and report2["version_b"] == 2
    assert item_map(report2)[(4, "G", "GTA")] == ([5, 6, 7], [5, 6])
    assert item_map(report2)[(4, "GTA", "G")] == ([8, 9], [7])

    # 列表包含两份报告摘要（不含明细）
    reports = client.get("/reports").json()
    assert [x["id"] for x in reports] == [report["id"], report2["id"]]
    assert "items" not in reports[0]
    assert reports[0]["common_count"] == 5


def test_report_rejections(client):
    a = create_dataset(client, "批次1")
    publish_batch(client, a, "batch1.vcf", 0)
    unpublished = create_dataset(client, "未发布")

    # 数据集不存在 -> 404
    assert make_report(client, a, 999).status_code == 404
    assert make_report(client, 999, a).status_code == 404

    # 尚未发布 -> 409
    r = make_report(client, a, unpublished)
    assert r.status_code == 409
    assert "尚未发布" in r.json()["detail"]
    r = make_report(client, unpublished, a)
    assert r.status_code == 409

    # 染色体不同 -> 422
    c = create_dataset(client, "chr2批次")
    publish_batch(client, c, "batch3_chr2.vcf", 0, fasta_name="chr2.fa")
    r = make_report(client, a, c)
    assert r.status_code == 422
    assert "染色体" in r.json()["detail"]

    # 报告号不存在 -> 404；以上拒绝均未留下任何报告
    assert client.get("/reports/999").status_code == 404
    assert client.get("/reports").json() == []


def test_report_empty_intersection(client):
    """染色体相同但无共同位点：正常生成，common_count=0。"""
    a = create_dataset(client, "批次1")
    publish_batch(client, a, "batch1.vcf", 0)
    d = create_dataset(client, "无交集批次")
    vcf = (
        b"##fileformat=VCFv4.2\n"
        b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        b"chr1\t25\t.\tA\tG\t50\tPASS\t.\n"
    )
    r = client.post(
        f"/datasets/{d}/uploads",
        files={
            "fasta": ("chr1.fa", (SAMPLES / "chr1.fa").read_bytes(), "text/plain"),
            "vcf": ("only.vcf", vcf, "text/plain"),
        },
    )
    job = r.json()
    assert job["status"] == "ready", job
    assert publish(client, d, job["id"], 0).status_code == 200

    r = make_report(client, a, d)
    assert r.status_code == 201, r.text
    assert r.json()["common_count"] == 0
    assert r.json()["items"] == []


def test_report_save_failure_leaves_no_partial_report(client, monkeypatch):
    """保存失败（模拟写入异常）整体回滚，不留下半份报告。"""
    a = create_dataset(client, "批次1")
    b = create_dataset(client, "批次2")
    publish_batch(client, a, "batch1.vcf", 0)
    publish_batch(client, b, "batch2.vcf", 0)

    from app.store import Store

    def boom(*args, **kwargs):
        raise RuntimeError("模拟保存失败")

    monkeypatch.setattr(Store, "_insert_report", boom)

    with pytest.raises(Exception):
        make_report(client, a, b)

    # 事务回滚：报告头与报告行都没有残留
    assert client.get("/reports").json() == []

    # 恢复后功能正常
    monkeypatch.undo()
    assert make_report(client, a, b).status_code == 201


def test_report_concurrent_replacement_snapshot_consistent(client):
    """报告生成期间并发替换已发布数据：每份报告必须来自同一时点的一致快照。

    数据集 A 被反复重新发布（batch1/batch2 交替，行号特征不同），
    并发生成的每份报告中 version_a 与 lines_a 必须属于同一个已发布版本，
    不允许出现版本号与位点内容跨版本撕裂。
    """
    a = create_dataset(client, "批次A")
    b = create_dataset(client, "批次B")
    publish_batch(client, a, "batch1.vcf", 0)  # v1 = batch1
    publish_batch(client, b, "batch2.vcf", 0)  # B 固定 v1 = batch2

    ROUNDS = 6
    barrier = threading.Barrier(4)
    errors = []
    reports = []

    def republish():
        # v2=batch2, v3=batch1, ... 交替替换 A 的已发布数据
        barrier.wait(timeout=10)
        for i in range(ROUNDS):
            vcf = "batch2.vcf" if i % 2 == 0 else "batch1.vcf"
            version = client.get(f"/datasets/{a}").json()["version"]
            job = upload(client, a, vcf).json()
            r = publish(client, a, job["id"], version)
            assert r.status_code == 200, r.text

    def generate():
        barrier.wait(timeout=10)
        for _ in range(8):
            r = make_report(client, a, b)
            assert r.status_code == 201, r.text
            reports.append(r.json())

    def run(fn):
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            errors.append(exc)

    threads = [threading.Thread(target=run, args=(republish,))] + [
        threading.Thread(target=run, args=(generate,)) for _ in range(3)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join(timeout=120)

    assert not errors, errors
    assert len(reports) == 24
    for rep in reports:
        # version_a 奇数 = batch1 内容，偶数 = batch2 内容；行号特征必须与版本号同源
        expected = BATCH1_LINES if rep["version_a"] % 2 == 1 else BATCH2_LINES
        got = {(p, r_, al): la for (p, r_, al), (la, _lb) in item_map(rep).items()}
        assert got == expected, f"报告 {rep['id']} 快照撕裂：version_a={rep['version_a']}"
        assert rep["version_b"] == 1
        for (p, r_, al), (_la, lb) in item_map(rep).items():
            assert lb == BATCH2_LINES[(p, r_, al)]

    # 并发结束后旧报告全部仍可按号读取
    for rep in reports:
        assert client.get(f"/reports/{rep['id']}").json() == rep
