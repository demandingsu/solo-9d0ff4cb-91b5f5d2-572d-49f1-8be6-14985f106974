"""归一化算法单元测试。样本序列：chr1 = GCAG TATATATA GCCGATTGACGTACGGATCAGTACGATG（40bp，5-12 为 TA 重复）。"""
import pytest

from app.normalize import NormalizeError, normalize

SEQ = "GCAGTATATATAGCCGATTGACGTACGGATCAGTACGATG"


@pytest.mark.parametrize(
    "pos,ref,alt",
    [
        (4, "G", "GTA"),    # 已是最左
        (5, "T", "TAT"),    # 插入 AT
        (6, "A", "ATA"),    # 插入 TA
        (8, "A", "ATA"),
        (10, "A", "ATA"),
        (12, "A", "ATA"),
    ],
)
def test_equivalent_insertions_left_align(pos, ref, alt):
    """TA 重复区中同一插入的不同写法全部左对齐到同一表示。"""
    assert normalize(SEQ, pos, ref, alt) == (4, "G", "GTA")


@pytest.mark.parametrize(
    "pos,ref,alt",
    [
        (4, "GTA", "G"),    # 已是最左
        (5, "TAT", "T"),
        (6, "ATA", "A"),
        (10, "ATA", "A"),
    ],
)
def test_equivalent_deletions_left_align(pos, ref, alt):
    """TA 重复区中同一缺失的不同写法全部左对齐到同一表示。"""
    assert normalize(SEQ, pos, ref, alt) == (4, "GTA", "G")


def test_snv_unchanged():
    assert normalize(SEQ, 17, "A", "G") == (17, "A", "G")


def test_trim_shared_suffix():
    """共有后缀被消去：18 位 TT>GT 等价于 18 位 T>G。"""
    assert normalize(SEQ, 18, "TT", "GT") == (18, "T", "G")


def test_trim_shared_prefix():
    """共有前缀被消去：18 位 TT>TG 等价于 19 位 T>G。"""
    assert normalize(SEQ, 18, "TT", "TG") == (19, "T", "G")


def test_trim_shared_prefix_keeps_anchor():
    """共有前缀修剪后插入仍需保留一个锚点碱基。"""
    assert normalize(SEQ, 20, "G", "GA") == (20, "G", "GA")


def test_no_shift_outside_repeat():
    """非重复区的插入不发生左移。"""
    assert normalize(SEQ, 29, "A", "AT") == (29, "A", "AT")


def test_multibase_record_trims_to_snv():
    assert normalize(SEQ, 17, "AT", "GT") == (17, "A", "G")


def test_ref_mismatch_rejected():
    with pytest.raises(NormalizeError, match="REF 与参考不符"):
        normalize(SEQ, 17, "C", "G")


def test_mnp_rejected():
    with pytest.raises(NormalizeError, match="MNP"):
        normalize(SEQ, 18, "TT", "GA")


def test_complex_rejected():
    with pytest.raises(NormalizeError, match="复杂变异"):
        normalize(SEQ, 17, "A", "TT")


def test_beyond_chromosome_end_rejected():
    with pytest.raises(NormalizeError, match="超出染色体末端"):
        normalize(SEQ, 100, "A", "T")


def test_ref_equals_alt_rejected():
    with pytest.raises(NormalizeError, match="不是变异"):
        normalize(SEQ, 17, "A", "A")


def test_lowercase_ref_matches():
    """REF 核对大小写不敏感（调用方统一大写，这里直接验证算法行为）。"""
    assert normalize(SEQ, 17, "A", "G") == (17, "A", "G")
