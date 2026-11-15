"""变异归一化核心算法。

流程（与 vt normalize / bcftools norm 的左对齐语义一致）：
1. 用上传的单染色体 FASTA 逐条核对 VCF 行的 REF；
2. 消去 REF/ALT 的共有后缀、共有前缀；
3. 纯插入/缺失保留一个锚点碱基（优先左锚点，染色体左端用右锚点）；
4. 串联重复中的插入/缺失向左对齐（左most 表示）；
5. 仅接受 SNV / 纯插入 / 纯缺失，其余（MNP、复杂替换等）判为坏行。
"""
from __future__ import annotations


class NormalizeError(ValueError):
    """单行 VCF 记录的校验/归一化错误（坏行）。"""


def normalize(seq: str, pos: int, ref: str, alt: str) -> tuple[int, str, str]:
    """归一化一条双等位变异。

    参数:
        seq: 参考序列（0-based 字符串，已大写）。
        pos: VCF 的 1-based POS。
        ref/alt: VCF 原始 REF/ALT（已大写，且通过字符集校验）。

    返回:
        (pos, ref, alt)：归一化后的 1-based 坐标与等位基因。

    异常:
        NormalizeError: REF 不符、越界、MNP/复杂变异等非 SNV/纯插入/纯缺失情形。
    """
    if ref == alt:
        raise NormalizeError("REF 与 ALT 相同，不是变异")

    pos0 = pos - 1
    if pos0 < 0:
        raise NormalizeError(f"POS 必须 >= 1（实际 {pos}）")

    # 1) 逐条核对 REF 与参考序列
    end = pos0 + len(ref)
    if end > len(seq):
        raise NormalizeError(f"REF 超出染色体末端（染色体长 {len(seq)}，POS+len(REF)-1={end}）")
    actual = seq[pos0:end]
    if actual != ref:
        raise NormalizeError(f"REF 与参考不符：参考为 {actual}，VCF 写法为 {ref}")

    # 2) 消去共有后缀，再消去共有前缀
    while ref and alt and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    while ref and alt and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos0 += 1

    if not ref and not alt:
        raise NormalizeError("REF 与 ALT 等价（共有碱基修剪后为空）")

    # 3) 分类：只允许 SNV / 纯插入 / 纯缺失
    if ref and alt:
        if len(ref) == len(alt):
            if len(ref) == 1:
                return pos0 + 1, ref, alt  # SNV，无需对齐
            raise NormalizeError(f"不支持 MNP（{ref}>{alt}），仅支持 SNV/插入/缺失")
        raise NormalizeError(f"复杂变异（{ref}>{alt}，替换伴随长度变化），仅支持纯 SNV/插入/缺失")

    # 4) 纯插入/缺失：保留一个锚点碱基
    if pos0 > 0:
        # 左锚点：取事件前一个参考碱基
        pos0 -= 1
        base = seq[pos0]
        ref = base + ref
        alt = base + alt
    else:
        # 染色体最左端无法左锚，改用右锚点
        idx = pos0 + len(ref)  # 插入时 ref 为空 -> 取 seq[0]
        if idx >= len(seq):
            raise NormalizeError("变异位于染色体末端，无法保留锚点碱基")
        base = seq[idx]
        ref = ref + base
        alt = alt + base

    # 5) 左对齐：串联重复中，若锚点碱基与差异序列末位相同，则整体左移一位并旋转差异序列
    if len(alt) > len(ref):
        diff = alt[len(ref):]  # 插入序列（此时 len(ref)==1）
        while pos0 >= 1 and seq[pos0] == diff[-1]:
            diff = diff[-1] + diff[:-1]
            pos0 -= 1
        ref = seq[pos0]
        alt = seq[pos0] + diff
    else:
        diff = ref[len(alt):]  # 缺失序列（此时 len(alt)==1）
        while pos0 >= 1 and seq[pos0] == diff[-1]:
            diff = diff[-1] + diff[:-1]
            pos0 -= 1
        ref = seq[pos0] + diff
        alt = seq[pos0]

    return pos0 + 1, ref, alt
