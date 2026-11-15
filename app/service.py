"""上传内容处理：逐条校验、归一化、归并等价位点。"""
from __future__ import annotations

from dataclasses import dataclass

from .normalize import NormalizeError, normalize
from .parse import VcfRow, parse_fasta, parse_vcf_rows

# VCF REF/ALT 允许的碱基字符
VCF_BASES = frozenset("ACGTN")


@dataclass
class RowError:
    """一条坏行的诊断信息。"""

    line: int
    reason: str
    chrom: str | None = None
    pos: str | None = None
    ref: str | None = None
    alt: str | None = None

    def to_dict(self) -> dict:
        return {
            "line": self.line,
            "reason": self.reason,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
        }


@dataclass
class NormVariant:
    """归一化后的位点，lines 为归并到该位点的所有原始行号。"""

    pos: int
    ref: str
    alt: str
    lines: list[int]


def _check_row(fasta_chrom: str, row: VcfRow) -> tuple[int, str, str]:
    """单行结构校验，返回 (pos, ref, alt)；不合法抛 NormalizeError。"""
    if row.parse_error:
        raise NormalizeError(row.parse_error)
    if row.chrom != fasta_chrom:
        raise NormalizeError(f"CHROM={row.chrom} 与 FASTA 染色体 {fasta_chrom} 不一致")
    try:
        pos = int(row.pos)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        raise NormalizeError(f"POS 不是整数：{row.pos!r}")
    ref = (row.ref or "").upper()
    alt = (row.alt or "").upper()
    if not ref or any(b not in VCF_BASES for b in ref):
        raise NormalizeError(f"REF 含非法字符：{row.ref!r}（仅支持 ACGTN）")
    if "," in alt:
        raise NormalizeError("ALT 含逗号：仅支持双等位变异")
    if not alt or any(b not in VCF_BASES for b in alt):
        raise NormalizeError(f"ALT 非法：{row.alt!r}（仅支持 ACGTN，不支持 <符号>/断点/*/.）")
    return pos, ref, alt


def process_upload(fasta_text: str, vcf_text: str) -> tuple[str, int, list[NormVariant], list[RowError]]:
    """处理一次上传。

    返回:
        (chrom, total_records, variants, errors)
        variants 已按 (pos, ref, alt) 归并等价位点并排序，lines 为升序原始行号。

    异常:
        ParseError: 文件级致命错误（由调用方转为 HTTP 400）。
    """
    chrom, seq = parse_fasta(fasta_text)
    rows = parse_vcf_rows(vcf_text)

    merged: dict[tuple[int, str, str], list[int]] = {}
    errors: list[RowError] = []
    for row in rows:
        try:
            pos, ref, alt = _check_row(chrom, row)
            npos, nref, nalt = normalize(seq, pos, ref, alt)
        except NormalizeError as exc:
            errors.append(
                RowError(line=row.line, reason=str(exc), chrom=row.chrom, pos=row.pos, ref=row.ref, alt=row.alt)
            )
            continue
        merged.setdefault((npos, nref, nalt), []).append(row.line)

    variants = [
        NormVariant(pos=p, ref=r, alt=a, lines=sorted(lines))
        for (p, r, a), lines in sorted(merged.items())
    ]
    return chrom, len(rows), variants, errors
