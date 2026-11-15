"""FASTA / VCF 文本的结构性解析。

文件级致命问题（FASTA 多条序列、VCF 缺 #CHROM 头等）抛 ParseError，
对应 HTTP 400；单行问题由调用方记为坏行。
"""
from __future__ import annotations

from dataclasses import dataclass

# 参考序列允许的 IUPAC 碱基符号
IUPAC_BASES = frozenset("ACGTRYSWKMBDHVN")


class ParseError(ValueError):
    """文件级致命解析错误（整个上传被拒绝，不产生任务）。"""


def parse_fasta(text: str) -> tuple[str, str]:
    """解析单染色体 FASTA，返回 (染色体名, 大写序列)。"""
    name: str | None = None
    chunks: list[str] = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name is not None:
                raise ParseError(f"第 {lineno} 行：检测到第二条序列，本服务要求单染色体 FASTA")
            header = line[1:].strip()
            if not header:
                raise ParseError(f"第 {lineno} 行：FASTA 头缺少序列名")
            name = header.split()[0]
            continue
        if name is None:
            raise ParseError(f"第 {lineno} 行：序列出现在 FASTA 头之前")
        seq_line = line.upper().replace(" ", "")
        invalid = sorted(set(seq_line) - IUPAC_BASES)
        if invalid:
            raise ParseError(f"第 {lineno} 行：非法碱基字符 {invalid}")
        chunks.append(seq_line)
    if name is None:
        raise ParseError("缺少 FASTA 头（>染色体名）")
    seq = "".join(chunks)
    if not seq:
        raise ParseError("FASTA 序列为空")
    return name, seq


@dataclass
class VcfRow:
    """一条 VCF 数据行。line 为文件内物理行号（从 1 起，含头部行）。"""

    line: int
    chrom: str | None
    pos: str | None
    ref: str | None
    alt: str | None
    parse_error: str | None = None


def parse_vcf_rows(text: str) -> list[VcfRow]:
    """解析 VCF，返回全部数据行（含行号）。

    缺少 #CHROM 头行等结构问题抛 ParseError；单列数不足记为坏行返回。
    """
    rows: list[VcfRow] = []
    header_seen = False
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith("##"):
            continue
        if line.startswith("#"):
            if header_seen:
                raise ParseError(f"第 {lineno} 行：重复的 #CHROM 头行")
            if not line.upper().startswith("#CHROM"):
                raise ParseError(f"第 {lineno} 行：头部行应以 #CHROM 开头")
            if len(line.split("\t")) < 8:
                raise ParseError(f"第 {lineno} 行：VCF 头列数不足 8")
            header_seen = True
            continue
        if not header_seen:
            raise ParseError(f"第 {lineno} 行：数据行出现在 #CHROM 头行之前")
        fields = line.split("\t")
        if len(fields) < 5:
            rows.append(
                VcfRow(
                    line=lineno,
                    chrom=None,
                    pos=None,
                    ref=None,
                    alt=None,
                    parse_error=f"列数不足（{len(fields)} 列，至少需要 CHROM/POS/ID/REF/ALT 5 列，制表符分隔）",
                )
            )
            continue
        rows.append(VcfRow(line=lineno, chrom=fields[0], pos=fields[1], ref=fields[3], alt=fields[4]))
    if not header_seen:
        raise ParseError("VCF 缺少 #CHROM 头行")
    return rows
