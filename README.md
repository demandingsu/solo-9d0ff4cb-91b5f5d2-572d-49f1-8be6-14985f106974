# VCF 变异归一化服务

面向测序中心数据员的跨批次变异归一化 API：不同批次 VCF 对同一变异常有不同写法
（锚点位置、共有碱基、重复区内的坐标漂移），本服务将其归一化为唯一表示并归并等价位点，
供跨批次比对。FastAPI + SQLite，Docker Compose 一键启动。

## 功能与规则

**输入**：单染色体 FASTA + 仅含双等位 SNV / 插入 / 缺失的 VCF。

**归一化流水线**（`app/normalize.py`，与 vt normalize / bcftools norm 左对齐语义一致）：

1. **逐条核对 REF**：VCF 每行的 REF 必须与上传 FASTA 对应位置完全一致（大小写不敏感）；
2. **消去共有碱基**：先修剪 REF/ALT 共有后缀，再修剪共有前缀；
3. **保留合法锚点**：纯插入/缺失保留一个锚点碱基（优先左锚点，染色体左端用右锚点）；
4. **左对齐**：串联重复中的插入/缺失左移至最左等价表示；
5. **归并等价位点**：归一化后 (POS, REF, ALT) 相同的行合并，保留**各原始行号**
   （VCF 文件物理行号，从 1 起，含头部行）。

**坏行**（任一即判坏行，含坏行的任务不可发布）：REF 与参考不符、CHROM 与 FASTA 不一致、
POS 非整数/越界、REF 超出染色体末端、多等位（ALT 含逗号）、符号/星号 ALT、MNP、
替换伴随长度变化的复杂变异、REF==ALT、列数不足。

**文件级致命错误**（HTTP 400，不产生任务）：FASTA 含多条序列、非法碱基字符、
VCF 缺 `#CHROM` 头、非 UTF-8 文件。

## 一致性模型

- 上传生成**暂存任务**（job），处理期间查询仍返回**完整旧版**；
- 发布按 `数据集号 + expected_version`（查询所得版本）在**单事务内整体替换**；
- **并发替换只允许版本匹配者成功**（乐观锁）：不匹配返回 409 及当前版本；
- **坏行任务或发布失败不改动旧版**：失败任务发布返回 422 及逐行原因，任何失败均回滚。

## 快速开始

```bash
docker compose up --build
```

- API 文档（Swagger UI）：http://localhost:8000/docs
- 健康检查：http://localhost:8000/healthz
- SQLite 数据持久化在 named volume `variant-data`（容器内 `/data/variants.db`）

### 一键演示（等价变异归并 + 版本乐观锁 + 坏行隔离）

```bash
python3 scripts/demo.py            # 默认 http://localhost:8000
```

### curl 示例

```bash
# 1. 创建数据集（version 从 0 起）
curl -s -X POST localhost:8000/datasets -H 'Content-Type: application/json' \
     -d '{"name":"my-dataset"}'

# 2. 上传 FASTA + VCF，生成暂存任务（不影响已发布数据）
curl -s -X POST localhost:8000/datasets/1/uploads \
     -F 'fasta=@samples/chr1.fa' -F 'vcf=@samples/batch1.vcf'

# 3. 查询当前版本（处理期间/发布前仍返回旧版）
curl -s localhost:8000/datasets/1

# 4. 按查询所得版本发布（乐观锁；不匹配则 409，坏行任务则 422）
curl -s -X POST localhost:8000/datasets/1/publish \
     -H 'Content-Type: application/json' -d '{"job_id":1,"expected_version":0}'

# 5. 查询任务状态（坏行明细 / 暂存结果）
curl -s localhost:8000/datasets/1/jobs/1
```

## API 一览

| 方法 | 路径 | 说明 |
|---|---|---|
| POST | `/datasets` | 创建数据集，返回 `version=0` |
| GET | `/datasets` | 列出全部数据集及版本 |
| GET | `/datasets/{id}` | 查询当前已发布版本（坐标、REF/ALT、各原始行号） |
| POST | `/datasets/{id}/uploads` | 上传 FASTA+VCF 并归一化，生成暂存任务 |
| GET | `/datasets/{id}/jobs/{job_id}` | 任务状态、坏行明细、暂存的归一化结果 |
| POST | `/datasets/{id}/publish` | 按 `expected_version` 乐观锁发布；409=版本冲突，422=任务含坏行 |
| GET | `/healthz` | 健康检查 |

## 演示样本（`samples/`）

`chr1.fa`：40 bp 演示染色体，5–12 位为 TA 串联重复（`GCAG TATATATA GCCG...`）。

`batch1.vcf` 与 `batch2.vcf` 是**同一批变异在两批数据中的不同写法**：

| 变异 | 批次 1 写法（行号） | 批次 2 写法（行号） | 归一化结果 |
|---|---|---|---|
| TA 重复区插入 TA | `6 A>ATA`(5)、`8 A>ATA`(6) | `10 A>ATA`(5)、`12 A>ATA`(6)、`5 T>TAT`(7) | **4 G>GTA** |
| TA 重复区缺失 TA | `6 ATA>A`(7) | `4 GTA>G`(8)、`5 TAT>T`(9) | **4 GTA>G** |
| SNV | `17 A>G`(8) | `17 A>G`(10) | **17 A>G** |
| SNV（带共有碱基写法） | `18 TT>GT`(9) | `18 T>G`(11) | **18 T>G** |
| 非重复区插入 | `20 G>GA`(10) | `20 G>GA`(12) | **20 G>GA** |

两批发布后归一化结果**完全一致**（仅原始行号不同），可直接跨批次比对。

`bad.vcf`：5 条坏行（REF 不符、多等位、MNP、越界、CHROM 不符），用于演示坏行隔离。

## 本地开发与测试

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/python -m pytest tests/ -q          # 27 个测试：算法 + API + 并发发布
DATABASE_PATH=data/variants.db .venv/bin/python -m uvicorn app.main:app --reload
```

## 目录结构

```
app/
  main.py        FastAPI 入口与路由
  normalize.py   归一化算法（核对/修剪/锚点/左对齐）
  parse.py       FASTA/VCF 结构解析
  service.py     逐行校验、归并等价位点
  store.py       SQLite 存储与乐观锁发布事务
  schemas.py     请求/响应模型
samples/         演示用 FASTA 与三份 VCF
scripts/demo.py  端到端演示脚本（仅标准库）
tests/           pytest 测试
Dockerfile / docker-compose.yml
```
