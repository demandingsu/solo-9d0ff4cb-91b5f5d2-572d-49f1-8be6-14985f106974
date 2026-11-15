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

## 批次交集报告

质控员核查两个批次的已发布变异是否重合：`POST /reports` 提交两个数据集号，
服务按**归一化后的染色体、坐标、REF/ALT** 找出共同位点，保存**双方各自的原始行号**
与**生成时点双方所用的已发布版本号**，返回报告号供后续复核（`GET /reports/{id}`）。

- **快照一致**：报告在单个 `BEGIN IMMEDIATE` 事务内读取两份已发布数据并写入——
  生成期间若有并发替换（发布同样需要写锁，被阻塞），报告仍来自**同一时点的两份完整
  已发布数据**，版本号与位点内容不会跨版本撕裂；
- **原子保存**：报告头与全部共同位点行同事务提交，保存失败整体回滚，**不留半份报告**；
- **持久复核**：报告与已发布数据解耦，替换已发布数据后旧报告仍可按报告号读取；
- **拒绝生成**：数据集不存在（404）、任一数据集尚未发布（409）、
  双方染色体不同或缺失（422）。

## 一致性模型

- 上传生成**暂存任务**（job），处理期间查询仍返回**完整旧版**；
- 发布按 `数据集号 + expected_version`（查询所得版本）在**单事务内整体替换**；
- **并发替换只允许版本匹配者成功**（乐观锁）：不匹配返回 409 及当前版本；
- **坏行任务或发布失败不改动旧版**：失败任务发布返回 422 及逐行原因，任何失败均回滚；
- **交集报告**在单事务快照内生成：并发替换不撕裂、保存失败不留半份、替换后旧报告仍可按号读取。

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

### 交集报告演示（共同位点 + 并发替换一致性 + 旧报告复核）

```bash
python3 scripts/demo_report.py     # 默认 http://localhost:8000
```

演示内容：两个批次发布后生成交集报告（共同位点、双方原始行号、所用版本）→
后台线程反复重新发布数据集 A 的同时连续生成报告（逐份校验版本号与位点内容同源、
无快照撕裂）→ 替换后按报告号复核旧报告（内容不变）→ 三类拒绝场景
（未发布 409 / 染色体不同 422 / 不存在 404）。

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

# 6. 生成两个数据集已发布变异的交集报告（返回报告号）
curl -s -X POST localhost:8000/reports -H 'Content-Type: application/json' \
     -d '{"dataset_a_id":1,"dataset_b_id":2}'

# 7. 按报告号复核（已发布数据被替换后仍可读）
curl -s localhost:8000/reports/1
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
| POST | `/reports` | 生成两数据集已发布变异的交集报告；404=数据集不存在，409=尚未发布，422=染色体不同 |
| GET | `/reports` | 列出全部交集报告（摘要，不含明细） |
| GET | `/reports/{id}` | 按报告号读取交集报告（含共同位点与双方行号；替换后仍可读） |
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

`chr2.fa` 与 `batch3_chr2.vcf`：chr2 演示批次，用于演示交集报告"染色体不同则拒绝"（422）。

## 本地开发与测试

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/python -m pytest tests/ -q          # 32 个测试：算法 + API + 并发发布 + 交集报告
DATABASE_PATH=data/variants.db .venv/bin/python -m uvicorn app.main:app --reload
```

## 目录结构

```
app/
  main.py        FastAPI 入口与路由
  normalize.py   归一化算法（核对/修剪/锚点/左对齐）
  parse.py       FASTA/VCF 结构解析
  service.py     逐行校验、归并等价位点
  store.py       SQLite 存储、乐观锁发布事务、交集报告事务
  schemas.py     请求/响应模型
samples/         演示用 FASTA 与 VCF（含 chr2 拒绝场景样本）
scripts/demo.py         端到端演示：归并/乐观锁/坏行隔离（仅标准库）
scripts/demo_report.py  交集报告演示：共同位点/并发替换/旧报告复核（仅标准库）
tests/           pytest 测试
Dockerfile / docker-compose.yml
```
