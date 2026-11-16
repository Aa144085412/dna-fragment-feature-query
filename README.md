# DNA Fragment Feature Query

**DNA Fragment Feature Query** 是一个可安装的 **Agent Skill**：面向合成 DNA 片段的突变评估，按「SnapGene/Excel 注释 → 本地 hg38 全基因组 BLAST → UCSC 参考碱基与功能注释 → 逐条修复优先级 → 审计报告」的顺序执行。

> **用途边界：** 这是生物信息学证据整理工具，不是突变真实性检测、DNA 自动修复、临床变异解释或实验有效性证明。`.dna` 序列可能是参考序列，报告突变可能仅出现在注释中；插件会保留这个区别。

## 能做什么

- 提取多个 SnapGene `.dna` 文件中的 `misc_feature` 突变标签、载体位置和周围约 150 bp DNA。
- 可选解析供应商 Excel：**B 列片段编号**、**E 列突变清单**（识别诸如 `1-10`、`ProjectABC-2-8` 形式的片段号），并提示 SnapGene/Excel 不一致。
- 使用**既有本地 NCBI BLAST+ hg38 索引**定位到染色体坐标，保留次优命中、覆盖率、方向、可能的多重匹配或 alt contig。
- 通过 UCSC 官方 REST API 查询 `getData/sequence`、`rmsk`（RepeatMasker）、`cCREregistry`（ENCODE4 候选调控元件）、`ncbiRefSeq`（转录本与外显子）。
- 对合适的 RefSeq 编码区 SNV 做**预测**密码子/氨基酸变化，识别可能的非同义替换。
- 输出包含**判定证据与未确认事项**的修复建议，不把“没有注释”当作“没有功能”。

## 安装

### 前置依赖

1. Python **3.10+**（本插件运行代码仅使用 Python 标准库，无需 `pip install`）。
2. [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)，找到 `blastn` 可执行文件。
3. 已经用 `makeblastdb -dbtype nucl` 建索引的 **hg38** 数据库，找到数据库前缀（例如路径 `.../blast/hg38` 对应 `hg38.nin/.nhr/.nsq`）。可以直接复用 **PrimerForge / SeqAtlas** 已有数据库，不要重复建库。
4. 只有启用 `--online` 才需要访问 `https://api.genome.ucsc.edu`。

### 安装为 Agent Skill（立即可用）

解压本项目发布的 `DNA-Fragment-Feature-Query-v0.1.0.zip`，得到 `dna-fragment-feature-query/` 文件夹。无需新建服务，直接复制整个目录到下列 Skills 位置。

**从 GitHub 安装：**

```bash
git clone https://github.com/Aa144085412/dna-fragment-feature-query.git
```

将克隆目录（保持名称 **`dna-fragment-feature-query`**）复制或克隆到所用 Agent 的 Skills 目录：

| Agent | 用户级安装目录 | 项目级目录 |
|---|---|---|
| Claude Code | `~/.claude/skills/dna-fragment-feature-query/` | `<project>/.claude/skills/dna-fragment-feature-query/` |
| Codex / 支持 Agent Skills 的 Agent | `~/.agents/skills/dna-fragment-feature-query/` | `<project>/.agents/skills/dna-fragment-feature-query/` |

目录根部必须存在 `SKILL.md`。在支持 Skill 的 Agent 中输入：

> 使用 DNA Fragment Feature Query 审查这个文件夹里的 SnapGene `.dna` 文件，使用现有 hg38 BLAST 数据库；如果批准访问 UCSC，请对所有变异出具 Excel 和 Markdown 逐条报告。

其他 Agent（如 OpenMausBot）可以在其已有的 Skill 注册配置中登记本仓库的 `SKILL.md`，**不要为了使用本工具新建 Agent 或授权网关**；具体入口取决于原系统的 Skill 发现机制。

### 不通过 Agent，直接运行 CLI

**Windows PowerShell 示例（复用 PrimerForge 的本地 BLAST）：**

```powershell
python scripts/dfq.py `
  --input "C:\path\to\SnapGene-files" `
  --xlsx "C:\path\to\mutation-manifest.xlsx" `
  --blastn "$env:USERPROFILE\Documents\PrimerForge-Inputs\tools\blast\bin\blastn.exe" `
  --db "$env:USERPROFILE\Documents\PrimerForge-Inputs\references\grch38\blast\hg38" `
  --online `
  --out "C:\path\to\dfq-results"
```

**Linux/macOS：**

```bash
python3 scripts/dfq.py \
  --input ./private-dna \
  --xlsx ./vendor.xlsx \
  --blastn blastn \
  --db /path/to/indexed/hg38 \
  --online \
  --out ./dfq-results
```

`--xlsx` 可省略。默认提取变异附近 DNA（`--query-min 1`）；若质粒载体前 400 bp 不是基因组插入片段，可使用 `--query-min 401` 限制窗口左端，避免将载体部分送入定位。**省略 `--online`** 则只运行本地 BLAST；由于没有查询在线功能注释，对所有成功定位位点的结论会显示**证据不足**，不会假装已核查调控元件。使用 `--online` 时会向 UCSC 查询染色体坐标，不会把完整 `.dna` 序列上传给 UCSC。

## 结果文件

| 输出 | 说明 |
|---|---|
| `review.xlsx` | 可在 Excel 打开的逐条核查表，包括位点、注释、参考碱基、建议与 UCSC 链接 |
| `review.tsv` | 带表头、适合后续处理的明细 |
| `result.json` | 输入 SHA-256、原始标签、BLAST 前十候选、UCSC 原始查询记录、核验/不确定性状态 |
| `report.md` | 修复分级统计及每条记录摘要 |

**修复结论不会自动修改任何实验设计。** `建议修复`、`优先评估修复`、`条件修复`、`需人工评估`、`核对原始注释`、`一般可暂不修复`、`暂不判定`、`证据不足` 均有显式规则。任何 SNV/indel 真实性仍需核对原始测序图谱。

## 方法和科学限制

1. **明确源数据。** SnapGene 注释、Excel 公司记录与 `.dna` 实际序列分别保留，不把报告 REF/ALT 误称为测出的变异。
2. **参考定位。** hg38 BLAST 使用 ±150 bp 窗口，默认高覆盖/高一致性筛选，保留多重命中信息。结果可能受到参考组装版本、重复序列、BLAST 参数影响。
3. **注释查询。** UCSC 官方接口返回 hg38 的 RepeatMasker、ENCODE4 cCRE 和 NCBI RefSeq；某个位点查不到 cCRE 并不证明它不是调控 DNA。
4. **编码改变计算。** 仅在 RefSeq CDS 外显子内，考虑链方向、转录本剪接结构预测密码子变化。不等于功能实验结果。
5. **重复区 indel。** A/T 多聚体、AC/AAAC 重复区在无原始 trace 时难以确认缺失/插入的具体标准化锚点，输出以条件修复/人工评估为主。
6. **拒绝伪验收。** `BLAST` 运行成功 ≠ 全局唯一定位；API 可访问 ≠ 每条完成；单元测试通过 ≠ 真实用户数据端到端通过。若 API 出错，报告明确保留失败。

此工具使用本地 NCBI BLAST+，**不自动操作 UCSC 网站的人机交互式 BLAT 页面**；本地比对与官方 UCSC API 注释组合是可重复的替代工作流。

## 测试

```bash
python -m unittest discover -s tests -v
```

单元测试不请求 UCSC、不包含私有 DNA。正式验收必须另跑一遍真实 `SnapGene/Excel → hg38 BLAST+ → UCSC API → 4 份报告`，核查所有坐标及未覆盖项。

## 目录

```text
dna-fragment-feature-query/
├── SKILL.md
├── README.md
├── scripts/dfq.py
├── tests/test_core.py
├── PUBLISH.md
└── .gitignore
```

## 数据隐私

`.gitignore` 排除真实 `.dna`、测序、Genome FASTA、BLAST 索引、Excel 和运行输出；**禁止将实验私有数据、完整公司文档或本地访问凭据提交到 GitHub**。本插件不内置新 Agent、模型、服务器、数据存储或管理员认证机制。

## 参考

- [UCSC REST API](https://genome.ucsc.edu/goldenPath/help/api.html)
- [Agent Skills Specification](https://agentskills.io/specification)
- [NCBI BLAST+ Documentation](https://blast.ncbi.nlm.nih.gov/doc/blast-help/)
