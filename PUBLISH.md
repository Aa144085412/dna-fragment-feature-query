# 维护者发布说明

仓库：`https://github.com/Aa144085412/dna-fragment-feature-query`，项目名 **DNA Fragment Feature Query**。

## 发布前检查

- 本地运行 `python -m unittest discover -s tests -v`；单元测试不能代替真实端到端验收。
- `git status --short` 和 `git diff --cached --stat` 不得包含实验 `.dna`、公司 `.xlsx`、测序峰图、比对索引、凭据、运行输出。
- 核对 README、SKILL.md、CLI 实际参数一致。
- 对正式科研用途，另跑真实 `SnapGene/Excel → hg38 BLAST → UCSC API → review.xlsx/TSV/JSON/MD`，记录是否确实通过。

## Git 命令（适用于现有仓库更新）

```bash
git add SKILL.md README.md PUBLISH.md .gitignore scripts/dfq.py tests/test_core.py
git diff --cached --stat
git commit -m "feat: release DNA Fragment Feature Query agent skill"
git push origin main
```

此说明不作为推送成功证据；请以实际 GitHub 提交和远端文件为准。当前仓库的可见性由 GitHub 设置决定；**公共仓库绝不可上传原始实验数据**。
