# Hulab Codex Skills

这是 Hulab 实验室整理的 Codex skills 集合，供团队成员和其他研究者直接安装使用。

## 包含的 skills

| 目录 | 用途 |
| --- | --- |
| `genome-annotation/` | 细菌基因组 Bakta 注释 |
| `genome-download/` | 根据物种名单或 GCF 名单批量下载 NCBI RefSeq 基因组文件 |
| `genome-search/` | 检索并审计 NCBI GCF/GCA 基因组名单 |
| `r-deseq2/` | R DESeq2 差异表达分析 |
| `r-wgcna/` | R WGCNA 共表达网络与候选基因挖掘 |
| `rnaseq-denovo/` | RNA-seq de novo 定量流程（Trinity + TransDecoder + Salmon） |
| `rnaseq-reference/` | RNA-seq reference-based 定量流程（HISAT2 + StringTie） |
| `structure-mapping/` | Foldseek / USalign / DEER 结构比对与检索 |
| `structure-prediction/` | PDB / AlphaFold DB 检索、同源建模与结构预测 |

所有 skill 都是可移植的：不依赖固定服务器路径，用户可以使用自己的 conda 环境、R 环境，以及自己的软件和数据库路径。具体安装和运行说明见各 skill 目录。

## 安装方法

### 方法一：直接复制到 Codex skills 目录

```bash
git clone https://github.com/Hulab-ECG/Skill.git
mkdir -p ~/.codex/skills
cp -r Skill/* ~/.codex/skills/
```

完成后，重新打开或进入下一个 Codex 会话即可使用。

### 方法二：只安装其中某个 skill

```bash
cp -r Skill/genome-download ~/.codex/skills/
```

### 方法三：通过 Codex 内的 skill-installer

如果目标机器已有 `skill-installer`，也可以让 Codex 从仓库自动安装单个 skill：

```bash
python3 ~/.codex/skills/.system/skill-installer/scripts/install-skill-from-github.py \
  --repo Hulab-ECG/Skill \
  --path genome-download
```

## 使用前需要注意

- 这些 skill 只负责提供流程和脚本说明，实际运行前请检查目标机器上的 conda/R/Python 环境。
- 每个 skill 中出现的 `PRJNA...`、`SRR...`、`GCA_...`、`/path/to/...` 等仅为示例，不是固定默认值。
- 涉及大型数据库下载、基因组下载、结构数据库构建或 GPU 运算时，请先确认服务器资源和权限。
- 使用 MODELLER 时，请遵守其许可协议并配置好自己的许可证。

## 更新

要更新已有安装，重新 `git pull` 仓库后重新复制，或替换对应 skill 目录即可：

```bash
git -C Skill pull
cp -r Skill/* ~/.codex/skills/
```
