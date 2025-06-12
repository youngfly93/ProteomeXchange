1. 环境准备
bash
复制
编辑
# conda 环境
conda create -n pxd_annot python=3.11 \
             pandas ppx requests pyyaml regex tqdm joblib
conda activate pxd_annot
ppx 用于高层封装 PRIDE / MassIVE API，带本地缓存功能 
ppx.readthedocs.io
。

joblib 的 Memory 可把 API JSON 缓存到磁盘，重复运行几乎零延迟 
joblib.readthedocs.io
。

2. 输入与输出
文件	说明
ids.txt	一行一个 accession（PXD/MSV/PASS…），置于版本控制。
diseases.yml	200+ 条疾病正则（如 `Lung_Cancer: "(lung cancer
输出 dataset_annotation.tsv	四列：all_accession, HLA, Scenario, Disease。
输出 needs_manual.csv	无法自动判定或冲突条目。

3. 拉取元数据
python
复制
编辑
import ppx, requests, itertools, joblib, yaml

mem = joblib.Memory(".cache", verbose=0)

@mem.cache
def fetch_json(acc):
    base = "https://www.ebi.ac.uk/pride/ws/archive/v2/projects"
    proj = requests.get(f"{base}/{acc}", timeout=15).json()         # /projects
    samp = requests.get(f"{base}/{acc}/samples", timeout=15).json() # /samples
    return proj, samp

def build_text(proj, samp):
    fields = [
        proj.get("projectTitle",""),
        proj.get("projectDescription",""),
        " ".join(proj.get("keywords", [])),
        " ".join(a.get("value","") for a in proj.get("additionalAttributes",[])),
        " ".join(s.get("attributes","") for s in samp or [])
    ]
    return " ".join(fields).lower()
/projects/{acc}、/projects/{acc}/samples 是官方文档给出的两个关键端点，可一次拿到 title/description/keywords 等核心字段 
ebi.ac.uk
ebi.ac.uk
。

additionalAttributes 常直接写明 “HLA Class I” 或疾病名称 
pmc.ncbi.nlm.nih.gov
。

4. 关键词库（外置 YAML/JSON）
4.1 HLA 正则
yaml
复制
编辑
HLA_I: "\\b(hla[- ]?i|class[- ]?i|mhc[- ]?i)\\b"
HLA_II: "\\b(hla[- ]?ii|class[- ]?ii|mhc[- ]?ii)\\b"
4.2 场景关键词
存放在 scenarios.yml，五大类（Cancer / Infection / Autoimmune / Normal / Mixed）；每类给出 5-10 个正则。关键词来自常见文献描述 �� cite turn0search2 turn0search9〉。

4.3 疾病正则
维护在 diseases.yml；示例：

yaml
复制
编辑
Lung_Cancer: "(lung cancer|pulmonary carcinoma|nsclc|sclc)"
Type_1_Diabetes: "(type ?1 diabetes|t1d)"
COVID-19: "(covid-19|sars[- ]cov[- ]2|coronavirus)"
易于社区协作、Pull Request 更新。

5. 分类函数
python
复制
编辑
import re

def classify(text, hla_rx, scene_dict, disease_rx):
    # 1) HLA
    has_I  = bool(re.search(hla_rx["HLA_I"],  text))
    has_II = bool(re.search(hla_rx["HLA_II"], text))
    hla = "I/II" if has_I and has_II else "I" if has_I else "II" if has_II else "Unspecified"

    # 2) 场景
    matches = [k for k,v in scene_dict.items() if re.search(v, text)]
    scenario = "Mixed" if len(matches) > 1 else matches[0] if matches else "Unspecified"

    # 3) 疾病
    for label, rx in disease_rx.items():
        if re.search(rx, text):
            return hla, scenario, label
    return hla, scenario, "Unspecified"
Python re 模块是官方推荐的正则实现 
w3schools.com
docs.python.org
。

6. 歧义与手工复核
python
复制
编辑
if hla == "Unspecified" or scenario == "Unspecified" or disease == "Unspecified":
    ambig_writer.writerow([acc, hla, scenario, disease])
needs_manual.csv 让生信同事只花 2–3 分钟浏览真正棘手的条目（通常 <5%）。

7. 性能：缓存 + 并行
joblib.Memory 把 JSON 保存在 .cache/，后续 run 秒级完成 
joblib.readthedocs.io
。

大清单可用 concurrent.futures.ThreadPoolExecutor ⬆︎ 并发拉取，PRIDE API 支持 20-40 QPS。

8. 结果写出
python
复制
编辑
df.to_csv("dataset_annotation.tsv", sep="\t", index=False)
与现有 anchor_mod/Snakemake 流水线对齐，后续规则直接 pandas.read_csv 过滤条件：df[df.HLA=="I"] 等。

9. GitHub Actions 自动化
yaml
复制
编辑
# .github/workflows/annot.yml
name: annotate
on:
  push:
    paths: [ids.txt, diseases.yml, scripts/**]
jobs:
  run:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: mamba-org/setup-micromamba@v1
        with:
          environment-file: env.yml
      - run: python scripts/annotate.py
      - run: git config user.email "bot@ci" && git config user.name "CI Bot"
      - run: git add dataset_annotation.tsv && git commit -m "auto: update annotations" || echo "no change"
      - run: git push
GitHub Actions 非常适合 bioinformatics CI/CD，示例博客可参照 
hugobottois.medium.com
github.com
。

10. 单元测试与验收
代表数据集	期望输出	来源验证
PXD014397	I / Cancer / Melanoma	项目描述“Melanoma HLA-I ligandome” 
proteomecentral.proteomexchange.org
PXD022930	II / Normal / Cell Line	描述含 “Class II ligandome” 
proteomexchange.org
PXD025499	I/II / Infection / COVID-19	PRIDE keywords “SARS-CoV-2” 
proteomecentral.proteomexchange.org

pytest 断言分类结果，回归升级更安心。

参考文献与资源
PRIDE REST API 文档 
ebi.ac.uk

ppx Python 库文档 
ppx.readthedocs.io

ProteomeXchange 主页 
proteomexchange.org

PRIDE REST “additionalAttributes” 介绍 
pmc.ncbi.nlm.nih.gov

PRIDE Sample 端点示例 
ebi.ac.uk

joblib.Memory 使用教程 
joblib.readthedocs.io

Python 正则官方教程 
w3schools.com
docs.python.org

GitHub Actions 在生信中的应用 
hugobottois.medium.com
github.com

PXD014397 数据集页面 
proteomecentral.proteomexchange.org

PXD022930 数据集页面 
proteomexchange.org