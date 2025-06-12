#!/usr/bin/env python3
import pandas as pd

# 读取预测结果和真实标签
df_pred = pd.read_csv('dataset_annotation.tsv', sep='\t')
df_true = pd.read_csv('meta.txt', sep='\t')
df_true.columns = ['all_accession', 'HLA_true', 'Scenario_true', 'Disease_true']
df_compare = pd.merge(df_pred, df_true, on='all_accession')

# 显示非Unspecified的预测结果
non_unspec = df_compare[df_compare['HLA'] != 'Unspecified']
print('成功分类的数据集:')
for _, row in non_unspec.iterrows():
    print(f'{row["all_accession"]}: {row["HLA"]} vs {row["HLA_true"]}, {row["Scenario"]} vs {row["Scenario_true"]}, {row["Disease"]} vs {row["Disease_true"]}')

print(f'\n成功获取API数据的数据集数量: {len(non_unspec)}/83')

# 计算准确率
accuracy_hla = (df_compare['HLA'] == df_compare['HLA_true']).mean()
accuracy_scenario = (df_compare['Scenario'] == df_compare['Scenario_true']).mean()
accuracy_disease = (df_compare['Disease'] == df_compare['Disease_true']).mean()

print(f'\n总体准确率:')
print(f'HLA准确率: {accuracy_hla:.2%}')
print(f'分析场景准确率: {accuracy_scenario:.2%}')
print(f'疾病类型准确率: {accuracy_disease:.2%}')

# 分析问题
print(f'\n问题分析:')
print(f'大部分数据集API返回404错误，无法获取样本信息')
print(f'可能的原因：')
print(f'1. 数据集已被撤回或移动')
print(f'2. API端点发生变化')
print(f'3. 需要更新的API访问方式')