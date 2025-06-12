#!/usr/bin/env python3
import pandas as pd

# 读取修复后的结果和真实标签
df_pred = pd.read_csv('dataset_annotation_fixed.tsv', sep='\t')
df_true = pd.read_csv('meta.txt', sep='\t')
df_true.columns = ['all_accession', 'HLA_true', 'Scenario_true', 'Disease_true']
df_compare = pd.merge(df_pred, df_true, on='all_accession')

# 计算准确率
accuracy_hla = (df_compare['HLA'] == df_compare['HLA_true']).mean()
accuracy_scenario = (df_compare['Scenario'] == df_compare['Scenario_true']).mean()
accuracy_disease = (df_compare['Disease'] == df_compare['Disease_true']).mean()

print("=== 修复后的结果 ===")
print(f'HLA准确率: {accuracy_hla:.2%}')
print(f'分析场景准确率: {accuracy_scenario:.2%}')
print(f'疾病类型准确率: {accuracy_disease:.2%}')

print(f'\n准确匹配的样本数:')
print(f'HLA: {(df_compare["HLA"] == df_compare["HLA_true"]).sum()}/{len(df_compare)}')
print(f'Scenario: {(df_compare["Scenario"] == df_compare["Scenario_true"]).sum()}/{len(df_compare)}')
print(f'Disease: {(df_compare["Disease"] == df_compare["Disease_true"]).sum()}/{len(df_compare)}')

# 显示一些正确分类的例子
print(f'\n=== 正确分类的例子 ===')
correct_hla = df_compare[df_compare['HLA'] == df_compare['HLA_true']]
for i, row in correct_hla.head(10).iterrows():
    print(f'{row["all_accession"]}: HLA {row["HLA"]} ✓, Scenario {row["Scenario"]} vs {row["Scenario_true"]}, Disease {row["Disease"]} vs {row["Disease_true"]}')

# 显示错误分类的例子
print(f'\n=== 错误分类的例子 ===')
wrong_hla = df_compare[df_compare['HLA'] != df_compare['HLA_true']]
for i, row in wrong_hla.head(5).iterrows():
    print(f'{row["all_accession"]}: HLA {row["HLA"]} vs {row["HLA_true"]} ✗')

print(f'\n=== 改进效果对比 ===')
print(f'原始准确率: HLA 6.02%, Scenario 4.82%, Disease 12.05%')
print(f'修复准确率: HLA {accuracy_hla:.2%}, Scenario {accuracy_scenario:.2%}, Disease {accuracy_disease:.2%}')
print(f'改进幅度: HLA +{accuracy_hla-0.0602:.2%}, Scenario +{accuracy_scenario-0.0482:.2%}, Disease +{accuracy_disease-0.1205:.2%}')