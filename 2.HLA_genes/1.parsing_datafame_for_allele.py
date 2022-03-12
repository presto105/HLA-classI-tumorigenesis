import numpy as np
import pandas as pd
from collections import Counter

df = pd.read_csv('TCGA_ann.tsv', sep='\t')
df = df.dropna(subset=['HLA-A','HLA-B','HLA-C']).reset_index(drop=True)

major_allele = {}
## classify major allele 
for gene in ['HLA-A','HLA-B','HLA-C'] :
    splited_allele = df[gene].str.split(',')
    unique_allele = splited_allele.map(set)

    merged_unique_allele = [ale for pt in unique_allele for ale in pt]
    allele_counter = Counter(merged_unique_allele).most_common()

    major_allele[gene] = [a for a,c in allele_counter if c > len(df) *0.01]


## append allele status by columns
for gene in ['HLA-A','HLA-B','HLA-C'] :
    splited_allele = df[gene].str.split(',')
    
    for ale_by_gene in major_allele[gene] :
        df[ale_by_gene] = 0
        df.loc[splited_allele.map(lambda x: ale_by_gene in x), ale_by_gene] = 1

df.to_csv('TCGA_ann_parsed_for_allele.tsv', sep='\t')
