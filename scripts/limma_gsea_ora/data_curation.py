# %%
import pandas as pd
import numpy as np

# %%
path = 'outpath/'

# %%
norm = pd.read_csv('outpath/GEX_matrix_merged_FINAL_n=352_AL.csv', index_col = 0)

# %%
annot = pd.read_csv('HeH_vs_t1221/annotateddata.csv', sep = ';')

# %%
annot2 = annot[annot.feature == 'gene']
annot2['gene_biotype'] = annot2.gene_biotype.str.strip(';')
annot2.rename(columns = {'gene_name': 'gene'}, inplace = True)
annot2

# %%
proteins = pd.read_csv(path + 'heatmap_top10_by_-abs(logFC).tsv', sep = '\t')
proteins

# %%
print(proteins.shape)

# %%
print(len(proteins.gene.unique()))

# %%
gene_list = proteins.gene.unique()

# %%
annot2[annot2.gene.isin(gene_list)].shape

# %%
proteins = proteins.merge(annot2[['gene', 'gene_id']], how = 'left')
proteins

# %%
print(len(proteins.gene_id.unique()))

# %%
norm_filtered = norm[proteins.gene_id]
norm_filtered

# %%
norm_filtered.columns = proteins.gene

# %%
norm_filtered

# %%
print(sum(norm_filtered.columns == proteins.gene))

# %%
norm_filtered.T.to_csv('outpath/Filtered_gex_all_patients.csv')
proteins.to_csv('outpath/protein_info.csv')


