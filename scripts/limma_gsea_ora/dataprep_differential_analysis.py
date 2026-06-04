# %%
import pandas as pd
import numpy as np

# %%
path = 'outpath/'

# %%
norm = pd.read_csv('outpath/GEX_matrix_merged_FINAL_n=352_AL.csv', index_col = 0)


# %%
proteins = pd.read_csv('outpath/Gene_protein_Olink_pairs.csv')
# The gene is the name from the Olink proteomics dataframe #Gene Updated is using translation with the GEX annotation file

# %%
len(proteins.gene_id.unique())

# %%
norm_filtered = norm[proteins.gene_id]
norm_filtered

# %%
norm_filtered.columns = proteins.gene

# %%
norm_filtered

# %%
sum(norm_filtered.columns == proteins.gene)

# %%
norm_filtered.T.to_csv('outpath/Filtered_gex_all_patients_Complete_Olink_Panel_AL.csv')


