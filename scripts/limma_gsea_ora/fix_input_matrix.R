library(data.table)


#d1 <- fread('translated_gex/GSE227832.filtered_genes.csv')
#d2 <- fread('translated_gex/stratmann2022.filtered_genes.csv')

#o1 <- fread('outpath/GEX_matrix_merged_FINAL_n=352.csv')
o2 <- fread('outpath/Gene_protein_Olink_pairs.csv')
#o3 <- fread('outpath/Filtered_gex_all_patients_Complete_Olink_Panel.csv')

p1 <- fread('outpath/Clinical_all_data_n=352.csv', header=T)

p_use <- p1[Final_Subtype %in% c('HeH', 't(12;21)') | group %in% c('T-ALL', 'AML', 'Control', 'T controls', 'B controls')]

d <- fread('translated_gex/normalized.csv')
#m_samples <- c(names(d1)[-1], names(d2)[-1])
#d <- cbind(d1, d2[,-1])

foo <- sub('^BM(\\d{4})$', 'BM\\1_I', p1[['V1']])
foo <- sub('^BM(\\d{4})-II$', 'BM\\1_II', foo)
foo <- sub('^AML(\\d{3})\\.D$','AML\\1-D', foo)

p1[,V1:=foo]

idx <- match(p1[['V1']], names(d))

p2 <- p1[!is.na(idx)]
dd <- d[, .SD, .SDcols=c(1, idx[!is.na(idx)])]
dd <- transpose(dd, make.names=1, keep.names='rn')

setnames(p2, c('', names(p2)[-1]))
setnames(dd, c('', names(dd)[-1]))
fwrite(dd, sep=',', quote=F, file='outpath/GEX_matrix_merged_FINAL_n=352_AL.csv')
fwrite(p2, sep=',', quote=F, file='outpath/Clinical_all_data_n=352_AL.csv')
