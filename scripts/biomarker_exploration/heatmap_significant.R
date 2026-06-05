# runs from/after umaps_2024-12-06.R
library(pheatmap)
library(matrixStats)

foo <- foo[l_type %in% c('ALL', 'AML', 'Control')]

# Make wide format expression matrix and matching annotation
dd_wide <- dcast(foo, UniProt + Assay + Panel ~ SampleID, value.var='NPX')
#dups <- dd_wide[duplicated(Assay), 1:2]
#dd_wide <- dd_wide[!duplicated(Assay), ]
gene_names <- paste(dd_wide$Assay, dd_wide$Panel, sep='.')
edata <- as.matrix(dd_wide[,.SD, .SDcols=!c('UniProt', 'Assay', 'Panel')])
rownames(edata) <- gene_names

anot <- unique(foo[, .(SampleID, l_type, immunopheno, subtype, sex, age, 
		age2, Status, RelapseYN, Dx_year)])[data.table(SampleID=colnames(edata)), on='SampleID']
anot$immunopheno <- revalue(anot$immunopheno, c('B-ALL'='BCP-ALL'))

model <- 'subtype.age.sex'
use_comps <- 1:4
top_n <- 10
order_by <- expression(-abs(logFC))

# Read top ten rows from each comparison result
hm_assays <- lapply(names(comps[use_comps]), function(n) fread(paste0('topTable_', n, '_.0.', model, '.tsv')))#, 
												#nrows=top_n))
hm_assays <- lapply(hm_assays, function(tt) tt[abs(logFC) > 1.5 & adj.P.Val < 0.05][order(eval(order_by))][1:10])

names(hm_assays) <- label_fun(sapply(comps[use_comps], function(x) x$tit), wrap=50, split=' ', index=T)
olgas_list <- rbindlist(hm_assays, idcol='comparison')
fwrite(olgas_list, file=paste0('heatmap_top10_by_', as.character(order_by), '.tsv'), sep='\t')

## Hemo test ##
if(T) {
hemo <- fread('hemolysis_model_results.tsv')
hemo <- hemo[1:10, .(gene=Assay, panel=Panel)]
hm_assays <- c(hm_assays, list(Hemo=hemo))
}
tmp2 <- fread('Uniprot_location.txt')
loc2 <- unique(tmp2[, c('UniProt', 'location_short')])[unique(dd_wide[, .(UniProt, Assay)]), on=c('UniProt'='UniProt')]
names(loc2) <- c('UniProt', 'loc_uniprot', 'Assay')
###############

# Pick out those assays from the npx data
hm_data_spl <- lapply(hm_assays, function(asy) edata[asy[, paste0(gene, '.', panel)],])

add_to_rownames <- function(m, n) {
	rownames(m) <- paste(n, rownames(m), sep='.')
	return(m)
}

hm_data_spl <- lapply(names(hm_data_spl), function(n) add_to_rownames(hm_data_spl[[n]], n))
names(hm_data_spl) <- names(hm_assays)

reorder_rows <- function(d) {
	hc <- hclust(dist(t(scale(t(d))), method='euclidean'), method='complete')
	return(d[hc$order,])
}
# cluster rows within comparison
hm_data <- do.call('rbind', lapply(hm_data_spl, reorder_rows))

# prepare column and row annotations for pheatmap, it's matched to data matrix by row and column names
ac <- anot[, as.data.frame(.SD[,.(immunopheno, subtype, sex)], row.names=SampleID)]
ar <- data.frame(comparison=sapply(strsplit(rownames(hm_data), '\\.'), '[', 1),
		gene=sapply(strsplit(rownames(hm_data), '\\.'), '[', 2),
		panel=sapply(strsplit(rownames(hm_data), '\\.'), '[', 3),
		row.names=rownames(hm_data))

ar$location_uniprot <- loc2[as.data.table(ar)[, .(gene)], on=c(Assay='gene')][, loc_uniprot]
ar$location_uniprot[nchar(ar$location_uniprot) < 1] <- NA

my.colors <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

my.cols3 <- my.cols[2:5]
names(my.cols3) <- paste0(names(my.cols3), ' vs controls') 

## Hemo check ##
if (F) {
my.cols3 <- c(my.cols3, Hemo='pink')
}
################

a.color <- list(immunopheno=my.cols, 
		subtype=setNames(brewer.pal(11, 'Set3'), levels(factor(ac$subtype))), 
		location_uniprot=setNames(brewer.pal(3, 'Set2'), levels(factor(ar$location_uniprot))),
		sex=setNames(brewer.pal(3, 'Set1')[1:2], levels(factor(ac$sex))),
		comparison=my.cols3,
		panel=setNames(brewer.pal(8, 'Paired'), levels(factor(ar$panel))))

# Hack
#debug(pheatmap:::lo)
use <- !ar$comparison %in% 'Hemo'
use2 <- !names(hm_data_spl) %in% 'Hemo'
pnu <- pheatmap(hm_data[use,], scale='none', 
					silent=T, 
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[use,c(1,4), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data)[use], 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl[use2], nrow)),
					cluster_cols=T, 
					cluster_rows=F)
ph_npx_unscaled <- wrap_elements(pnu[['gtable']])

ctrls <- rownames(ac)[ac$immunopheno == 'Control']
ctrl_mean <- rowMeans(hm_data[,ctrls])

hm_data_fc <- hm_data - ctrl_mean 

pfu <- pheatmap(hm_data_fc[use,], scale='none', 
					silent=T,
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[use,c(1,4), drop=F],
					#annotation_row=ar[use,c(1), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data_fc)[use], 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl[use2], nrow)),
					cluster_cols=T, 
					cluster_rows=F)
ph_fc_unscaled <- wrap_elements(pfu[['gtable']])

## Hack
use <- ar$comparison %in% 'Hemo'
use2 <- names(hm_data_spl) %in% 'Hemo'
a.color$comparison <- c('Hemo'='pink')

hemo_pnu <- pheatmap(hm_data[use,], scale='none', 
					silent=T, 
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[use,c(1,4), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data)[use], 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl[use2], nrow)),
					cluster_cols=pnu$tree_col, 
					cluster_rows=F)
hemo_ph_npx_unscaled <- wrap_elements(hemo_pnu[['gtable']])

hemo_pfu <- pheatmap(hm_data_fc[use,], scale='none', 
					silent=T, 
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[use,c(1,4), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data_fc)[use], 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl[use2], nrow)),
					cluster_cols=pfu$tree_col, 
					cluster_rows=F)
hemo_ph_fc_unscaled <- wrap_elements(hemo_pfu[['gtable']])
