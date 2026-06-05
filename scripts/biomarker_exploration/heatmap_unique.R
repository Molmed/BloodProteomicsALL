# NOT anymor anymoree runs from/after umaps_2024-12-06.R
library(OlinkAnalyze)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(patchwork)
library(RColorBrewer)
library(umap)
library(ggrepel)
library(ComplexUpset)

label_fun <- function(id, wrap=20, split='_', index=-1) {
	foo <- lapply(strsplit(as.character(id), split), '[', index)
	wrap_fun <- function(parts) {
		s <- parts[1]
		for (p in parts[-1]) {
			nc <- nchar(s) + nchar(p) + 1
			if (nc > wrap) {
				sep <- '\n'
			} else {
				sep <- ' '
			}
			s <- paste0(s, sep, p)
		}
		return(s)
	}
	ret <- sapply(foo, wrap_fun)
	return(ret)
}
tmp <- read_NPX('WE-3707_NPX_2023-09-15.csv')

pheno <- fread('pheno_2023-09-26.txt')

additional_pheno <- fread('ALL_AML_merged_pheno_2024-01-09.csv')
additional_pheno[, sample_id:=gsub('/', '_', prov_id)]

keep_cols <- c(names(pheno), 'Status', 'RelapseYN', 'Dx_year')

pheno <- additional_pheno[pheno, on=c(sample_id='sample_id')]
pheno <- pheno[, ..keep_cols]

# left_join is from dplyr
foo1 <- as.data.table(left_join(x=tmp, y=pheno, by=join_by(SampleID == sample_id)))
foo1 <- foo1[QC_Warning == 'PASS' & Assay_Warning == 'PASS',]
foo1 <- foo1[!SampleID %in% c('07_252', '07_154', 'K-023')]
foo1 <- foo1[, l_type:=plyr::revalue(immunopheno, c('T-ALL'='ALL', 'B-ALL'='ALL'))]
foo1[,dup:=sapply(strsplit(nopho_nr, '_'), '[', 2)]
foo1[is.na(dup), dup:=0]
#foo1[, immunopheno:=factor(immunopheno, levels=names(my.cols1[['immunopheno']]))]
foo1[, age2:=cut(age, breaks=c(0,6,10,15,Inf), right=F)]
trans <- c('[0,6)'='<=5', '[6,10)'='6-9', '[10,15)'='10-14', '[15,Inf)'='>=15')
foo1[, age2:=revalue(age2, trans)]
foo1[immunopheno == 'Control' ,subtype:='Control']
foo <- foo1

med_npx <- foo[,.(med_npx=median(NPX)), by='SampleID']
foo[,SampleID:=factor(SampleID, levels=med_npx[['SampleID']][order(med_npx[['med_npx']])])]

	
use <- foo[,!is.na(immunopheno)]

#my.cols <- c(
#		#Leukemia='#2a9d8f', 
#		AML = '#e9c46a',
#		'BCP-ALL'='#f4a261', 
#		'T-ALL'='#e76f51', 
#		Control='#264653' 
#	) 

my.cols <- c(
		Control='#a6a6a6' ,
		#Leukemia='#2a9d8f', 
		AML = '#efb2d1',
		'T-ALL'='#87ccb2', 
		#BCP-ALL'='#447597',
		'HeH'='#447597',
		'ETV6:RUNX1'='#88c2eb'
	) 


comps <- list(
		'AML-Control'=list(use=expression(immunopheno %in% c('AML', 'Control')),
					tit='AML vs controls',
					col=my.cols,
					groups=list(
							AML=expression(immunopheno %in% 'AML'), 
							Control=expression(immunopheno %in% 'Control')
					)),
		'B-Control'=list(use=expression(immunopheno %in% c('B-ALL', 'Control')),
					tit='BCP-ALL vs controls',
					col=my.cols,
					groups=list(
							'BCP-ALL'=expression(immunopheno %in% 'B-ALL'), 
							Control=expression(immunopheno %in% 'Control')
					)),
		'T-Control'=list(use=expression(immunopheno %in% c('T-ALL', 'Control')),
					tit='T-ALL vs controls',
					col=my.cols,
					groups=list(
							'T-ALL'=expression(immunopheno %in% 'T-ALL'), 
							Control=expression(immunopheno %in% 'Control')
					))
		)

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
anot[immunopheno == 'BCP-ALL', immunopheno:=subtype]
anot[, immunopheno:=revalue(immunopheno, c('t1221'='ETV6:RUNX1'))]

a.ord <- order(factor(anot$immunopheno, levels=names(my.cols)), anot$subtype)
edata <- edata[, a.ord]
anot <- anot[a.ord]

top_n <- 10
order_by <- expression(-abs(logFC))

# Read top ten rows from each comparison result
hm_assays <- lapply(names(comps), function(n) fread(paste0('DE_unique_to_', n, '.tsv')))#, 
												#nrows=top_n))
hm_assays <- lapply(hm_assays, function(tt) tt[order(eval(order_by))][1:10])

names(hm_assays) <- label_fun(sapply(comps, function(x) x$tit), wrap=50, split=' ', index=T)

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
hm_data_spl <- hm_data_spl[c(1,3,2)]

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

my.cols3 <- c(
		AML = '#efb2d1',
		'T-ALL'='#87ccb2', 
		'BCP-ALL'='#447597'
	) 

names(my.cols3) <- paste0(names(my.cols3), ' vs controls') 


a.color <- list(immunopheno=my.cols, 
		subtype=setNames(brewer.pal(11, 'Set3'), levels(factor(ac$subtype))), 
		location_uniprot=setNames(brewer.pal(3, 'Set2'), levels(factor(ar$location_uniprot))),
		sex=setNames(brewer.pal(3, 'Set1')[1:2], levels(factor(ac$sex))),
		comparison=my.cols3,
		panel=setNames(brewer.pal(8, 'Paired'), levels(factor(ar$panel))))

#debug(pheatmap:::lo)
pnu <- pheatmap(hm_data, scale='none', 
					silent=T, 
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[,c(1), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data), 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl, nrow)),
					cluster_cols=T, 
					cluster_rows=F)
ph_npx_unscaled <- wrap_elements(pnu[['gtable']])

ctrls <- rownames(ac)[ac$immunopheno == 'Control']
ctrl_mean <- rowMeans(hm_data[,ctrls])

hm_data_fc <- hm_data - ctrl_mean 

pfu <- pheatmap(hm_data_fc, scale='none', 
					silent=T,
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[,c(1), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data_fc), 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl, nrow)),
					cluster_cols=T, 
					cluster_rows=F)
ph_fc_unscaled <- wrap_elements(pfu[['gtable']])

scl <- 0.8
ggsave(paste0('unique_proteins_npx_unscaled_', as.character(order_by), '.pdf'), ph_npx_unscaled,  width=3.5*4*scl, height=4.5*scl*2)
ggsave(paste0('unique_proteins_fc_unscaled_', as.character(order_by), '.pdf'), ph_fc_unscaled,  width=3.5*4*scl, height=4.5*scl*2)

pnu <- pheatmap(hm_data, scale='none', 
					silent=T, 
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[,c(1), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data), 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl, nrow)),
					cluster_cols=F, 
					cluster_rows=F)
ph_npx_unscaled <- wrap_elements(pnu[['gtable']])

pfu <- pheatmap(hm_data_fc, scale='none', 
					silent=T,
					color=my.colors, 
					annotation_colors=a.color, 
					show_colnames=F, 
					legend=T, 
					fontsize=8, 
					annotation_col=ac[,1, drop=F],
					annotation_row=ar[,c(1), drop=F],
					annotation_names_col=F,
					annotation_names_row=F,
					labels_row=ar[rownames(hm_data_fc), 'gene'],
					gaps_row=cumsum(sapply(hm_data_spl, nrow)),
					cluster_cols=F, 
					cluster_rows=F)
ph_fc_unscaled <- wrap_elements(pfu[['gtable']])

ggsave(paste0('unique_proteins_npx_unscaled_unclustered_', as.character(order_by), '.pdf'), ph_npx_unscaled,  width=3.5*4*scl, height=4.5*scl*2)
ggsave(paste0('unique_proteins_fc_unscaled_unclustered_', as.character(order_by), '.pdf'), ph_fc_unscaled,  width=3.5*4*scl, height=4.5*scl*2)
