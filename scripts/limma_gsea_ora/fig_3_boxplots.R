library(OlinkAnalyze)
library(plyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(patchwork)
library(RColorBrewer)
library(pheatmap)
library(matrixStats)

tmp <- read_NPX('WE-3707_NPX_2023-09-15.csv')

pheno <- fread('pheno_2023-09-26.txt')

additional_pheno <- fread('ALL_AML_merged_pheno_2024-01-09.csv')
additional_pheno[, sample_id:=gsub('/', '_', prov_id)]

keep_cols <- c(names(pheno), 'Status', 'RelapseYN', 'Dx_year')

pheno <- additional_pheno[pheno, on=c(sample_id='sample_id')]
pheno <- pheno[, ..keep_cols]

my.cols1 <- list(immunopheno=c('Control'='#FD9873', 'B-ALL'='#2F84AE', 'T-ALL'='#62BE9C', 'AML'='#8782D6'), 
		sex=brewer.pal(2, 'Set2'),
		age2=brewer.pal(5, 'RdPu')[2:5]) 

# left_join is from dplyr
foo1 <- as.data.table(left_join(x=tmp, y=pheno, by=join_by(SampleID == sample_id)))
foo1 <- foo1[QC_Warning == 'PASS' & Assay_Warning == 'PASS',]
foo1 <- foo1[!SampleID %in% c('09_228', '07_252', '07_154', 'K-023')]
foo1 <- foo1[, l_type:=plyr::revalue(immunopheno, c('T-ALL'='ALL', 'B-ALL'='ALL'))]
foo1[,dup:=sapply(strsplit(nopho_nr, '_'), '[', 2)]
foo1[is.na(dup), dup:=0]
foo1[, immunopheno:=factor(immunopheno, levels=names(my.cols1[['immunopheno']]))]
foo1[, age2:=cut(age, breaks=c(0,6,10,15,Inf), right=F)]
trans <- c('[0,6)'='<=5', '[6,10)'='6-9', '[10,15)'='10-14', '[15,Inf)'='>=15')
foo1[, age2:=revalue(age2, trans)]
foo1[immunopheno == 'Control' ,subtype:='Control']
foo <- foo1

med_npx <- foo[,.(med_npx=median(NPX)), by='SampleID']
foo[,SampleID:=factor(SampleID, levels=med_npx[['SampleID']][order(med_npx[['med_npx']])])]

foo <- foo[l_type %in% c('ALL', 'AML', 'Control')]
foo[, immunopheno:=revalue(immunopheno, c('B-ALL'='BCP-ALL'))]

#my.cols <- c(
#		Leukemia='#2a9d8f', 
#		AML = '#e9c46a',
#		'BCP-ALL'='#f4a261', 
#		'T-ALL'='#e76f51', 
#		Control='#264653' 
#	)
my.cols <- c(
		Control='#a6a6a6' ,
		Leukemia='#2a9d8f', 
		AML = '#efb2d1',
		'T-ALL'='#87ccb2', 
		'BCP-ALL'='#447597' 
	) 

genes <- fread('boxplot_genes.txt', header=F)[[1]]

comps <- list(
		'AML-Control'=list(use=expression(immunopheno %in% c('AML', 'Control')),
					tit='AML vs controls',
					col=my.cols
					),
		'B-Control'=list(use=expression(immunopheno %in% c('B-ALL', 'Control')),
					tit='BCP-ALL vs controls',
					col=my.cols
					),
		'T-Control'=list(use=expression(immunopheno %in% c('T-ALL', 'Control')),
					tit='T-ALL vs controls',
					col=my.cols
					)
)

read_tt <- function(n) {
	model <- 'subtype.age.sex'
	de <- fread(paste0('topTable_', n, '_.0.', model, '.tsv'))
	dups <- duplicated(de[,gene])	
	de <- de[!dups]
	de[,comparison:=n]
	return(de)
}

tts <- rbindlist(lapply(names(comps), read_tt))
tt <- tts[gene %in% genes]
tt[,immunopheno:=revalue(sapply(strsplit(comparison, '-'), '[', 1), c('T'='T-ALL', 'B'='BCP-ALL'))]

dd <- foo[Assay %in% genes, .(SampleID, Assay, Panel, NPX, sex, age2, immunopheno, subtype)] 
dd <- tt[dd, on=c(gene='Assay', panel='Panel', immunopheno='immunopheno')]
dd[, immunopheno:=factor(immunopheno, levels=c('Control', 'AML', 'T-ALL', 'BCP-ALL'))]

plot_fun <- function(d, mult=c(0, 0.1), p.size=3) {
	d_pv <- unique(d[,.(gene, immunopheno, P.Value)])
	d_pv[, P.Value:=ifelse(P.Value <= 0.05, as.character(signif(P.Value, 2)), 'ns')]
	p <- ggplot(d, aes(x=immunopheno, y=NPX)) + 
		geom_boxplot(aes(fill=immunopheno), outlier.shape=NA) + 
		geom_jitter(aes(fill=immunopheno), width=0.4, height=0, shape=21, 
						size=2, alpha=0.7, show.legend=F) + 
		scale_fill_manual(values=my.cols) + 
		scale_y_continuous(expand=expansion(mult=mult)) + 
		#geom_label(aes(x=immunopheno, y=Inf, label=P.Value), 
		#		data=d_pv, vjust=1.2, size=2.8, show.legend=F) + 
		geom_text(aes(x=immunopheno, y=Inf, label=P.Value), 
				data=d_pv, vjust=1.5, size=p.size, show.legend=F) + 
		facet_wrap(~gene) + 
		OlinkAnalyze::set_plot_theme() + 
		theme(axis.line.y=element_blank(), 
			axis.line.x=element_blank(), 
			strip.text.x=element_text(size=10, margin=margin(1, 0, 1, 0)), 
			axis.title.x=element_blank(),
			legend.title=element_blank(),
			#axis.text.x=element_text(angle=45, vjust=1, hjust=1),
			axis.text.x=element_blank(),
			axis.ticks.x=element_blank(),
			panel.border=element_rect(fill=NA, color="#333333")
			)
	return(p)
}
pl <- lapply(genes, function(g) plot_fun(dd[gene %in% g]))

wrapper <- function(plist, design, ax='collect', eq=T, mult=c(0, 0.1)) {
	rr <- lapply(plist, function(p) layer_scales(p)$y$range$range)
	rr <- rr[!sapply(rr, is.null)]
	yl <- c(min(sapply(rr, '[', 1)), max(sapply(rr, '[', 2)))
	print(yl)
	yl[1] <- yl[1] - (yl[2] - yl[1]) * mult[1]
	yl[2] <- yl[2] + yl[2] * mult[2]
	print(yl)
	if (eq) plist <- lapply(plist, function(p) p + ylim(yl))
	wp <- wrap_plots(plist, design=design, 
				guides='collect', axes=ax, axis_titles='collect')
	return(wp)
}

p0 <- readRDS('fig3a_gsea.rds')
p0 <- p0 + theme(strip.text.x=element_text(size=10, margin=margin(1, 0, 1, 0)),
			axis.text.y=element_text(size=8),
			axis.text.x=element_text(size=10),
		)

lo1 <- "
	ABC
	DEF
	G##
	"
pl1 <- c(pl[1:6], list(guide_area()))
#pl1[[3]] <- pl1[[3]] + theme(axis.text.x=element_text(margin=margin(2.2,0,-45,0)))
pl1[[1]] <- pl1[[1]] + ggtitle('Angiogenesis')

panel_b <- wrapper(pl1, design=lo1)
ggsave('fig4_boxplots_b_2025-05-15.pdf', panel_b , width=7.5, height=7.5)

pl11 <- c(pl[7:8], list(guide_area()))
pl11[[1]] <- pl11[[1]] + ggtitle('E2F Targets')
panel_c <- wrapper(pl11, design="ABC", eq=T)
ggsave('fig4_boxplots_c_2025-05-15.pdf', panel_c , width=7.5, height=3*0.9)

scl <- 0.8
pll <- lapply(genes[c(9,10,15,12,17)], function(g) plot_fun(dd[gene %in% g], p.size=2))
pl12 <- c(pll, list(guide_area()))
pl12 <- lapply(pl12, function(p) p + scale_y_continuous(expand=expansion(mult=c(.05, .2))))
panel_e <- wrapper(pl12, design="ABCDEF", eq=F)
ggsave('fig2_boxplots_e-j_one_row_2025-10-10.pdf', panel_e , width=3.5*4*scl, height=4.5*scl/1.5)

pll <- lapply(genes[9:14], function(g) plot_fun(dd[gene %in% g], p.size=3))
pl12 <- c(pll, list(guide_area()))
pl12 <- lapply(pl12, function(p) p + scale_y_continuous(expand=expansion(mult=c(.05, .2))))
d <-	"
	ABCD
	EFGH
	"
panel_e <- wrapper(pl12, design=d, eq=F)
ggsave('fig2_boxplots_e-j_two_rows_2025-05-15.pdf', panel_e , width=3.5*4*scl, height=4.5*scl)

lo2 <- "
	ABCD
	EFGH
	" 
pl2 <- lapply(pl[7:14], function(p) p + guides(fill='none'))
pl2[[1]] <- pl2[[1]] + ggtitle('Proteins that regulate E2F pathway')
pl2[[3]] <- pl2[[3]] + ggtitle('Proteins involved in innate immunity')
pl2[[5]] <- pl2[[5]] + ggtitle('Unique to AML')
pl2[[6]] <- pl2[[6]] + ggtitle('Unique to BCP-ALL')
pl2[[7]] <- pl2[[7]] + ggtitle('Unique to T-ALL')

panel_x <- wrapper(pl2, design=lo2, ax='collect_x', eq=F)
ggsave('fig3_boxplots_rest_2025-05-15.pdf', panel_x , width=14, height=5)

lo3 <- "
	ABB
	CCC
	"
ggsave('fig3_2025-05-15.pdf', wrap_plots(p0, panel_b, free(panel_x, side='l'), design=lo3), width=14, height=10)

mapply(function(x, y){print(y); ggsave(paste0('boxplot_', y, '_2025-05-15.pdf'), x, width=4.3, height=3)},  
									pl[14:18], genes[14:18])
