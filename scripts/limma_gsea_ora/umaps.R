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

proj_plot <- function(d, group, type='pca', tit='', col) {
	d <- d[!d$dup == '1',]
	if (type == 'pca') {
		p <- olink_pca_plot(df=d, color_g=group, byPanel=F, quiet=T)[[1]]
		dd <- as.data.table(p$data)
	} else if (type == 'umap') {
		p <- olink_umap_plot(df=d, color_g=group, byPanel=F, quiet=T)[[1]]
		dd <- as.data.table(p$data, keep.rownames='SampleID')
	}
	dd[d, on=.(SampleID), dup:=i.dup]
	p$data <- as.data.frame(dd)
	p$layers[[1]]$show.legend=T
	p <- p + aes(shape=dup) + scale_shape_manual(values = c('1' = 8, '2' = 16, '0' = 16)) + 
				scale_color_manual(values=col, limits=names(col), drop=F) + 
				guides(color=guide_legend(override.aes=list(size=3)), shape='none') + 
				theme(legend.title=element_blank()) 
	return(p + ggtitle(tit))
}

	
use <- foo[,!is.na(immunopheno)]
#my.cols <- c(
#		Leukemia='#8fce00', 
#		AML = '#8782D6',
#		'BCP-ALL'='#2F84AE', 
#		'T-ALL'='#62BE9C', 
#		Control='#FD9873' 
#	) 
my.cols <- c(
		Control='#a6a6a6' ,
		Leukemia='#2a9d8f', 
		AML = '#efb2d1',
		'T-ALL'='#87ccb2', 
		'BCP-ALL'='#447597' 
	) 

foo[, immunopheno2:=revalue(immunopheno, c('B-ALL'='BCP-ALL'))]

p0 <- proj_plot(as.data.frame(foo[use]), 'immunopheno2', 'pca', col=my.cols, tit='All samples')
#p2 <- proj_plot(as.data.frame(foo[use]), 'immunopheno2', 'umap', col=my.cols[-1])

#ggsave('fig1_B_pca.pdf', p1, width = 5.5, height = 3.5)
#ggsave('fig1_B_umap.pdf', p2, width = 5.5, height = 3.5)

comps <- list(
		'Leukemia-Control'=list(use=expression(l_type %in% c('ALL', 'AML', 'Control')),
					tit='Leukemia vs controls',
					col=my.cols,
					groups=list(
							Leukemia=expression(l_type %in% c('ALL', 'AML')), 
							Control=expression(l_type %in% 'Control')
					)),
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
					)),
		'B-AML'=list(use=expression(immunopheno %in% c('B-ALL', 'AML')),
					tit='BCP-ALL vs AML',
					col=my.cols[2:4],
					groups=list(
							'BCP-ALL'=expression(immunopheno %in% 'B-ALL'), 
							AML=expression(immunopheno %in% 'AML')
					)),
		'T-AML'=list(use=expression(immunopheno %in% c('T-ALL', 'AML')),
					tit='T-ALL vs AML',
					col=my.cols[2:4],
					groups=list(
							'T-ALL'=expression(immunopheno %in% 'T-ALL'), 
							AML=expression(immunopheno %in% 'AML')
					)),
		'B-T'=list(use=expression(immunopheno %in% c('T-ALL', 'B-ALL')),
					tit='BCP-ALL vs T-ALL',
					col=my.cols[2:4],
					groups=list(
							'BCP-ALL'=expression(immunopheno %in% 'B-ALL'), 
							'T-ALL'=expression(immunopheno %in% 'T-ALL')
					))
		)


volcano_plot <- function(d, p_lim=0.05, lfc_lim=1.5, tit) {

	my.cols <- c(Upregulated="#EEAF61", Downregulated="#78258F", "Not significant"='grey')

	d[, Significant:=factor(ifelse(adj.P.Val < p_lim & logFC < -lfc_lim, "Downregulated", 
                                 ifelse(adj.P.Val < p_lim & logFC > lfc_lim, "Upregulated", "Not significant")),
				levels=names(my.cols))]

	# Get the top 20 proteins with the lowest adjusted p-value
	top20_proteins <- d[!Significant %in% 'Not significant',][order(adj.P.Val), ][1:20, ]

	# Create the volcano plot
	p <- ggplot(d, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
  	geom_point() +
  	scale_color_manual(values=my.cols) +
	guides(color=guide_legend(override.aes=list(size=3))) + 
  	labs(
    		title = tit,
    		x = "Estimate",
    		y = "-log10(Adjusted p-value)"
  	) +
	OlinkAnalyze::set_plot_theme() + 
	theme(legend.title=element_blank()) + 
  	# Remove grid lines, background, and add axis lines
  	#theme_minimal() +
  	#theme(panel.grid.major = element_blank(), 
  	#      panel.grid.minor = element_blank(),
  	#      panel.background = element_blank(),  # Remove grey background
  	#      plot.background = element_blank()) +   # Remove plot background
  	#      panel.border = element_blank(),  
  	#      axis.line = element_line(color = "black")) +
  	# Add the labels for the top 20 proteins
  	geom_text_repel(data = top20_proteins, aes(label = gene), 
  	                box.padding = 0.35, point.padding = 0.3, size=3,
  	                segment.color = 'grey50', max.overlaps = 20, show.legend=F)
	return(p)
}

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
				
dotplot_mod <- function(d, tit) {
	d[,label:=label_fun(ID, 30)]
	d[,Count:=sapply(strsplit(core_enrichment, '/'), length)]
	d[,GeneRatio:=Count / setSize]
	setorder(d, GeneRatio)
	d[, label:=factor(label, levels=rev(unique(label)))]
	p <- ggplot(d, aes(x=label, y=GeneRatio, size=Count, fill=p.adjust)) + 
				geom_point(shape=21, colour='black') + 
				#coord_cartesian(ylim=c(0.1, 0.8)) + 
				scale_fill_viridis_c(name='Adjusted\np-value', breaks=seq(0.01, 0.03, 0.01), 
										limit=c(0, 0.04)) +
				scale_size(breaks=seq(10, 40, 10), limit=c(0, 50)) + 
				OlinkAnalyze::set_plot_theme()
	p <- p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=6), 
			axis.title.x=element_blank()) + 
	#		panel.background=element_blank(),
	#		panel.grid.minor=element_blank(),
	#		panel.grid.major.y=element_blank(),
	#		panel.grid.major.x=element_blank(),
	#		#panel.grid.major.x=element_line(colour='grey80', linewidth=0.5),
	#		panel.border=element_rect(fill=NA, colour='grey20')) +
		ggtitle(tit)
	return(p)
}
	


model <- 'subtype.age.sex'

plist <- list()
for (n in names(comps)) {
	print(n)
	cm <- comps[[n]]
	bar <- foo[eval(cm$use),, drop=T]
	for (vn in names(cm$groups)) {
		bar[eval(cm$groups[[vn]]), (n):=vn]
	}
	bar[[n]] <- factor(bar[[n]], levels=names(cm$groups))
	p1 <- proj_plot(as.data.frame(bar), n, type='pca', tit=cm$tit, col=cm$col)
	#p2 <- proj_plot(as.data.frame(bar), n, type='umap', tit=cm$tit, col=cm$col)

	tt <- fread(paste0('topTable_', n, '_.0.', model, '.tsv'))
	p3 <- volcano_plot(tt, tit=cm$tit)

	g <- fread(paste0('GSEA_', n, '_H.tsv'))
	if (nrow(g) > 0) {
		nr <- min(10, nrow(g))
		p4 <- dotplot_mod(g[1:nr,], cm$tit)
	} else {
		p4 <- ggplot() + theme_void()
	}

	#plist[[n]] <- list(pca=p1, umap=p2, volcano=p3, dot=p4)
	plist[[n]] <- list(pca=p1, volcano=p3, dot=p4)
}



exp_bar_plot <- function(comps, p_lim=0.05, lfc_lim=1.5, tit='') {
	dl <- lapply(names(comps), function(n) fread(paste0('topTable_', n, '_.0.', model, '.tsv')))
	names(dl) <- label_fun(sapply(comps, function(x) x$tit), wrap=10, split=' ', index=-(2:3))
	d <- rbindlist(dl, idcol='comparison')
	d[, comparison:=factor(comparison, level=names(dl))]

	my.cols <- c(Upregulated="#EEAF61", Downregulated="#78258F", "Not significant"='grey')

	d[, Significant:=factor(ifelse(adj.P.Val < p_lim & logFC < -lfc_lim, "Downregulated", 
                                 ifelse(adj.P.Val < p_lim & logFC > lfc_lim, "Upregulated", "Not significant")),
				levels=names(my.cols))]

	d <- d[!Significant == "Not significant"]
	# Drop duplicates from different panels when counting genes, keep the most significant
	d <- d[!duplicated(d[,c('comparison', 'gene')])]
	d2 <- d[, .(count=.N), by=c('Significant', 'comparison')]
	setorder(d2, comparison, -Significant)
	d3 <- d2[,.(count=count, csum=cumsum(count), Significant=Significant), by='comparison']
	p1 <- ggplot(d, aes(x=comparison, fill=Significant)) + 
			geom_bar(position='stack') +
			geom_text(mapping=aes(y=csum+5, x=comparison, label=count), data=d3, size=3, vjust=0) +
  			scale_fill_manual(values=my.cols) +
			guides(color=guide_legend(override.aes=list(size=3)), fill='none') + 
			labs(
				title = tit,
				x = NULL,
				y = "Number of genes"
			) +
			OlinkAnalyze::set_plot_theme() + 
			theme(legend.title=element_blank(), axis.text.x=element_text(angle=45, vjust=1, hjust=1))

	upset_data <- function(direction) {	
		ret <- dcast(d[Significant == direction, c(1,2,4)], 
					gene~comparison, 
					fun.aggregate=function(x) as.numeric(any(!is.na(x))), 
					value.var='logFC')
	}
	down <- upset_data('Downregulated')
	up <- upset_data('Upregulated')

	my.cols2 <- c(
		Leukemia='#2a9d8f', 
		AML = '#e9c46a',
		'BCP-ALL'='#f4a261', 
		'T-ALL'='#e76f51', 
		Control='#264653' 
	)

	upset_plot <- function(dd, tit, ss_lim) {
		qrs <- lapply(names(my.cols2), function(n) upset_query(set=n, fill=my.cols2[n]))
		p <- upset(dd, names(dd)[-1], sort_sets=F, sort_intersections='descending',
			base_annotation=list('Intersection size'=intersection_size(
									text_colors=c(on_background='black', 
											on_bar='black'),
									#text=list(size=2.5), 
									text=list(size=3.5), 
									fill='grey') + 
								OlinkAnalyze::set_plot_theme() + 
								theme(panel.grid=element_blank(),
									#axis.ticks.y=element_line(colour='grey30',
									#			linewidth=0.5),
									axis.title.x=element_blank(),
									axis.text.x=element_blank(),
									axis.ticks.x=element_blank(),
									axis.line.x=element_blank(),
									axis.line.y=element_blank()
									)
						), 
			#stripes=upset_stripes(
			#		geom=geom_segment(size=5),
			#		colors=unname(my.cols2)
			#	), 
			stripes='white',
			queries=qrs, 
			set_sizes=F, 
			#set_sizes=(upset_set_size() + 
			#	geom_text(aes(label=..count.., hjust=ifelse(..count.. < ss_lim, 1, 0)), 
			#					stat='count', size=2.5, colour='black') + 
			#	OlinkAnalyze::set_plot_theme() + 
			#	theme(#axis.text.x=element_text(angle=45, vjust=1, hjust=1), 
			#		#axis.title.x=element_blank(),
			#		#axis.ticks.x=element_line(colour='grey30',
			#		#			linewidth=0.5),
			#		#axis.line.x=element_line()
			#		axis.line.x=element_blank(),
			#		axis.title.y=element_blank(),
			#		axis.text.y=element_blank(),
			#		axis.text.x=element_blank(),
			#		axis.ticks.y=element_blank(),
			#		axis.ticks.x=element_blank(),
			#		axis.line.y=element_blank()
			#		)),
			matrix=intersection_matrix(geom=geom_point(size=2.5)),
			name=tit)# + ggtitle(tit)
	}
	p2 <- upset_plot(up, 'Upregulated', 100)
	p3 <- upset_plot(down, 'Downregulated', 1)

	
	return(list(p1, p2, p3))
} 
			
pl2 <- exp_bar_plot(comps[c(1:2, 4, 3)])


tr.pca <- wrap_plots(c(list(plot_spacer(), plot_spacer(), p0, guide_area()),
			lapply(plist[1:4], '[[', 'pca')), 
				nrow=2, 
				guides='collect', axis_titles='collect_y') & theme(legend.position='right') 
									#legend.justification.bottom='left')

#tr.umap <- wrap_plots(lapply(plist[1:4], '[[', 'umap'), nrow=1, 
#				guides='collect', axis_titles='collect_y') & theme(legend.position='bottom', 
#									legend.justification.bottom='left')
#lo <- "
#	AABBCCDD
#	EEFFFGGG
#	"
mr1 <- wrap_plots(lapply(plist[c(1:2, 4, 3)], '[[', 'volcano'), design='ABCD',
				guides='collect', axis_titles='collect_y') & theme(legend.position='bottom', 
									legend.justification.bottom='left')
mr2 <- wrap_plots(list(pl2[[1]], pl2[[2]], pl2[[3]]), design='AABBBCCC')

mr <- wrap_plots(mr1, mr2, ncol=1)
# heatmap (ph)
source('heatmap_significant_2024-12-10.R')

br <- wrap_plots(lapply(plist[1:4], '[[', 'dot'), nrow=1, 
				guides='collect', #axes='collect_y', 
						axis_titles='collect_y') & theme(legend.position='bottom', 
									legend.justification.bottom='left')

scl <- 0.8

ggsave('fig1_pca.pdf', tr.pca, width=3.6*4*scl, height=3*scl*2)
#ggsave('fig2_A_umap.pdf', tr.umap, width=3.5*4*scl, height=4*scl)
ggsave('fig2_A.pdf', mr, width=3.5*4*scl, height=4.5*scl*2)
ggsave(paste0('fig2_E_npx_unscaled_', as.character(order_by), '.pdf'), ph_npx_unscaled,  width=3.5*4*scl, height=4.5*scl*2)
ggsave(paste0('fig2_E_fc_unscaled_', as.character(order_by), '.pdf'), ph_fc_unscaled,  width=3.5*4*scl, height=4.5*scl*2)
ggsave(paste0('fig2_hemo_E_npx_unscaled_', as.character(order_by), '.pdf'), hemo_ph_npx_unscaled,  width=3.5*4*scl, height=4.5*scl*0.7)
ggsave(paste0('fig2_hemo_E_fc_unscaled_', as.character(order_by), '.pdf'), hemo_ph_fc_unscaled,  width=3.5*4*scl, height=4.5*scl*0.7)
ggsave('fig2_A.png', mr, width=3.5*4*scl, height=4.5*scl*2)
ggsave('fig2_C.pdf', br, width=3.5*4*scl, height=5*scl)

ggsave('fig2.pdf', wrap_plots(tr.pca, mr, br, ncol=1, heights=c(0.3, 0.3, 0.4)), width=3*length(plist)*scl, height=4*3*scl + 1)


tr.pca <- wrap_plots(lapply(plist[5:7], '[[', 'pca'), nrow=1, 
				guides='collect', axis_titles='collect_y') & theme(legend.position='bottom', 
									legend.justification.bottom='left')

#tr.umap <- wrap_plots(lapply(plist[5:7], '[[', 'umap'), nrow=1, 
#				guides='collect', axis_titles='collect_y') & theme(legend.position='bottom', 
#									legend.justification.bottom='left')
br <- wrap_plots(lapply(plist[5:7], '[[', 'volcano'), nrow=1, 
				guides='collect', axis_titles='collect_y') & theme(legend.position='bottom', 
									legend.justification.bottom='left')

ggsave('fig3_C_pca.pdf', tr.pca, width=3.5*3*scl, height=4*scl)
#ggsave('fig3_C_umap.pdf', tr.umap, width=3.5*3*scl, height=4*scl)
ggsave('fig3_D.pdf', br, width=3.5*3*scl, height=4.5*scl)
