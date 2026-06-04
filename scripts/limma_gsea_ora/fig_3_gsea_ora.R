library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(RColorBrewer)

my.cols <- c(
		Leukemia='#2a9d8f', 
		AML = '#e9c46a',
		'BCP-ALL'='#f4a261', 
		'T-ALL'='#e76f51', 
		Control='#264653' 
	) 

comps <- list(
	#	'Leukemia-Control'=list(use=expression(l_type %in% c('ALL', 'AML', 'Control')),
	#				tit='Leukemia vs controls',
	#				col=my.cols
	#				),
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

dotplot_prep_data <- function(d) {
	d[,label:=label_fun(ID, 30)]
	if (!'Count' %in% names(d)) {
		d[,Count:=sapply(strsplit(core_enrichment, '/'), length)]
	}
	if (!'GeneRatio' %in% names(d)) {
		d[,GeneRatio:=Count / setSize]
	} else {
		d[,GeneRatio:=sapply(strsplit(GeneRatio, '/'), function(x) as.numeric(x[1]) / as.numeric(x[2]))]
	}
	setorder(d, p.adjust)
	#if (nrow(d) > 0) {
	#	nr <- ifelse(nrow(d) > 10, 10, nrow(d))
	#}
	nr <- nrow(d)
	return(d)
	#return(d[1:nr, c('label', 'Count', 'GeneRatio', 'p.adjust')])
}
		
	
library(scales)
dotplot_mod <- function(d, tit) {
	#browser()
	d <- d[method == 'GSEA', .(.N, comb=paste(sort(comparison), collapse='-')), by=c('label')][d,on='label']
	d_na <- d[is.na(N), .(.N, comb=paste(sort(comparison), collapse='-')), 
					by=c('label')][d[is.na(N), .SD, .SDcols=!c('N', 'comb')], on='label']
	d_na[, N:=N+3]
	d <- rbind(d[!is.na(N)], d_na)
	setorder(d, N, comb, p.adjust)
	lvl <- d[, unique(label)]
	d[,label:=factor(label, levels=lvl)]
	#p <- ggplot(d, aes(x=comparison, y=label, size=Count, fill=p.adjust)) + 
	my.col <- rev(brewer.pal(n = 7, name = "RdYlBu"))
	browser()	
	p <- ggplot(d, aes(x=comparison, y=label, size=Count, fill=avg)) + 
				geom_point(shape=21) + 
				#geom_point(shape=21, colour='black') + 
				#coord_cartesian(ylim=c(0.1, 0.8)) + 
				#scale_color_manual(values=c('green', 'black')) + 
				scale_fill_gradient2(name='Average log2\nfold change', 
								low=my.col[1], 
								mid=my.col[4],
								high=my.col[7],
								midpoint=0,
								limits=c(-3,3),
								oob=scales::squish
										) +  
				#scale_fill_viridis_c(name='Adjusted\np-value', breaks=seq(0.01, 0.04, 0.01), 
				#						limit=c(0, 0.05)) +
				scale_size(breaks=seq(10, 40, 10), limit=c(0, 50)) + 
				OlinkAnalyze::set_plot_theme()
	p <- p + theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1, size=12), 
			axis.title.x=element_blank(),
			axis.title.y=element_blank()) + 
	#		panel.background=element_blank(),
	#		panel.grid.minor=element_blank(),
	#		panel.grid.major.y=element_blank(),
	#		panel.grid.major.x=element_blank(),
	#		#panel.grid.major.x=element_line(colour='grey80', linewidth=0.5),
	#		panel.border=element_rect(fill=NA, colour='grey20')) +
		ggtitle(tit)
	return(p)
}


dd <- list()
for (n in names(comps)) {
	print(n)
	cm <- comps[[n]]

	methods <- c('GSEA', 'ORA', 'ORA downregulated', 'ORA upregulated')
	ll <- lapply(gsub(' ', '_', methods), function(m) fread(paste0(m, '_', n, '_H.tsv')))
	ll <- lapply(ll, dotplot_prep_data)
	names(ll) <- methods
	d <- rbindlist(ll, idcol='method', fill=T)
	d[, comparison:=label_fun(cm$tit, 100, ' ', 1)]
	### Fold change hack ###
	model <- 'subtype.age.sex'
	de <- fread(paste0('topTable_', n, '_.0.', model, '.tsv'))
	dups <- duplicated(de[,gene])	
	de <- de[!dups]
	p_genes <- ifelse(d$method == 'GSEA', d$core_enrichment, d$geneID)
	gene_avg <- function(g) {
		avg <- de[gene %in% g, .(avg=mean(logFC), min=min(logFC), max=max(logFC))]
		return(avg)
	}
	gn <- rbindlist(lapply(strsplit(p_genes, '/'), gene_avg ))
	gn[, equal:= (min > 0) == (max > 0)]
	d <- cbind(d, gn)
	########################
	dd[[n]] <- d
}

dd <- rbindlist(dd)
dd[, comparison:=factor(comparison, levels=c('AML', 'T-ALL', 'BCP-ALL'))]

if(F) {
spl <- split(dd, by='method')
plist <- lapply(names(spl)[1:2], function(n) dotplot_mod(spl[[n]], n))
p <- wrap_plots(plist, nrow=1, guides='collect', 
				axis_titles='collect_y') #& theme(legend.position='bottom', 
							#legend.justification.bottom='left')

ggsave('fig3_pathway_1.pdf', p, width=10, height=5)
}

p2 <- dotplot_mod(dd[method %in% c('GSEA', 'ORA')], '') + 
				facet_wrap(~method, nrow=1) + 
				theme(axis.line=element_blank(),
					panel.border=element_rect(fill=NA , color='black'))
saveRDS(p2, file='fig3a_gsea.rds')
ggsave('fig3_pathway_2_2025-02-03.pdf', p2, width=6.5, height=6)
