library(clusterProfiler)
library(msigdbr)
library(data.table)
library(enrichplot)
library(patchwork)
library(ggplot2)

my.cols <- 'black'

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
					)),
		'HeH-t1221'=list(use=expression(subtype %in% c('HeH', 't1221')), 
				groups=list(
						'HeH'=expression(subtype %in% 'HeH'), 
						't1221'=expression(subtype %in% 't1221')
					))
		)

model <- 'subtype.age.sex'
cats <- c('H', paste0('C', 1:8))

for (n in names(comps)) {
	print(n)

	de <- fread(paste0('topTable_', n, '_.0.', model, '.tsv'))
	dups <- duplicated(de[,gene])	
	de <- de[!dups]

	#order.var <- 'logFC'
	order.var <- 't'
	setorderv(de, order.var, order=-1)


	for (ca in cats) {
		#t2g <- msigdbr(species='Homo sapiens', category='C5', subcategory='GO:BP') %>% 
		#						dplyr::select(gs_name, gene_symbol)
		t2g <- msigdbr(species='Homo sapiens', category=ca) %>% dplyr::select(gs_name, gene_symbol)
		
		#GSEA
		geneList <- de[[order.var]]
		names(geneList) <- de[['gene']]
		res_gsea <- GSEA(geneList, TERM2GENE = t2g)

		#ORA
		genes <- de[abs(logFC) > 1.5 & adj.P.Val < 0.05, gene]
		genes_up <- de[logFC > 1.5 & adj.P.Val < 0.05, gene]
		genes_down <- de[logFC < -1.5 & adj.P.Val < 0.05, gene]
		res_ora <- enricher(genes, universe=de$gene, TERM2GENE = t2g)
		res_ora_up <- enricher(genes_up, universe=de$gene, TERM2GENE = t2g)
		res_ora_down <- enricher(genes_down, universe=de$gene, TERM2GENE = t2g)

		res <- list(res_gsea, res_ora, res_ora_up, res_ora_down)
		names(res) <- c('GSEA', 'ORA', 'ORA upregulated', 'ORA downregulated')

		for (i in names(res)) {
			fwrite(as.data.table(res[[i]]), file=paste0(gsub(' ', '_', i), '_', n, '_', ca, '.tsv'))
		}

		res <- lapply(res, function(r) if (nrow(r)) pairwise_termsim(r) else r)
		#pl1 <- lapply(res, function(r) if (nrow(r)) treeplot(r, showCategory=30, 
		#					cex_category=0.01, fontsize=3) else ggplot() + theme_void())
		#pl2 <- lapply(res, function(r) if (nrow(r)) dotplot(r, showCategory=30, 
		#					, font.size=3) else ggplot() + theme_void())
		#plist <- list()
		#for (i in names(res)) {
		#	 plist[[i]] <- wrap_elements(wrap_plots(pl1[[i]], pl2[[i]], ncol=2) + plot_annotation(title=i))
		#}
		plist <- lapply(names(res), function(n) if (nrow(res[[n]])) dotplot(res[[n]], font.size=5, 
									showCategory=20) + 
									ggtitle(n) else ggplot() + theme_void())
		ggsave(paste0('gsea_ora_', n, '_', ca, '.pdf'), wrap_plots(plist, ncol=2) + plot_annotation(title=n),
										width=11, height=8)
	}
}
