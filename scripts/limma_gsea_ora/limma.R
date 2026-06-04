library(limma)
library(OlinkAnalyze)
library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)
library(patchwork)
library(RColorBrewer)
library(plyr)
library(umap)
library(ggrepel)
library(ComplexUpset)
#library(ComplexHeatmap)

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


comps <- list(one=list(
			subset=expression(l_type %in% c('ALL', 'AML', 'Control')),
			contrasts=list('Leukemia-Control'=list(
						Leukemia=expression(unique(subtype[l_type %in% c('ALL', 'AML')])), 
						Control=expression(unique(subtype[l_type %in% 'Control']))
						),
					'AML-Control'=list(
						T=expression(unique(subtype[immunopheno %in% 'AML'])), 
						Control=expression(unique(subtype[immunopheno %in% 'Control']))
						),
					'B-Control'=list(
						B=expression(unique(subtype[immunopheno %in% 'B-ALL'])), 
						Control=expression(unique(subtype[immunopheno %in% 'Control']))
						),
					'T-Control'=list(
						T=expression(unique(subtype[immunopheno %in% 'T-ALL'])), 
						Control=expression(unique(subtype[immunopheno %in% 'Control']))
						),
					'ALL-AML'=list(
						B=expression(unique(subtype[l_type %in% 'ALL'])), 
						AML=expression(unique(subtype[immunopheno %in% 'AML']))
						),
					'B-AML'=list(
						B=expression(unique(subtype[immunopheno %in% 'B-ALL'])), 
						AML=expression(unique(subtype[immunopheno %in% 'AML']))
						),
					'T-AML'=list(
						T=expression(unique(subtype[immunopheno %in% 'T-ALL'])), 
						AML=expression(unique(subtype[immunopheno %in% 'AML']))
						),
					'B-T'=list(
						B=expression(unique(subtype[immunopheno %in% 'B-ALL'])), 
						T=expression(unique(subtype[immunopheno %in% 'T-ALL']))
						),
					'HeH-t1221'=list(
						HeH=expression('HeH'), 
						t1221=expression('t1221')
						)
					)
		)
)

# Limma  test
for (n in names(comps)) {
	cm <- comps[[n]]
	bar <- foo[eval(cm$subset),, drop=T]
	dd_wide <- dcast(bar, Assay + Panel ~ SampleID, value.var='NPX')
	gene_names <- paste(dd_wide$Assay, dd_wide$Panel, sep='.')
	dups <- dd_wide[duplicated(Assay), 1:2]
	edata <- as.matrix(dd_wide[,.SD, .SDcols=!c('Assay', 'Panel')])
	rownames(edata) <- gene_names

	anot <- unique(bar[, .(SampleID, l_type, immunopheno, subtype, sex, age, 
			age2, Status, RelapseYN, Dx_year)])[data.table(SampleID=colnames(edata)), on='SampleID']

	anot[, subtype:=gsub('-', '.', subtype)]

	#pp <- prcomp(t(edata[,colSums(is.na(edata)) < 1]), scale.=T, rank.=20, retx=T)
	pp <- prcomp(t(edata), scale.=T, rank.=20, retx=T)
	pd <- as.data.table(pp$x, keep.rownames='SampleID')
	dd <- anot[pd, on='SampleID']

	p1 <- ggplot(dd, aes(x=age, y=immunopheno)) + geom_violin() + geom_jitter(width=0.2)
	p2 <- ggplot(dd, aes(x=age, y=subtype)) + geom_violin() + geom_jitter(width=0.2)
	ggsave('age_vs_class.pdf', wrap_plots(p1, p2, nrow=1), width=8, height=4)

	use_vars <- c('immunopheno',
			'subtype', 
			'sex', 
			'age', 
			'age2',
			'Dx_year')

	plot_pc <- function(v) {
		pl <- list()
		for (pc in paste0('PC', 1:20)) {
			p <- ggplot(dd, aes(x=.data[[pc]], y=.data[[v]]))
			if (is.numeric(dd[[v]])) {
				ct <- cor.test(dd[[v]], dd[[pc]])
				p <- p + geom_point() + ggtitle(pc, 
								subtitle=paste(paste0('r=', signif(ct$estimate, 3)), 
										paste0('p=', signif(ct$p.value, 3)))) 
			} else {
				p <- p + geom_violin() + ggtitle(pc)
			}
			pl[[pc]] <- p + theme(plot.subtitle=element_text(size=5))
		}
		return(pl)
	}

	p_list <- lapply(use_vars, plot_pc)
	mapply(function(pl, tit) ggsave(paste0(tit, '_pca_scores.pdf'), wrap_plots(pl, ncol=4, nrow=5) + 
									plot_annotation(title=tit), 
										,width=8, height=11), 
										p_list, use_vars)
	ct_str <- vector(mode='character')
	for (ctn in names(cm$contrasts)) {
		ct <- cm$contrasts[[ctn]]
		g1 <- anot[,eval(ct[[1]])]
		g2 <- anot[,eval(ct[[2]])]
		s <- paste0('(', paste(paste0('subtype', g1), collapse='+'), ')/', length(g1), '-', 
				'(', paste(paste0('subtype', g2), collapse='+'), ')/', length(g2))
		names(s) <- ctn
		ct_str <- c(ct_str, s)
	}
	
	do_fit <- function(f) {
		mm <- model.matrix(as.formula(f), data=anot)
		ctr <- makeContrasts(contrasts=ct_str, levels=colnames(mm))
		colnames(ctr) <- names(cm$contrasts)
		fit_mm <- lmFit(edata, design=mm)
		fit_mm <- eBayes(fit_mm)
		fit_ctr <- contrasts.fit(fit_mm, ctr)
		fit_ctr <- eBayes(fit_ctr)
		return(list(fit_mm, fit_ctr))
	}
	
	formulas <- c('~0+subtype+age+sex',
			'~0+subtype+age',
			'~0+subtype')
	

	fits <- lapply(formulas, do_fit)



	p_lim <- 0.05
	lfc_lim <- 1.5

	dt <- lapply(fits, function(f) decideTests(f[[2]], p.value=p_lim, lfc=lfc_lim))

	mts <- lapply(names(cm$contrasts), function(cn) {
							res <- do.call('cbind', lapply(dt, function(d) d[,cn])); 
							colnames(res) <- formulas; 
							return(res)
							})
	names(mts) <- names(cm$contrasts)

	do_upset <- function(n) {
		m <- mts[[n]]

		tt <- fread(paste0(n, '_ttest.txt'))
		tmp <- tt[, .(assay=paste0(Assay, '.', Panel), 
				olink_ttest=ifelse(estimate < -lfc_lim & Adjusted_pval < p_lim, -1, 
					ifelse(estimate > lfc_lim & Adjusted_pval < p_lim, 1, 0)))]
		m <- merge(x=as.data.table(m, keep.rownames='assay'), y=tmp, by='assay', all=T, sort=F)
		m <- as.matrix(m[,.SD, .SDcols=!c('assay')])

		up <- as.data.frame(m > 0)
		down <- as.data.frame(m < 0)
		up <- up[rowSums(up) > 0,]
		down <- down[rowSums(down) > 0,]

		p1 <- upset(up, names(up), sort_sets=F, sort_intersections=F,
			base_annotation=list('Intersection size'=intersection_size(text=list(size=2.5))), 
			set_sizes=(upset_set_size() + 
				geom_text(aes(label=..count..), hjust=0, stat='count', size=3, colour='white') + 
				theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))),
						name='Linear model') + ggtitle('Upregulated')

		p2 <- upset(down, names(up), sort_sets=F, sort_intersections=F, 
			base_annotation=list('Intersection size'=intersection_size(text=list(size=2.5))), 
			set_sizes=(upset_set_size() + 
				geom_text(aes(label=..count..), hjust=0, stat='count', size=3, colour='white') + 
				theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))),
							name='Linear model') + ggtitle('Downregulated')
		panel <- wrap_plots(p1, p2, nrow=1) + plot_annotation(title=n)

		return(wrap_elements(panel))
	}

	plist <- lapply(names(mts), do_upset)
	ggsave(paste0('upset_', n, '_P_', p_lim, '_LFC_', lfc_lim, '.pdf'), wrap_plots(plist, ncol=1) + 
			plot_annotation(title=paste0('Adjusted P-value < ', p_lim, 
					', Log fold change > ', lfc_lim)), width=11, height=24)
	tts <- lapply(names(cm$contrasts), function(n) lapply(1:length(fits), 
						function(i) as.data.table(topTable(fits[[i]][[2]], coef=n, 
							number=nrow((fits[[i]][[2]]))), keep.rownames='gene.panel')))
	for (ci in 1:length(names(cm$contrasts))) {
		for (mi in 1:length(formulas)) {
			dd <- tts[[ci]][[mi]]
			dd[, c('gene', 'panel', 'gene.panel'):=list(sapply(strsplit(gene.panel, '\\.'), '[', 1), 
								sapply(strsplit(gene.panel, '\\.'), '[', 2), 
								NULL)]
			setcolorder(dd, c((ncol(dd)-1):ncol(dd), 1:(ncol(dd)-2)))
			write.table(dd, file=paste0('topTable_', 
							names(cm$contrasts)[ci], '_', 
							gsub('[+~]', '.', formulas[mi]), '.tsv'),
					sep='\t', quote=F, row.names=F)
		}
	}
}

