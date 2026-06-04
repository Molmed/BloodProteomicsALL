library(data.table)
library(readxl)
library(patchwork)
library(ggplot2)

init_path = '../Data_for_Olga_2025-03-17/'
prot = fread(file.path(init_path , 'olink_proteins.txt'), sep = '\t')
uniprot = fread(file.path(init_path , 'Uniprot_location.txt'), sep = '\t')
uniprot = prot[uniprot, on='UniProt']

path = '../DEP_results_2024-12-19/'

npx = fread(file.path(init_path, 'WE-3707_NPX_2023-09-15.csv'), sep = ';')
npx = npx[, .(NPX=mean(NPX)), by=c('SampleID', 'Assay')]

qc_data = as.data.table(read_excel('../Summary_Pats-QC-pheno_JN-2024-06-04.xlsx'))
npx = qc_data[,.(sample_id, public_id)][npx, on=.(sample_id=SampleID)]

exclude_pid = c('AML_101', 'AML_139', 'ALL_920', 'K-023')
npx = npx[!public_id %in% exclude_pid]

clinical_olink = fread(file.path(init_path, 'pheno_2023-09-26.txt'), sep = '\t')


npx = clinical_olink[, .(sample_id, immunopheno)][npx, on=.(sample_id=sample_id), nomatch=NULL]

# Load gex early to check and remove overlap
gex = fread('../outpath/Filtered_gex_all_patients_Complete_Olink_Panel_AL.csv')
gex = melt(gex, id.vars='gene', variable.name='public_id', value.name='GEX')
# Don't remove overlap
#gex = gex[!public_id %in% npx[,public_id]]

npx_avg = npx[immunopheno != 'Control', .(NPX=mean(NPX)), by=c('immunopheno', 'Assay')] 
npx_avg[immunopheno == 'B-ALL', immunopheno:='BCP-ALL'] 

npx_avg <- uniprot[,.(Assay, location_short)][npx_avg, on=.(Assay=Assay)]

npx_avg[ is.na(location_short) | nchar(location_short) == 0, location_short:= 'Unclassified']


tt_path = dir(path, pattern='(AML|B|T)-Control.+\\.tsv', full.names=T)
df_final = lapply(tt_path, fread)
names(df_final) <- plyr::revalue(sapply(strsplit(basename(tt_path), '_'), '[', 2), 
				c('T-Control'='T-ALL', 'B-Control'='BCP-ALL', 'AML-Control'='AML'))
df_final = rbindlist(df_final, idcol='comparison')[abs(logFC)>1.5 & adj.P.Val < 0.05]
df_final <- df_final[,.(logFC=logFC[which.max(abs(logFC))]), by=.(gene, comparison)]

df_final = unique(npx_avg[,.(Assay, location_short)])[df_final, on=.(Assay=gene)]

# Select top 20 from combiantions of immunopheno and uniprot class
deps20 = df_final[, .SD[order(abs(logFC), decreasing=T)[1:min(20, .N)]], by=c('comparison', 'location_short')]

deps20 = npx_avg[deps20[,.SD, .SDcols=!c('location_short')], on=.(immunopheno=comparison, Assay=Assay)]

# All deps
deps  = npx_avg[df_final[,.SD, .SDcols=!c('location_short')], on=.(immunopheno=comparison, Assay=Assay)]

# Gene expression
clinical = fread('../outpath/Clinical_all_data_n=352_AL.csv', header=T)
gex = clinical[gex, .(public_id=V1, subtype=Final_Subtype, group, gene, GEX) , on=.(V1=public_id)]
gex = gex[subtype %in% c('HeH', 't(12;21)') | group %in% c('T-ALL', 'AML')]

gex_avg = gex[, .(GEX=mean(GEX)), by=c('group', 'gene')]


group_col = c('AML'='#e9c46a',
             'T-ALL'='#e76f51',
  		'BCP-ALL'='#f4a261' )

uniprot_col = c('Extracellular'='#66C2A599', 
		'Membrane'='#8DA0CB99', 
		'Intracellular'='#FC8D6299', 
		'Unclassified'='#80808099')


plot_fun <- function(d, n) {
	ct <- cor.test(~ NPX + GEX, method='spearman', data=d)
	pv <- signif(ct$p.value, 2)
	rho <- signif(ct$estimate, 2)
	p <- ggplot(d, aes(x=GEX, y=NPX, fill=location_short)) + 
		geom_point(shape=21, color='white', size=3, stroke=0.3) + 
		scale_fill_manual(values=uniprot_col) + 
		geom_smooth(aes(fill=NULL), method='lm', se=T, color='#5D8AA8', fill='#5D8AA8', alpha=0.2, linewidth=0.5) + 
		guides(fill='none') + 
		geom_text(x=Inf, y=-Inf, label=paste0('r=', rho, ', p=', pv), fontface='italic',
								 vjust=-1, hjust=1.1, size=3) + 
		labs(x='Mean GEX', y='Mean NPX', title=n) + 
		theme_bw() + 
		theme(panel.grid=element_blank(),
			axis.ticks=element_line(linewidth=0.3),
			panel.border=element_rect(fill=NA, linewidth=0.3),
			axis.title=element_text(size=10),
			plot.title=element_text(face='bold.italic', vjust=0, hjust=1, margin=margin(r=5, b=-15)))
	return(p)
}

group_order = c('AML', 'T-ALL', 'BCP-ALL')

dd = deps20[gex_avg, on=.(immunopheno=group, Assay=gene), nomatch=NULL]
spl = split(dd, by='immunopheno')
pl <- lapply(group_order, function(n) plot_fun(spl[[n]], n))
ggsave('test.pdf', wrap_plots(pl, ncol=3), width=8, height=3.2)

dd = deps[gex_avg, on=.(immunopheno=group, Assay=gene), nomatch=NULL]
spl = split(dd, by='immunopheno')
pl <- lapply(group_order, function(n) plot_fun(spl[[n]], n))
pl <- mapply(function(p, xul) p + scale_y_continuous(breaks=seq(-4, 10, 2)) + 
				scale_x_continuous(breaks=seq(-4, xul, 2)) + 
				coord_cartesian(ylim=c(-4, 10), xlim=c(-4, xul)), pl, c(12, 10, 8),
				SIMPLIFY=F)
#ggsave('Figure2_our_data_non_overlapping.pdf', wrap_plots(pl, ncol=3), width=8, height=3.2)
ggsave('Figure2_our_data_original_overlapping.pdf', wrap_plots(pl, ncol=3), width=8, height=3.2)
