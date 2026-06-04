library(data.table)

comps <- list(
		'AML-Control'=list(use=expression(immunopheno %in% c('AML', 'Control')),
					tit='AML vs controls'
					),
		'B-Control'=list(use=expression(immunopheno %in% c('B-ALL', 'Control')),
					tit='BCP-ALL vs controls'
					),
		'T-Control'=list(use=expression(immunopheno %in% c('T-ALL', 'Control')),
					tit='T-ALL vs controls'
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

p_lim <- 0.05
lfc_lim <- 1.5

tts[, Significant:=factor(ifelse(adj.P.Val < p_lim & logFC < -lfc_lim, "Downregulated", 
                  ifelse(adj.P.Val < p_lim & logFC > lfc_lim, "Upregulated", "Not significant")))]

tts[, Significant_bool:=as.numeric(Significant %in% c('Upregulated', 'Downregulated'))]

ret <- dcast(tts, 
			gene+panel~comparison, 
			#fun.aggregate=function(x) paste(x, collapse=':'), 
			value.var=c('logFC', 'adj.P.Val', 'Significant_bool'))

cn <- paste('Significant_bool', names(comps), sep='_')
ret[, unique:=apply(.SD, 1, sum) == 1, .SDcols=cn]

for (n in names(comps)) {
	use <- ret[,apply(.SD, 1, all), .SDcols=c('unique', paste('Significant_bool', n, sep='_'))]
	fwrite(ret[use, 1:8], file=paste0('DE_unique_to_', n, '_all_comparisons_included.tsv'), sep='\t')
	foo <- tts[ret[use, 1:2][, comparison:=n], on=c(gene='gene', panel='panel', comparison='comparison')]
	fwrite(foo[, .SD, .SDcols=!c('Significant_bool')], file=paste0('DE_unique_to_', n, '.tsv'), sep='\t')
}
