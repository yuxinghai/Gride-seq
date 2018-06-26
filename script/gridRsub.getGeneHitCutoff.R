## gridRsub.getGeneHitCutoff.R

## R code for chromatin-enriched genes


usage = 'Rscript gridRsub.getGeneHitCutoff.R <genecvg.bed>'
args = commandArgs(TRUE)

if(length(args) == 1) {
  
  x = sort(read.table(args[1], as.is = T, sep = '\t')[,5], decreasing = T)
  n = as.integer(length(x)/100)
  y = x[1:(length(x)-n)] - x[(n+1):length(x)]
  cat(x[max(which(y >= n))])
} else {
  cat(usage)
}
