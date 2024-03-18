powerAnalysis <- function(data, meta, betas=c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0),
                          n.genes=100, n.reps=50, anno=NULL, expgenes=NULL,
                          parallel=TRUE) {
  # power analysis
  #
  # Parameters
  # ----------
  # data (dataframe) : genes x samples dataframe
  # meta (dataframe) : samples x variables dataframe
  # betas (numeric) : numeric vector of log2 effect sizes
  # n.genes (numeric) : number of genes to simulate as having an effect in 
  #   each replicate
  # n.reps (numeric) : number of replicates to perform
  # anno (dataframe) : gene to chromosome mapping; if given, only genes to 
  #   to be analyzed/recorded are autosomal genes
  # expgenes (character) : character vector of genes that are expressed; 
  #   effects will only be simulated for these genes
  
  # data <- data.fib
  # meta <- meta.fib
  # betas <- BETAS
  # n.genes <- N.GENES
  # n.reps <- N.REPS
  # anno <- anno
  # expgenes <- expgenes.fib
  # parallel=TRUE
  
  
  # set up labels <variable>|<beta>|<k>
  labels <- c()
  for (v in c("x_count", "y_count")) {
    for (b in betas) {
      labels <- c(labels, paste(v, paste(b, 1:n.reps, sep="|"), sep="|"))
    }
  }
  res <- mclapply(labels, calculatePower, data=data, meta=meta, ngenes=n.genes,
                  anno=anno, expgenes=expgenes, mc.cores=numCores)
  # reconstruct labels
  varcodes <- c("xcount", "ycount")
  labels <- paste(sapply(res, function (x) varcodes[x[["var_code"]]+1]),
                  sapply(res, function (x) x[["beta"]]),
                  sapply(res, function (x) x[["k"]]), sep="_")
  names(res) <- labels
  res <- t(data.frame(res))
  return(res)
}