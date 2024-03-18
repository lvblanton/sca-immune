calculatePower <- function(label=NULL, data=NULL, meta=NULL, x=NULL, beta=NULL, 
                           ngenes=100, anno=NULL, expgenes=NULL, 
                           deseq.parallel=FALSE) {
  
  
  # label <- labels
  # data <- data
  # meta <- meta
  # ngenes <- n.genes
  # anno <- anno
  # expgenes <- expgenes
  # deseq.parallel = FALSE
  
  if (!is.null(label)) {
    items <- strsplit(label, "|", fixed=TRUE)[[1]]
    x <- items[1]
    beta <- as.numeric(items[2])
    k <- as.numeric(items[3])
  } else {
    k <- 0
  }
  
  # estimate power to detect simulated effects
  if (!(x %in% c("x_count", "y_count"))) {
    stop("x must be one of `x_count`, `y_count`")
  }
  
  # shuffle X, Y
  meta <- meta[colnames(data), ]
  cols <- c("x_count", "y_count")
  mc <- meta[, cols]
  meta <- meta[, !(colnames(meta) %in% cols)]
  # add back original
  mc1 <- mc
  colnames(mc1) <- paste(colnames(mc), "orig", sep="_")
  meta <- cbind(meta, mc1)
  # shuffle x/ys
  rownames(mc) <- sample(rownames(mc))
  meta <- cbind(meta, mc[rownames(meta), ])
  # record correlation with original x_count, y_count
  x_count_cor <- cor(meta$x_count, meta$x_count_orig)
  y_count_cor <- cor(meta$y_count, meta$y_count_orig)
  
  # select genes
  genes0 <- rownames(data)
  if (!is.null(expgenes)) {
    genes0 <- intersect(genes0, expgenes)
  } 
  if (!is.null(anno)) {
    # use only autosomal genes
    genes0 <- intersect(genes0, filter(anno, chrom=="auto")$gene)
  }
  genes <- sample(genes0, ngenes, replace=FALSE)
  
  # set up effect sizes
  mod <- model.matrix(as.formula(paste("~", x, "- 1")), meta)
  gene_effs <- data.matrix(rep(0, dim(data)[1]))
  rownames(gene_effs) <- rownames(data)
  gene_effs[genes, ] <- 1
  X <- gene_effs %*% t(mod)
  BX <- 2^(beta*X)
  
  # add effects
  data0 <- data
  data <- data * BX
  # record fractional increase in reads
  reads0 <- apply(data0, 2, sum)
  reads <- apply(data, 2, sum)
  frac_increase <- 2^mean(log2(reads / reads0))
  
  # run DESeq
  message("running deseq")
  dds <- DESeqDataSetFromMatrix(round(data), meta, ~x_count+y_count+batch_libprep)
  dds <- DESeq(dds, parallel=deseq.parallel)
  xres <- results(dds, name="x_count", alpha=0.05)
  yres <- results(dds, name="y_count", alpha=0.05)
  
  # count X and Y responsive genes
  xresp.genes <- rownames(xres[(!is.na(xres$padj)) & (xres$padj < 0.05), ])
  yresp.genes <- rownames(yres[(!is.na(yres$padj)) & (yres$padj < 0.05), ])
  if (!is.null(anno)) {
    # use only autosomal genes
    autogenes <- filter(anno, chrom=="auto")$gene
    xresp.genes <- intersect(xresp.genes, autogenes)
    yresp.genes <- intersect(yresp.genes, autogenes)
  }
  n.true.x <- length(intersect(xresp.genes, genes))
  n.true.y <- length(intersect(yresp.genes, genes))
  n.other.x <- length(setdiff(xresp.genes, genes))
  n.other.y <- length(setdiff(yresp.genes, genes))
  if (x == "x_count") {
    tpr <- n.true.x / ngenes
    #fdr <- n.other.x / (n.other.x + n.true.x)
    var_code <- 0
  } else if (x == "y_count") {
    tpr <- n.true.y / ngenes
    #fdr <- n.other.y / (n.other.y + n.true.y)
    var_code <- 1
  } else {
    tpr <- NA
    #fdr <- NA
    var_code <- 2
  }
  
  res <- c(var_code=var_code, k=k,
           beta=beta,  n_genes=ngenes, n_true_X=n.true.x, n_true_Y=n.true.y,
           n_other_X=n.other.x, n_other_Y=n.other.y, TPR=tpr)
  #FDR=fdr, x_count_cor=x_count_cor, y_count_cor=y_count_cor, 
  #frac_increase=frac_increase)
  return(res)
}