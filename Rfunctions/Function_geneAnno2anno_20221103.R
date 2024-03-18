geneAnno2anno <- function(geneAnno) {
  # reformat geneAnno dataframe for use in analyses here
  anno <- data.frame(gene=geneAnno$Gene, 
                     chrom=gsub("chr[0-9]+", "auto", geneAnno$chr))
  return(anno)
}