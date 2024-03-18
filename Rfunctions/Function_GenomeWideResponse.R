GenomeWideResponse <- function(DESeqObject, Variable, expressedgenes){
  responsive <- results(DESeqObject, name = Variable, alpha=0.05)
  sigGenes <- rownames(subset(responsive, padj < 0.05))
  res.expressed <- responsive[rownames(responsive) %in% expressedgenes,]
  res.anno_expressed <- associateWrite(res.expressed)
  res.anno_expressed <- na.omit(res.anno_expressed)
  res.anno_expressed <- res.anno_expressed[res.anno_expressed$gene_type.107 == "protein_coding" | res.anno_expressed$gene_type.107 == "lncRNA",]
}
