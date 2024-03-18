SignificantResponse <- function(DESeqObject, Variable, expressedgenes){
  responsive <- results(DESeqObject, name = Variable, alpha=0.05)
  sigGenes <- rownames(subset(responsive, padj < 0.05))
  res.expressed <- responsive[rownames(responsive) %in% expressedgenes,]
  res.anno_expressed <- associateWrite(res.expressed)
  res.anno_expressed <- na.omit(res.anno_expressed)
  res.expressed <- res.expressed[rownames(res.expressed) %in% rownames(res.anno_expressed),]
  res.expressed_sig <- res.expressed[rownames(res.expressed) %in% sigGenes,]
  res.anno_expressed_sig <- associateWrite(res.expressed_sig)
  res.anno_expressed_sig <- na.omit(res.anno_expressed_sig)
  res.anno_expressed_sig <- res.anno_expressed_sig[res.anno_expressed_sig$gene_type.107 == "protein_coding" | res.anno_expressed_sig$gene_type.107 == "lncRNA",]
}
