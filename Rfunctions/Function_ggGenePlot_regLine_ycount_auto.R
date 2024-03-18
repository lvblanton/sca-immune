ggGenePlot_regLine_ycount_auto <- function(gene,all_data, normBy, prefix, myWidth, myHeight){
  #all_data is dataframe of normalized counts with associated metadata; rows are samples and columns are genes/metadata
  all_data[,"Y_count"] <- as.numeric(all_data[,"Y_count"])
  all_data[,gene] <- as.numeric(all_data[,gene] /normBy)


  data <- all_data[,c("Y_count", gene)]
  colnames(data) <- c("Y_count", "gene_expression")
  pdf(file=paste0(myPath, "/Figures/",prefix,"_",gene,"_YlinearPlot.pdf"), width=myWidth, height=myHeight)
  p <- ggplot(data = data, aes(x = Y_count, y = gene_expression)) +
    geom_jitter(color="#00000080", position=position_jitter(0.2), size=1.5, stroke=0) +
    geom_smooth(method = "lm", color=myPurple, size=0.25, fill=myPurple_light, fullrange=TRUE) +
    xlim(-0.3,2.3) +
    labs(
      x= "Number of Chr Y",
      y= paste0("Normalized read\ncounts (x",normBy,")"),
      title = gene
      ) +
    #labs(title = gene, x="Number of X chromosomes", y = "Normalized read counts") +
    #geom_vline(xintercept=0, col="black") +
    #geom_hline(yintercept=0, col="black")+
    theme_classic(base_size = 7)
  print(p)
  dev.off()
}
