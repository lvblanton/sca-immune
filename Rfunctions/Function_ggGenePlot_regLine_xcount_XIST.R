ggGenePlot_regLine_xcount_XIST <- function(Gene,all_data, normBy, prefix, myWidth=1.525, myHeight=1.75, myMin = -0.5, myMax = 10){
  ## Uncomment to troubleshoot
  # Gene <- "XIST"
  # normCounts <- lcl_normCounts
  # metadata <- lcl_metadata
  # resTable <- lcl_res
  # prefix<-"LCL"
  # adjustLabel <-"bottomRight"
  # normBy <- 1000
  # region<-"NPX"
  
  all_data$gene_expression <- all_data$gene_expression / normBy
  all_data$X_count <- all_data$X_count + 1
  all_data_xist <- all_data[all_data$X_count > 1,]
  
  pdf(file=paste0("Figures/",prefix,"_",Gene,"_linearPlot.pdf"), width=myWidth, height=myHeight)
  p <- ggplot() +
    geom_jitter(data=all_data, mapping = aes(x=X_count, y=gene_expression), color="#00000080", position=position_jitter(0.2), size=1.5, stroke=0) +
    geom_smooth(data=all_data_xist, mapping = aes(x=X_count, y=gene_expression), method = "lm", color=myOrange, size=0.25, fill=myOrange_light) +
    # annotate(
    #   geom = "text", 
    #   x = xloc, 
    #   y = yloc,
    #   label = paste0('dEX=', myDelta, '\nFDR=', myP),
    #   hjust = myH, 
    #   vjust = myV,
    #   size=2
    #   ) +
    labs( 
      x= "Number of Chr X",
      y= paste0("Normalized read\ncounts (x",normBy,")"),
      title = paste0("*",Gene,"*")
      ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = ggtext::element_markdown(hjust = 0.5),
      axis.text = element_text(color = "black"), 
      axis.ticks = element_line(color = "black", size = 0.25), 
      axis.line = element_line(size = 0.25),
      # aspect.ratio=1
      ) +
    ylim(myMin,myMax) + xlim(0.7,4.3)
  print(p)
  dev.off()
}
