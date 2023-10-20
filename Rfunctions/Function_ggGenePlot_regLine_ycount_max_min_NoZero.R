ggGenePlot_regLine_ycount_max_min_noZero <- function(Gene,male_data,all_data, normBy, prefix, myWidth=1.75, myHeight=1.75, myMin = -0.01, myMax=10){

  all_data$gene_expression <- all_data$gene_expression / normBy

  male_data$gene_expression <- male_data$gene_expression / normBy
    male_data$Y_count <- male_data$Y_count + 1


  pdf(file=paste0("Figures/",prefix,"_",Gene,"_linearPlot.pdf"), width=myWidth, height=myHeight)
  p <- ggplot() +
    geom_jitter(data=male_data, mapping = aes(x=Y_count, y=gene_expression), color="#00000080", position=position_jitter(0.2), size=1.5, stroke=0) +
    geom_smooth(data=male_data, mapping = aes(x=Y_count, y=gene_expression), method = "lm", color= myPurple, size=0.25, fill=myPurple_light) +
    labs( 
      x= "Number of Chr Y",
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
    ylim(myMin,myMax) + xlim(0.7,2.3)
  print(p)
  dev.off()
}
