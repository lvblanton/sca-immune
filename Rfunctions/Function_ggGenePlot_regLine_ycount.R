ggGenePlot_regLine_ycount <- function(Gene,male_data,all_data, normBy, prefix, myWidth=1.525, myHeight=1.75){

  all_data$gene_expression <- all_data$gene_expression / normBy
  male_data$gene_expression <- male_data$gene_expression / normBy
  male_data$Y_count <- male_data$Y_count + 1

  pdf(file=paste0("Figures/",prefix,"_",Gene,"_linearPlot.pdf"), width=myWidth, height=myHeight)
  p <- ggplot() +
    geom_jitter(data=all_data, mapping = aes(x=Y_count, y=gene_expression), color="#00000080", position = position_jitter(0.2), size=1.5, stroke=0) +
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
    ylim(-0.01,max(all_data$gene_expression)+ 0.5) + xlim(-0.3,2.3)
  print(p)
  dev.off()
}
