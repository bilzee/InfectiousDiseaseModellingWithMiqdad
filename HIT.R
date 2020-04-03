library("tidyverse")
library("scales")

calculate_HIT = function(R0){
  
  HIT = ((R0 - 1) / R0)
  
  return(HIT)
}

plot_HIT = function(){
  R0 = seq(from = 1, to = 10, by = 0.1)
  HIT = calculate_HIT(R0)
  
  hit_plot = ggplot(data=tibble(R0,HIT)) +
    geom_line(aes(x=R0, y=HIT)) +
    ggtitle("Herd immunity threshold as a function of R0") +
    xlab("R0") +
    ylab("HIT (%)") + 
    scale_y_continuous(labels=scales::percent) +
    theme_bw(base_size = 10) + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  return(hit_plot)
}