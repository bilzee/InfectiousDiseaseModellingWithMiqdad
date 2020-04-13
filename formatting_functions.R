library("tidyverse")
library("plotly")
library("stringr")

maths_to_html = function(maths){
  html = str_squish(paste0(gsub("~","",maths) %>% gsub("\"","",.), collapse="")) %>% 
    gsub("\\[0\\]","<sub>0</sub>",.) %>%
    gsub("\\[","",.) %>%
    gsub("delta \\* t","dt",.) %>%
    gsub(" = ","=",.) %>%
    gsub(" ,",",",.) %>%
    gsub("R0","R<sub>0</sub>",.)
  
  return(html)
}

convert_ggplot_to_plotly = function(plot){
  main_title = maths_to_html(plot$labels$title)
  subtitle = maths_to_html(plot$labels$subtitle)
  caption = maths_to_html(plot$labels$caption)
  plotly = ggplotly(plot) %>%
    layout(title = list(text = paste0(main_title,
                                      '<br>',
                                      '<sup>',subtitle,'</sup>')),
           xaxis = list(title = maths_to_html(plot$labels$x)),
           yaxis = list(title = maths_to_html(plot$labels$y)),
           annotations = 
             list(x = 1, y = -0.4,
                  text = caption, 
                  showarrow = F, xref='paper', yref='paper', 
                  xanchor='right', yanchor='auto', xshift=0, yshift=0,
                  font=list(size=10, color="darkgreen", face="bold")),
           legend = list(orientation = "h", xanchor = "center", x = 0.5, y=-0.4))

  return(plotly)
}