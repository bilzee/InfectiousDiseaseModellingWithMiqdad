###################################################################################################
# Relationship exploring the properties of tests
# 
# Author: Miqdad Asaria
# Date: April 2020
###################################################################################################

library("tidyverse")
library("scales")

##################################################################################################
# Function to explore relationship between prevalence and test characteristics
#
# sensitivity = sensitivity of the test (how often it predicts true cases correctly)
# specificity = specificity of the test (how often it predicts false cases correctly)
# max_prevalence = range of disease prevalence over which to calculate test characteristics (zoom)
##################################################################################################
calculate_test_ppv_npv = function(sensitivity, specificity, max_prevalence){
  # calculate the positive predictive values and negative predictive values by prevalance
  test_reliability = tibble(prevalence = seq(from=0, to=max_prevalence, by=0.005),
         TruePositive = prevalence * sensitivity,
         FalseNegative = prevalence * (1 - sensitivity),
         FalsePositive = (1 - specificity) * (1 - prevalence),
         TrueNegative = specificity * (1 - prevalence),
         PPV = TruePositive / (TruePositive + FalsePositive),
         NPV = TrueNegative / (TrueNegative + FalseNegative))
  
  # reshape the data into tidy form for plotting
  graph_data = test_reliability %>% 
    select(prevalence, PPV, NPV) %>% 
    gather(criteria, value, c(PPV, NPV)) %>%
    mutate(criteria = factor(criteria, c("PPV", "NPV"),c("Positive Predictive Value", "Negative Predictive Value")))
  
  # plot PPVs and NPVs against prevalence
  test_plot = ggplot(graph_data, aes(x=prevalence, y=value, group=criteria, colour=criteria)) +
    geom_line(aes(linetype=criteria), size=1) +
    ggtitle(paste0("Test characteristics: sensitivity=",round(sensitivity*100,1),"% specificity=",round(specificity*100,1),"%")) +
    xlab("Prevalence of disease (%)") +
    ylab("Predictive accuracy (%)") + 
    scale_x_continuous(labels=scales::percent) +
    scale_y_continuous(labels=scales::percent) +
    scale_colour_manual(name="Criteria", values=c("green","red")) +
    scale_linetype_manual(name="Criteria", values=c(1,2)) +
    theme_bw(base_size = 10) + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  return(test_plot)
}