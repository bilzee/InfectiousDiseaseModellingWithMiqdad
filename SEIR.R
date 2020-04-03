####################################################################################################
# Simple SEIR model implementation
#
# Miqdad Asaria
# April 2020
###################################################################################################

library("tidyverse")

##################################################################################################
# Function to simulate simple SEIR model
#
# N = population
# R0 = basic reproductive number
# D_pre = average duration of pre-infectious period (days)
# D = average duration of infectious period (days)
# delta_t = time step of model (days)
# times = how many time steps to run the model for
##################################################################################################
simulate_SEIR = function(N = 100000, R0 = 2.2, D = 2.9, D_pre = 5.2, delta_t = 0.1, days = 100) {
  # start by using inputs to calculate key model prameters
  f = delta_t / D_pre
  r = delta_t / D
  beta = (R0 * r) / N 
  iterations = round(days/delta_t)
  
  # set up a markov trace to capture the number of people in each state (S,E,I,R) over time
  markov_trace = matrix(data = NA, nrow = iterations, ncol = 6, dimnames = list(NULL,c("time","S","E","I","R","Total_I")))
  # intial population in each state at time 0
  markov_trace[1,"time"] = 0
  markov_trace[1,"S"] = N-1
  markov_trace[1,"E"] = 1
  markov_trace[1,"I"] = 0
  markov_trace[1,"R"] = 0
  markov_trace[1,"Total_I"] = 0
  
  # run simulation for the desired number of iterations saving results into the markov trace
  for (t in 2:iterations) {
    # populations in each state in previous iteration
    time = markov_trace[t-1, "time"]
    S = markov_trace[t-1, "S"]
    E = markov_trace[t-1, "E"]
    I = markov_trace[t-1, "I"]
    R = markov_trace[t-1, "R"]
    Total_I = markov_trace[t-1, "Total_I"]
      
    # apply difference equations formulas to calculate new populations in each state
    markov_trace[t, "time"] = time + delta_t
    markov_trace[t, "S"] = S - (beta * I * S)
    markov_trace[t, "E"] = E + (beta * I * S) - (f * E)   
    markov_trace[t, "I"] = I + (f * E) - (r * I) 
    markov_trace[t, "R"] = R + (r * I) 
    markov_trace[t, "Total_I"] = Total_I + (f * E)
  }
  
  # convert matrix to tibble to make it easier to manipulate
  markov_trace = as_tibble(markov_trace)
  
  # calculate the herd immunity threshold
  HIT = (R0-1)/R0
  
  # look at markov trace to see when herd immunity achieved if at all
  if(markov_trace %>% filter(R>=HIT*N) %>% select(time) %>% nrow() > 0){
    HIT_date = markov_trace %>% filter(R>=HIT*N) %>% select(time) %>% min() %>% round()
  } else {
    HIT_date = NA
  }
  
  # reshape the data into tidy form for plotting
  graph_data = markov_trace %>%
    gather("state", "population", c("S","E","I","R")) %>%
    mutate(state = factor(state, 
                          c("S","E","I","R"), 
                          c("S = susceptible","E = pre-infections","I = infectious","R = recovered")))
  
  # plot the SEIR graph
  seir_plot = ggplot(graph_data, aes(x=time, y=population, group=state)) +
    xlab("Time (days since first case)") +
    ylab("Number of people") +
    geom_line(aes(colour=state, linetype=state), size=1) +
    labs(title = "Simple infectious disease model",
         subtitle = bquote("Input parameters:" ~ R[0] ~ " = " ~ .(R0) ~ ", D' = " ~ .(D_pre) ~ ", D = " ~ .(D) ~ ", " ~ delta*t ~ " = " ~ .(delta_t) ~ ", N = " ~ .(N) ~ " Derived parameters: " ~ beta ~ " = " ~ .(beta) ~ ", f = " ~ .(f) ~ ", r = " ~ .(r)),
         caption = paste0("Herd immunity threshold of ", round(100*HIT) ,"%",ifelse(is.na(HIT_date)," not achieved in this simulation", paste0(" achieved after ",HIT_date," days")))) +
    scale_colour_manual(name="States", values=c("red","orange","blue","green")) +
    scale_linetype_manual(name="States", values=c(1,2,1,1)) +
    theme_bw(base_size = 10) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  # return the entire markov trace the SIR plot and the HIT details as outputs of the function
  return(list("seir_markov_trace"=markov_trace,
              "seir_plot"=seir_plot,
              "HIT"=HIT,
              "HIT_date"=HIT_date))
}

