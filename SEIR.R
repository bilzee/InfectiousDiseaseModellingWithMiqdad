####################################################################################################
# Simple SEIR model implementation
#
# Miqdad Asaria
# April 2020
###################################################################################################

library("tidyverse")
library("deSolve")

# calculate the herd immunity threshold
calculate_HIT = function(R0){
  HIT = ((R0 - 1) / R0)
  return(HIT)
}

# look up date that the herd immunity threshold is reached
calculate_HIT_date = function(markov_trace, R0, N){ 
  HIT = calculate_HIT(R0)
  # look at markov trace to see when herd immunity achieved if at all
  if(markov_trace %>% filter(S<=(1-HIT)*N) %>% select(time) %>% nrow() > 0){
    HIT_date = markov_trace %>% filter(S<=(1-HIT)*N) %>% select(time) %>% min() %>% round()
  } else {
    HIT_date = NA
  }
  return(HIT_date)
}

plot_HIT = function(){
  R0 = seq(from = 1, to = 10, by = 0.1)
  HIT = calculate_HIT(R0)
  
  hit_plot = ggplot(data=tibble(R0,HIT)) +
    geom_line(aes(x=R0, y=HIT)) +
    labs(title=bquote("Herd immunity threshold as a function of " ~ R[0])) +
    xlab(bquote(R[0])) +
    ylab("HIT (%)") + 
    scale_y_continuous(labels=scales::percent) +
    theme_bw(base_size = 10) + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  return(hit_plot)
}

##################################################################################################
# Function to simulate simple SEIR model using differential equations
#
# N = population
# R0 = basic reproductive number
# D_pre = average duration of pre-infectious period (days)
# D = average duration of infectious period (days)
# days = how many days to run the model for
# S_init = inital number of people susceptible 
# E_init = inital number of people pre-infectious 
# I_init = inital number of people infectious 
# R_init = inital number of people immune 
# Total_I_init = inital cumulative number of people infectious
##################################################################################################
simulate_SEIR_differential = function(N = 100000, R0 = 2.2, D = 2.9, D_pre = 5.2, days = 150, 
                              S_init = N-1, E_init = 0, I_init = 1, R_init = 0, Total_I_init = 0) {
  # start by using inputs to calculate key model prameters
  f = 1 / D_pre
  r = 1 / D
  beta = (R0 * r) / N 
  
  # time sequence to produce output for 
  time = seq(0, days)
  
  # parameters
  parameters = c("r"=r, "f"=f, "beta"=beta)
  
  # initial condition
  initial_state = c("S"=S_init, "E"=E_init, "I"=I_init, "R"=R_init, "Total_I"=Total_I_init)
  
  # calculate the value of the derivatives at each time value
  seir = function(t, state, parameters){
    with(as.list(c(state, parameters)), {
      dS = -(beta * I * S)
      dE = (beta * I * S) - (f * E)   
      dI = (f * E) - (r * I) 
      dR = (r * I) 
      dTotal_I = (f * E)
      return(list(c(dS, dE, dI, dR, dTotal_I)))
    })
  }
  
  # run ode to numerically integrate the system of equations
  markov_trace = ode(y=initial_state, times=time, func=seir, parms=parameters)

  return(as_tibble(markov_trace))
}

probability_to_rate = function(p,t){
  r = -log(1-p)/t
  return(r)
}

rate_to_probability = function(r,t){
  p = 1-exp(-r*t)
  return(p)
}

##################################################################################################
# Function to simulate simple SEIR model using difference equations
#
# N = population
# R0 = basic reproductive number
# D_pre = average duration of pre-infectious period (days)
# D = average duration of infectious period (days)
# delta_t = time step of model (days)
# days = how many days to run the model for
##################################################################################################
simulate_SEIR_difference = function(N = 100000, R0 = 2.2, D = 2.9, D_pre = 5.2, delta_t = 0.1, days = 150, 
                         time_init = 0, S_init = N-1, E_init = 0, I_init = 1, R_init = 0, Total_I_init = 0) {
  
  # start by using inputs to calculate key model prameters
  f = delta_t / D_pre
  r = delta_t / D
  beta = (R0 * r) / N 
  iterations = round(days/delta_t) + 1
  
  # set up a markov trace to capture the number of people in each state (S,E,I,R) over time
  markov_trace = matrix(data = NA, nrow = iterations, ncol = 6, dimnames = list(NULL,c("time","S","E","I","R","Total_I")))
  # intial population in each state at time 0
  markov_trace[1,"time"] = time_init
  markov_trace[1,"S"] = S_init
  markov_trace[1,"E"] = E_init
  markov_trace[1,"I"] = I_init
  markov_trace[1,"R"] = R_init
  markov_trace[1,"Total_I"] = Total_I_init
  
  if(iterations > 1){
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
  }
  # convert matrix to tibble to make it easier to manipulate
  markov_trace = as_tibble(markov_trace)
 
  return(markov_trace) 
}

suppress_release_SEIR = function(N, R0, R0_sup, D, D_pre, D_sup, sup_start, total_days){
  if(total_days < sup_start + D_sup){
    total_days = sup_start + D_sup
  }
  
  # let disease spread in population initially
  first_wave = simulate_SEIR_differential(N, R0, D, D_pre, 
                                          days = sup_start, 
                                          S_init = N-1, 
                                          E_init = 1, 
                                          I_init = 0, 
                                          R_init = 0, 
                                          Total_I_init = 0)
  

  # suppress disease after sup_start days for sup_dur days with R0 falling immediately to R0_sup
  initial_values = tail(first_wave, 1)
  suppression_period = simulate_SEIR_differential(N, R0_sup, D, D_pre, 
                                                 days = D_sup, 
                                                 S_init = initial_values[["S"]], 
                                                 E_init = initial_values[["E"]], 
                                                 I_init = initial_values[["I"]], 
                                                 R_init = initial_values[["R"]], 
                                                 Total_I_init = initial_values[["Total_I"]])
  suppression_period = suppression_period %>% mutate(time=time+sup_start)
  
  # release suppression and allow second wave where original R0 applies immediately after sup_start + D_sup days
  initial_values = tail(suppression_period, 1)
  second_wave = simulate_SEIR_differential(N, R0, D, D_pre, 
                                           days = total_days-(sup_start + D_sup), 
                                           S_init = initial_values[["S"]], 
                                           E_init = initial_values[["E"]], 
                                           I_init = initial_values[["I"]], 
                                           R_init = initial_values[["R"]], 
                                           Total_I_init = initial_values[["Total_I"]])
  second_wave = second_wave %>% mutate(time=time+sup_start+D_sup)
  
  markov_trace = suppressWarnings(bind_rows(first_wave, suppression_period, second_wave) %>% distinct())
  
  return(markov_trace)
}

plot_SEIR = function(markov_trace, N, R0, D, D_pre){
  HIT = calculate_HIT(R0)
  HIT_date = calculate_HIT_date(markov_trace, R0, N)
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
         subtitle = bquote("Input parameters:" ~ R[0] ~ " = " ~ .(R0) ~ ", D' = " ~ .(D_pre) ~ ", D = " ~ .(D) ~ ", N = " ~ .(format(N, big.mark=",")) ),
         caption = paste0("Herd immunity threshold of ", round(100*HIT) ,"%",ifelse(is.na(HIT_date)," not achieved in this simulation", paste0(" achieved after ",HIT_date," days")))) +
    scale_colour_manual(name="States: ", values=c("red","orange","blue","green")) +
    scale_linetype_manual(name="States: ", values=c(1,2,1,1)) +
    theme_bw(base_size = 10) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  # return SEIR plot
  return(seir_plot)
}

plot_infection_curves = function(markov_trace_1, markov_trace_2, label_1, label_2){
  total_infections_1 = round(max(markov_trace_1$Total_I))
  infections_1 = markov_trace_1 %>% select(time, I) %>% mutate(state=label_1)
  peak_infections_1 = round(max(infections_1$I))
  total_infections_2 = round(max(markov_trace_2$Total_I))
  infections_2 = markov_trace_2 %>% select(time, I) %>% mutate(state=label_2)
  peak_infections_2 = round(max(infections_2$I))
  
  graph_data = suppressWarnings(bind_rows(infections_1, infections_2)) %>%
    mutate(state = factor(state,c(label_1,label_2),c(label_1,label_2))) 
  
  # plot the infection curve graph
  infection_plot = ggplot(graph_data, aes(x=time, y=I, group=state)) +
    xlab("Time (days since first case)") +
    ylab("Number of people") +
    geom_line(aes(colour=state, linetype=state), size=1) +
    labs(title = paste0("Infected over time"), 
         subtitle = paste0(label_1," total infections: ",format(total_infections_1, big.mark=","),"; ",label_2," total infections: ", format(total_infections_2, big.mark=",")),
         caption = paste0(label_1," peak infections: ",format(peak_infections_1, big.mark=","),"; ",label_2," peak infections: ", format(peak_infections_2, big.mark=","))) +
    scale_colour_manual(name="States", values=c("blue","blue")) +
    scale_linetype_manual(name="States", values=c(1,2)) +
    theme_bw(base_size = 10) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  return(infection_plot)
}

