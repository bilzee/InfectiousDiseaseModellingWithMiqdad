library("tidyverse")
library("lubridate")
library("RcppRoll")


get_country_list = function(){
  data = read_csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv")
  data = data %>% select(countriesAndTerritories, geoId) %>% distinct() %>% arrange(countriesAndTerritories)
  countries = as.list(data$geoId)
  names(countries) = gsub("_", " ", data$countriesAndTerritories)
  return(countries)
}

#####################################################################
# This function downloads data on all COVID-19 cases and deaths 
# from the ECDC, smooths these data by calculating moving averages
# and discards data prior to a set number of deaths being observed
# in the moving average
#
# moving_average_days = the number of days over which to calculate
# the moving average e.g. 7
# deaths_cut_off = the number of deaths after which the epidemic is
# established in the country, data prior to this is discarded e.g. 3
#####################################################################
get_tidy_data_set = function(moving_average_days, deaths_cut_off){
  # get the latest COVID-19 data from the ECDC
  data = read_csv("https://opendata.ecdc.europa.eu/covid19/casedistribution/csv")
  
  # local function to calculate the start date of the analysis per country
  get_start_date = function(data, cut_off){
    data %>% filter(deaths>cut_off) %>% summarise(min(dateRep, na.rm=TRUE)) %>% pull()
  }
  
  tidy_data = data %>%
    # convert the dateRep to a date datatype
    mutate(dateRep = dmy(dateRep)) %>%
    # for each country
    group_by(geoId) %>%
    # arrange by date
    arrange(dateRep) %>%
    # convert deaths and cases into x day moving average to deal with reporting blips
    mutate(deaths = roll_meanr(deaths, moving_average_days), cases = roll_meanr(cases, moving_average_days)) %>%
    # start from when daily deaths pass 3
    filter(dateRep >= get_start_date(tibble(dateRep, deaths), deaths_cut_off)) %>%
    # calculate days from this start date
    mutate(days = as.numeric(dateRep-min(dateRep, na.rm=TRUE))) %>%
    # now that we have finished doing analysis by country ungroup
    ungroup() %>%
    # select the data fileds that we want to return
    select(dateRep, days, geoId, deaths, cases) %>%
    # convert our data into tidy format for further analysis
    gather(variable, value, deaths, cases)
  
  return(tidy_data)
}

#####################################################################
# This function plots moving average COVID-19 cases and deaths 
# fon a log scale against days from the first time a given number of 
# deaths (moving average) are observed 
#
# moving_average_days = the number of days over which to calculate
# the moving average e.g. 7
# deaths_cut_off = the number of deaths after which the epidemic is
# established in the country, data prior to this is discarded
# countries = the geoIds (2 character country codes) of the countries 
# that you want to plot data for 
#####################################################################
plot_cases_and_deaths = function(moving_average_days, deaths_cut_off, countries){
  tidy_data = suppressWarnings(get_tidy_data_set(moving_average_days, deaths_cut_off))
  
  plot = ggplot(tidy_data %>% filter(geoId %in% countries), 
                aes(x=days, y=log(value), group=geoId, colour=geoId)) + 
    geom_point(size=0.5) +
    geom_smooth(size=1) +
    xlab(paste0("days from ",deaths_cut_off," deaths")) +
    ylab(paste0(moving_average_days," day moving average (log scale)")) +
    facet_wrap(.~variable, scales="free") +
    labs(
      title = "COVID-19 deaths and cases over time",
      subtitle = "Measured in logs",
      caption = "Source ECDC: https://opendata.ecdc.europa.eu/covid19/casedistribution/csv"
    ) +
    theme_bw(base_size = 15) + 
    theme(plot.margin = unit(c(1, 1, 1, 1), "lines"),
          legend.position = "bottom")
  
  return(plot)
}
