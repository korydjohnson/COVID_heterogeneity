###########################################3
# Create Graphs: master
###########################################3

library(tidyverse)
library(magrittr)
library(foreach)
options(dplyr.summarise.inform = FALSE)
library(gridExtra)
library(grid)
library(ggthemes)

Sys.setlocale("LC_TIME", 'en_GB.UTF-8')

source("scripts/simulation_final.R")
source("scripts/processing_functions_final.R")
source("scripts/plotting_functions.r")

make_plots = function(option="Delta", toReturn=F) {
  
  nu = 13
  data_folder = "samples/"
  plot_folder = "plots/"
  ifelse(!dir.exists(file.path(data_folder)), dir.create(file.path(data_folder)), FALSE)
  ifelse(!dir.exists(file.path(plot_folder)), dir.create(file.path(plot_folder)), FALSE)
  
  # "Omega1_delay",
  optionList = c("Delta", "Omega1", "Omega2_delay", "DeltaWane", "Omega2Wane_delay", 
                 "Omega1Wane", "DeltaWaneBooster", "Omega1WaneBooster")
  optionList = c(optionList, paste0(optionList, "_rhoV"), paste0(optionList,"_prevalence"),
                 paste0(optionList, "_compare"))
  if (! (option %in% optionList)) {
    stop(cat("Option must be in:", optionList))
  }
  
  
  # Delta -------------------------------------------------------------------
  
  data_list_Delta = list(
    "low" = list(
      file_left = list(date="2021-06-12", control="proactive", bounds="25-50", variant="Delta"),
      file_right = list(date="2021-06-12", control="reactive", bounds="25-50", variant="Delta"),
      forecasting_date = as.Date("2021-06-12"),
      lower = 318 / 2,
      upper = 318,
      upper_R = 1.2,
      variants_to_plot = c("Alpha", "Beta", "Gamma", "Delta"),
      plot_file = paste0(plot_folder, "Delta_25-50.png"),
      info_left = "proactive",
      info_right = "reactive",
      plot_cutoff = 20
    ),
    "high" = list(
      file_left = list(date="2021-06-12", control="proactive", bounds="25-150", variant="Delta"),
      file_right = list(date="2021-06-12", control="reactive", bounds="25-150", variant="Delta"),
      forecasting_date = as.Date("2021-06-12"),
      lower = 318 / 2,
      upper = 318 * 3,
      upper_R = 1.2,
      variants_to_plot = c("Alpha", "Beta", "Gamma", "Delta"),
      plot_file = paste0(plot_folder, "Delta_25-150.png"),
      info_left = "proactive",
      info_right = "reactive",
      plot_cutoff = 50
    )
  )
  
  # Omega -------------------------------------------------------------------
  
  data_list_Omega = list(
    "high" = list(
      file_left = list(date="2021-08-08", control="proactive", bounds="25-150", variant=""),
      file_right = list(date="2021-08-08", control="reactive", bounds="25-150", variant=""),
      forecasting_date = as.Date("2021-08-08"),
      lower = 318 / 2,
      upper = 318 * 3,
      upper_R = 1.2,
      variants_to_plot = c("Delta", "Omega"),
      plot_file = paste0(plot_folder, "Omega_25-150.png"),
      info_left = "proactive",
      info_right = "reactive",
      plot_cutoff = 50
    ),
    "low" = list(
      file_left = list(date="2021-08-08", control="proactive", bounds="25-50", variant=""),
      file_right = list(date="2021-08-08", control="reactive", bounds="25-50", variant=""),
      forecasting_date = as.Date("2021-08-08"),
      lower = 318 / 2,
      upper = 318,
      upper_R = 1.2,
      variants_to_plot = c("Delta", "Omega"),
      plot_file = paste0(plot_folder, "Omega_25-50.png"),
      info_left = "proactive",
      info_right = "reactive",
      plot_cutoff = 50
    )
  )
  
  # Cases -------------------------------------------------------------------
  
  realData_cutoff = ifelse(str_detect(option, "Delta"), "2021-09-01", "2021-08-08")
  # realData_cutoff = "2021-0"
  
  # cases_who = read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv") %>%
  cases = read.csv("data/cases_who.csv") %>% 
    mutate(date = as.Date(date)) %>% 
    rename(cases_total = cases) %>% 
    mutate(cases_ra = round(stats::filter(cases_total, rep(1 / 7, 7), sides = 1))) %>% 
    filter(date <= "2021-10-31")
  # %>% 
  #   filter(date <= realData_cutoff)
  
  
  # Plots -------------------------------------------------------------------
  
  variant = str_extract(option, "^[^_[:digit:]W]*")
  data_list = get(paste0("data_list_", variant))
  
  simName_variant = str_extract(option, "^[^_]*")
  for (i in 1:length(data_list)) {
    data_list[[i]]$plot_file = str_replace(data_list[[i]]$plot_file, variant, option)
    data_list[[i]]$file_left$variant = simName_variant
    data_list[[i]]$file_right$variant = simName_variant
  }
  
  if (str_detect(option, "delay")) {
    stopifnot(str_detect(simName_variant, "Omega2"))
    for (i in 1:length(data_list)) {
      data_list[[i]]$file_left$delay = 21
      data_list[[i]]$file_right$delay = 21
      data_list[[i]]$plot_cutoff = 60
    }
  }
  
  if(str_detect(option, "prevalence")) {
    walk(data_list, make_prevalence_plots, data_folder)
  } else if(str_detect(option, "rhoV")) {
    plots = map(data_list, make_rho_plot, toReturn=toReturn)
  } else if(str_detect(option, "compare")) {
    plots = map(data_list, make_comparison_plots, toReturn=toReturn)
    # walk(data_list, make_comparison_plots)
  } else {
    walk(data_list, make_side_by_side_plot, cases, data_folder, as.Date("2022-05-01"))
  }
  if (toReturn) {
    plots
  }
}
