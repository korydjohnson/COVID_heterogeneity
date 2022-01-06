##################################3
# Generating Samples from Nu; Separate File as a "one-off"
##################################3

library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
library(truncnorm)
options(dplyr.summarise.inform = FALSE)

source("scripts/simulation_final.R")
source(paste0("scripts/szenarios_list_Omega1.r"))
source("scripts/plotting_functions.r")
source("scripts/processing_functions_final.R")

# input -------------------------------------------------------------------

strain_names = c("W", "A", "B", "G", "D", "O", "N")
data_folder = "samples/"
plot_folder = "plots/"
nRuns = 7
sim_end = as.Date("2022-05-01")
toWane = toBoost = T
q = c(0.025, 0.25, 0.5, 0.75, 0.975)
nu = 13
history=T


# script ------------------------------------------------------------------

vaccination_data_cutoff = as.Date("2021-08-08")
maximum_vaccine_uptake = 0.85 # percentage of whole population which will get a vaccine

files = list.files(data_folder)
if (sum(str_detect(files, "Nu3")) == 0) {
  kPrev=0
} else {
  kPrev = files %>% 
    str_subset("Nu3") %>% 
    str_extract("_s.*") %>% 
    str_extract("[:digit:]+") %>% 
    as.integer() %>% max()
}

### RESISTANCE PARAMETERS
vaccine_names = c("1", "2", "3")
vaccine_names_full = c("mRNA", "vector", "booster")
wane_categories = c("V","O")

resistance_ranges = 
  read.csv("data/resistance_parameters.csv", stringsAsFactors = F) %>%
  #filter(bound == "middle") %>% select(-bound) %>%
  select(c("name", "type", "bound", all_of(strain_names))) %>%
  filter((type != "infection") | (name %in% strain_names))

wane_ranges = 
  read.csv("data/resistance_parameters.csv", stringsAsFactors = F) %>%
  select(c(name, type, bound, V=Vm6, O=Om6)) %>%
  filter((type != "infection") | (name %in% strain_names))

forecasting_date = szenarios_list[[1]]$forecasting_date
simLength = as.integer(sim_end - forecasting_date)


R_observed = szenarios_list[[1]]$R_observed

geometric_growth_model = szenarios_list[[1]]$geometric_growth_model
vaccination_plan_with_past_data = 
  szenarios_list[[1]]$vaccination_plan_with_past_data

delta_prev = szenarios_list[[1]]$delta_prev
import = szenarios_list[[1]]$import
import$O = 0
import$N = 1

I_import_start = (as.Date("2021-12-01") - as.Date(szenarios_list[[1]]$forecasting_date)) %>% 
  as.numeric("days")
source("scripts/multi_variant_population_parameters_imports.r")

nuCores = detectCores() - 1
for (index in 1:2) {
  pars_controls = controls_list[[index]]
  
  # run with standard generation interval -----------------------------------
  # fileName = paste(szenarios_list[[1]]$name, "Nu46", 
  #                  pars_controls$name, sep="_")
  # 
  # registerDoParallel(nuCores)  # set up cluster and register
  # foreach(k = seq_len(nRuns)) %dopar% {
  #   source("scripts/gather_data.r")
  #   pars_compartments$strains$pars[[7]]$I_import_start = I_import_start
  #   # "omega" is actually "nu"; this is done for plotting
  #   pars_compartments$strains$pars[[7]]$name = "O"
  #   pars_compartments$strains$pars[[6]] = pars_compartments$strains$pars[[7]]
  #   pars_compartments$strains$pars[[7]] = NULL
  #   pars_compartments$strains$names = pars_compartments$strains$names[-7] 
  #   
  #   variableName = paste0(fileName, paste0("_s", kPrev + k))
  #   file_name = paste0(
  #     data_folder,
  #     variableName,
  #     ".Rdata"
  #   )
  #   assign(variableName,
  #          sim_multistrain(
  #            parsList = pars_compartments,
  #            pars_controls = pars_controls,
  #            simLength = simLength, wane=toWane, history=history
  #          )
  #   )
  #   save(list = variableName, file=file_name)
  #   rm(list=variableName)  # all names are unique, so remove
  # }
  # stopImplicitCluster()
  
  
  # run with shorter generation interval ------------------------------------
  fileName = paste(szenarios_list[[1]]$name, "Nu3", 
                   pars_controls$name, sep="_")
  
  registerDoParallel(nuCores)  # set up cluster and register
  foreach(k = seq_len(nRuns)) %dopar% {
    source("scripts/gather_data.r")
    pars_compartments$strains$pars[[7]]$I_import_start = I_import_start
    pars_compartments$strains$pars[[7]]$w = 
      generationInt(iTime = 13, E = E_gT[["N"]], Var = 2.63 ^ 2)
    # "omega" is actually "nu"; this is done for plotting
    pars_compartments$strains$pars[[7]]$name = "O"
    pars_compartments$strains$pars[[6]] = pars_compartments$strains$pars[[7]]
    pars_compartments$strains$pars[[7]] = NULL
    pars_compartments$strains$names = pars_compartments$strains$names[-7] 
    
    variableName = paste0(fileName, paste0("_s", kPrev + k))
    file_name = paste0(
      data_folder,
      variableName,
      ".Rdata"
    )
    assign(variableName,
           sim_multistrain(
             parsList = pars_compartments,
             pars_controls = pars_controls,
             simLength = simLength, wane=toWane, history=history
           )
    )
    save(list = variableName, file=file_name)
    rm(list=variableName)  # all names are unique, so remove
  }
  stopImplicitCluster()  
}


# plotting ----------------------------------------------------------------

data_list = list(
  list(
    file_left = list(date="2021-08-08", control="proactive", bounds="25-150", variant="Nu46"),
    file_right = list(date="2021-08-08", control="reactive", bounds="25-150", variant="Nu46"),
    forecasting_date = as.Date("2021-08-08"),
    lower = 318 / 2,
    upper = 318 * 3,
    upper_R = 1.2,
    variants_to_plot = c("Delta", "Omega"),
    plot_file = paste0(plot_folder, "Nu46_25-150.png"),
    info_left = "proactive",
    info_right = "reactive",
    plot_cutoff = 50
  ),
  list(
    file_left = list(date="2021-08-08", control="proactive", bounds="25-150", variant="Nu3"),
    file_right = list(date="2021-08-08", control="reactive", bounds="25-150", variant="Nu3"),
    forecasting_date = as.Date("2021-08-08"),
    lower = 318 / 2,
    upper = 318 * 3,
    upper_R = 1.2,
    variants_to_plot = c("Delta", "Omega"),
    plot_file = paste0(plot_folder, "Nu3_25-150.png"),
    info_left = "proactive",
    info_right = "reactive",
    plot_cutoff = 100
  ),  
  list(
    file_left = list(date="2021-08-08", control="proactive", bounds="25-50", variant="Nu46"),
    file_right = list(date="2021-08-08", control="reactive", bounds="25-50", variant="Nu46"),
    forecasting_date = as.Date("2021-08-08"),
    lower = 318 / 2,
    upper = 318,
    upper_R = 1.2,
    variants_to_plot = c("Delta", "Omega"),
    plot_file = paste0(plot_folder, "Nu46_25-50.png"),
    info_left = "proactive",
    info_right = "reactive",
    plot_cutoff = 30
  ),
  list(
    file_left = list(date="2021-08-08", control="proactive", bounds="25-50", variant="Nu3"),
    file_right = list(date="2021-08-08", control="reactive", bounds="25-50", variant="Nu3"),
    forecasting_date = as.Date("2021-08-08"),
    lower = 318 / 2,
    upper = 318,
    upper_R = 1.2,
    variants_to_plot = c("Delta", "Omega"),
    plot_file = paste0(plot_folder, "Nu3_25-50.png"),
    info_left = "proactive",
    info_right = "reactive",
    plot_cutoff = 30
  )
)

cases = read.csv("data/cases_who.csv") %>% 
  mutate(date = as.Date(date)) %>% 
  rename(cases_total = cases) %>% 
  mutate(cases_ra = round(stats::filter(cases_total, rep(1 / 7, 7), sides = 1))) %>% 
  filter(date <= "2021-10-31")

walk(data_list, make_rho_plot)
plots = map(data_list, make_comparison_plots, toReturn=T)
walk(data_list, make_side_by_side_plot, cases, data_folder, as.Date("2022-05-01"))

# plt1 = plots[[1]]$mitPlot
# plt2 = plots[[2]]$mitPlot
# mergePlots(plt1, plt2, captions = c("test1", "test2"))








