###########################################
# Get Samples from Simulation: master
###########################################


# input -------------------------------------------------------------------
settingVec = c("Delta", "Omega1", "Omega2", "DeltaWane", "Omega1Wane",
               "DeltaWaneBooster", "Omega1WaneBooster")
strain_names = c("W", "A", "B", "G", "D", "O")

if (!exists("setting")) {
  index = 1  # 1-5 to select setting from vector above
} else {
  index = which(settingVec==setting)
  if(length(index) == 0) {
    stop("Setting in:'Delta', 'Omega1', 'Omega2', 'DeltaWane', 'Omega1Wane',
         'DeltaWaneBooster', 'Omega1WaneBooster'.")
  }
}

data_folder = "samples/"

if (!exists("nRuns")) {
  nRuns = 7
}
addToPrevious = F  # set to true to not replace previous files, just add

# either
sim_end = as.Date("2022-05-01")
simLength_with_end = T
# or
# simLength_with_end = F
# simLength = 30

# code --------------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
library(truncnorm)
options(dplyr.summarise.inform = FALSE)

source("scripts/simulation_final.R")
# source("scripts/plotting_functions.r")
source(paste0("scripts/szenarios_list_", 
              str_extract(settingVec[index], "[^W]*"), ".r"))

ifelse(!dir.exists(file.path(data_folder)), dir.create(file.path(data_folder)), FALSE)

run_safely = ifelse(nRuns < 5, T, F)
history = !run_safely  # if false, outputs entire object (compartments etc)
# whoa; this is actually a promise of sorts, or at least stored as a function
# it overwrites the "wane" function defined previously
toWane = ifelse(str_detect(settingVec[index], "Wane"), T, F)
toBoost = ifelse(str_detect(settingVec[index], "Booster"), T, F)  # used in multi_variant_population_parameters.r
if (addToPrevious) {
  files = list.files(data_folder)
  if (sum(str_detect(files, settingVec[index])) == 0) {
    kPrev=0
  } else {
    kPrev = files %>% 
      str_subset(settingVec[index]) %>% 
      str_extract("_s.*") %>% 
      str_extract("[:digit:]") %>% 
      as.integer() %>% max()
  }
} else {
  kPrev = 0
}
q = c(0.025, 0.25, 0.5, 0.75, 0.975)
nu = 13

vaccination_data_cutoff = as.Date("2021-08-08")
maximum_vaccine_uptake = 0.85 # percentage of whole population which will get a vaccine

### RESISTANCE PARAMETERS
vaccine_names = c("1", "2")
vaccine_names_full = c("mRNA", "vector")
wane_categories = c("V","O")
if (toWane) {
  vaccine_names = c(vaccine_names, "3")
  vaccine_names_full = c("mRNA", "vector", "booster")
}

resistance_ranges = 
  read.csv("data/resistance_parameters.csv", stringsAsFactors = F) %>%
  #filter(bound == "middle") %>% select(-bound) %>%
  select(c("name", "type", "bound", all_of(strain_names))) %>%
  filter((type != "infection") | (name %in% strain_names))

wane_ranges = 
  read.csv("data/resistance_parameters.csv", stringsAsFactors = F) %>%
  #filter(bound == "middle") %>% select(-bound) %>%
  select(c(name, type, bound, V=Vm6, O=Om6)) %>%
  filter((type != "infection") | (name %in% strain_names))

for (j in 1:length(szenarios_list)) {
  print(paste0("Do: ", szenarios_list[[j]]$name))
  forecasting_date = szenarios_list[[j]]$forecasting_date
  if (simLength_with_end) {
    simLength = as.integer(sim_end - forecasting_date)
  }
  
  R_observed = szenarios_list[[j]]$R_observed
  
  geometric_growth_model = szenarios_list[[j]]$geometric_growth_model
  vaccination_plan_with_past_data = 
    szenarios_list[[j]]$vaccination_plan_with_past_data
  
  delta_prev = szenarios_list[[j]]$delta_prev
  import = szenarios_list[[j]]$import
  source("scripts/multi_variant_population_parameters_imports.r")
  ### DO SIMULATIONS
  for (i in 1:length(controls_list)) {
    pars_controls = controls_list[[i]]
    print(paste0("Do: ",szenarios_list[[j]]$name, ", ",
        settingVec[index], ", ", pars_controls$name))
    fileName = paste(szenarios_list[[j]]$name, settingVec[index], 
                     pars_controls$name, sep="_")

    tic = Sys.time()
    print(tic)
    print("Simulate...")
    
    if (!run_safely) {
      nuCores = detectCores() - 1
      registerDoParallel(nuCores)  # set up cluster and register
      foreach(k = seq_len(nRuns)) %dopar% {
        source("scripts/gather_data.r")
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
    } else {
      print("..safely..")
      sim_multistrain_safely = safely(sim_multistrain)
      nuCores = detectCores() - 1
      registerDoParallel(nuCores)  # set up cluster and register
      simData_safely = foreach(k = seq_len(nRuns)) %dopar% {
        #1.draw sample of the strain resistances
        #resistance_data = sample_variant_data_truncnormal(resistance_ranges, variant_names = strain_names)
        source("scripts/gather_data.r")
        sim_multistrain_safely(
          parsList = pars_compartments,
          pars_controls = pars_controls,
          simLength = simLength, wane=toWane, history=history
        )
      }
      stopImplicitCluster()
      
      sim_okay = sapply(map(simData_safely, 2), is_null)
      simData = map(simData_safely[sim_okay], 1)
      error_messages = map(simData_safely[!sim_okay], 2)
      
      #save errorlog
      error_file =  paste0(
        data_folder,
        fileName,
        "-errorlog.txt"
      )
      cat(paste0("There were ", length(error_messages), " errors.", "\n\n"), file = error_file)
      cat(
        error_messages %>% sapply(toString),
        file = error_file,
        sep = "",
        append = T
      )
      
      file_name = paste0(
        data_folder,
        fileName,
        ".Rdata"
      )
      
      print(paste0("Save:" , file_name))
      save(simData, file = file_name)
    }
    print(paste0("...done. total runtime: ", Sys.time() - tic))
    
  }
}
