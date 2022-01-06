######################################################
# Script for Multistrain Simulation
######################################################

require(magrittr, tiyverse)

# Helper Functions Init ---------------------------------------------------

seasonality_day = function(date, e = 0.4, t_peak = as.Date("2021-01-01")) {
  time_to_peak = -as.numeric(date - t_peak)
  (1 - e * (cos(2 * pi * (time_to_peak / 365.25 + 1 / 2)) + 1) / 2)
}

# compute_rolling_averages <- function(I, ra = 7) {
#   #Given a time series "I" this function computes the rolling averages over "ra" many days.
#   #The first (ra-1) days the rolling averages will be computed over resp. shorter periods.
#   #OUTPUT: RA is a vector of same length as I.
#   RA <- rep(NA, length(I))
#   for (t in 1:(ra - 1)) { RA[t] <- sum(I[1:t]) / t }
#   for (t in ra:length(I)) { RA[t] <- sum(I[(t - ra + 1):t]) / ra }
#   
#   RA
# }

#Probability generating function of the serial interval
G_w = function(x, w = generationInt()) {
  nu = length(w)
  return(sum(w * x ^ (1:nu)))
}

#(Numerical) inverse of probability generating function of the serial interval
G_w_inv = function(x, w = generationInt()) {
  uniroot(
    f = function(b) {
      x - G_w(x = b, w)
    },
    interval = c(0.001, 10)
  )$root
}

#Estimation of growth rate given R(_eff)
estimate_beta = function(R = 1, w = generationInt()) {
  return(1 / G_w_inv(x = 1 / R, w = w))
}

estimate_growth_rate = function(g, R) {
  w = g$w
  return(1 / G_w_inv(x = 1 / R, w = w))
}

compute_initial_mit = function(pars_compartments,
                               R_observed = 1,
                               L = 1) {
  N = pars_compartments$population$N
  strain_names = pars_compartments$strains$names
  
  groups = c(pars_compartments$strains$pars,
             #strain groups
             pars_compartments$vaccines$pars,
             #vaccination groups
             # add uninfected group
             list(
               list(
                 S0 = pars_compartments$population$pSuscept * N,
                 res = setNames(as.list(rep(
                   0, length(strain_names)
                 )), strain_names)
               )
             ))
  W = with(pars_compartments$population, sum(rev(w) * I0))
  
  group_sizes = map_dbl(groups, ~.$S0)
  
  #compute Gamma_g and Lambda_g, only relevant for infectious strain!
  Gamma_g = sapply(strain_names, function(strain) {
    sum((1 - map_dbl(groups, function(g) {
      g$res[[strain]]
    })) * group_sizes) / N
  })
  
  Lambda_g = sapply(pars_compartments$strains$pars,
                    function(g) {
                      w = g$w
                      I0 = g$I0
                      return(g$R * sum(I0 * rev(w)))
                    })
  1 - R_observed * W / (L * sum(Gamma_g * Lambda_g))
}

get_resistances = function(resistance_data, group_names = c("mRNA", "vector")) {
  strain_names = resistance_data %>% 
    filter(type == "infection") %>% pull(name) %>% unique
  resistances = foreach(group = group_names) %do% {
    resistance_data %>% filter(name == group) %>% 
      select(all_of(strain_names)) %>% as.list
  }
  setNames(c(list(setNames(
    as.list(rep(0, length(strain_names))), strain_names
  )),
  resistances), c("_", group_names))
}

# Sample truncnorm regions
sample_variant_data_truncnormal = function(resistance_ranges, variant_names) {
  foreach(name_ = unique(resistance_ranges$name), .combine = rbind) %do% {
    upper  = resistance_ranges %>% 
      filter(name == name_, bound == "upper") %>% 
      select(all_of(variant_names)) %>% unlist
    middle = resistance_ranges %>% 
      filter(name == name_, bound == "middle") %>% 
      select(all_of(variant_names)) %>% unlist
    lower  = resistance_ranges %>% 
      filter(name == name_, bound == "lower") %>% 
      select(all_of(variant_names)) %>% unlist
    
    sd = (upper - lower) / 4  # approx 95% of samples within 2sd of mean
    
    lower[sd == 0] = -Inf
    upper[sd == 0] = Inf
    
    resistance_ranges %>%
      filter(name == name_) %>%
      select(name, type) %>%
      head(1) %>%
      cbind(as_tibble(setNames(as.list(
        rtruncnorm(
          length(variant_names),
          a = lower,
          b = upper,
          mean = middle,
          sd = sd
        )
      ), variant_names)))
  }
}

read_samples = function(folder="samples", date = "2021-06-12", variant="Delta", 
                        control="reactive", bounds = "25-50", delay=7) {
  files = list.files(folder)
  indexes = enframe(files) %>% 
    separate(value, sep="_", into = c("date", "variant_", "control_", "bounds_",
                                      "delay_", NA)) %>% 
    filter(variant_ == variant, 
           control_ == control,
           bounds_ == bounds,
           delay_ == paste0("delay", delay)) %>% 
    pull(name)
  
  env = new.env()
  walk(paste(folder, files[indexes], sep="/"), load, envir = env)
  as.list(env)
}

# Helper Functions, Simulation --------------------------------------------

seasonality = function(pars_population, tPeak = as.Date("2021-01-01")) {
  dateRange = with(as.list(pars_population),
                   (date - iTime + 1):(date + simLength) - as.numeric(tPeak))
  timeToPeak = -as.numeric(dateRange)
  with(pars_population,
       1 - e_seasonality*(cos(2*pi*(timeToPeak/365.25 + 1/2)) + 1) / 2)
}

generationInt = function(iTime = 13, E = 4.46, Var = 2.63 ^ 2) {
  w = diff(pgamma(0:iTime, E ^ 2 / Var, E / Var))
  w / sum(w)
}

get_strainNames = function(compartments) {
  names(compartments) %>% 
    .[which(str_detect(., "[:upper:]$"))]
}

# Constructor Functions ---------------------------------------------------

make_strain = function(name, strainPars, pars_population, simLength, day = 1) {
  seasonality = pars_population$seasonality
  N = pars_population$N
  mitig = NA
  infName = str_sub(name, start = -1)
  w = strainPars$w
  theta0 = with(strainPars, rgamma(length(I0), shape = I0 * k, rate = k / R))
  theta = c(theta0, rep(0, simLength))
  I = with(strainPars, c(I0, rep(0, simLength)))
  
  if (!is.null(strainPars$I_import)) {  # imports
    if(is.null(strainPars$I_import_start)){
      I_i = with(strainPars, 
                 c(rep(0, length(I0)), rep(I_import, length.out = simLength)))
    } else {
      I_i = with(strainPars, 
                 c(rep(0, length(I0) + I_import_start-1), 
                   rep(I_import, length.out = simLength-I_import_start+1)))
    }
  } else{  # no imports
    I_i = with(strainPars, rep(0, simLength + length(I0)))
  }
  theta_i = with(strainPars, rgamma(length(I_i), shape = I_i * k, rate = k / R))
  S = with(strainPars, c(rep(S0, iTime), rep(0, simLength)))
  
  with(strainPars, stopifnot(length(I0) == iTime)) # C needs to be correct length
  C = with(strainPars, c(S0 + cumsum(I0), rep(0, simLength)))
  t = length(strainPars$I0) + day  # next time period for simulation
  S[t] = S[t - 1]  # remove susceptibles "during the day" as happens repeatedly
  C[t] = C[t - 1]  # remove susceptibles "during the day" as happens repeatedly
  W = GammaTilde = Theta = Lambda = NA
  Gamma_Tilde_vec = R_e_t = rep(NA, simLength) # group R_eff
  
  list(
    state = function() {
      list(
        name = name, infName = infName, day = t, theta = theta[t], I = I[t],
        S = S[t], C = C[t], mitig = mitig, W = W, GammaTilde = GammaTilde,
        Theta = Theta, Lambda = Lambda
      )
    },
    history = function() {
      list(theta = theta, I = I, S = S, C = C, I_i = I_i,
           R_e_t = R_e_t, Gamma_Tilde_vec = Gamma_Tilde_vec)
    },
    pars = function() {
      strainPars
    },
    create_new_infs = function(susceptible) {
      resist = susceptible$pars()$res[[infName]]
      relHist = (t - strainPars$iTime):(t - 1)
      INu = rpois(1, (susceptible$state()$S / N) * (1 - resist) * (1 - mitig) *
                    seasonality[t] * crossprod(theta[rev(relHist)], w)) %>% 
        min(susceptible$state()$S)  # can't infect more than is in group
      thetaNu = with(strainPars, rgamma(1, shape = INu * k, rate = k / R))
      return(list(I = INu, theta = thetaNu))
    },
    add_new_infs = function(newInfs) {
      I[t] <<- I[t] + newInfs$I
      theta[t] <<- theta[t] + newInfs$theta
    },
    next_day = function() {
      S[t] <<- S[t] + I[t - strainPars$iTime]  # recovered become susceptible
      C[t] <<- C[t] + I[t]
      theta[t] <<- theta[t] + theta_i[t]
      t <<- t + 1
      S[t] <<- S[t - 1]  # remove susceptibles (repeatedly) "during the day"
      C[t] <<- C[t-1]
    },
    compute_summary = function(compartments, yesterday=F) {
      # rev(.) argument goes to t to include today; new_day() hasn't been called
      if (yesterday) { t <<- t - 1 }
      W <<- crossprod(I[rev((t - strainPars$iTime + 1):t)], w)
      GammaTilde <<- N^{-1} * 
        sum(map_dbl(compartments, ~ .$state()$S * (1-.$pars()$res[[infName]])))
      Gamma_Tilde_vec[t] <<- GammaTilde
      Theta <<- crossprod(theta[rev((t - strainPars$iTime + 1):t)], w)
      Lambda <<- strainPars$R * W
      R_e_t[t] <<- strainPars$R * seasonality[t] * (1 - mitig) * GammaTilde
      if (yesterday) { t <<- t + 1 }
    },
    remove_susceptibles = function(toRemove) {
      # called when moving groups, ie reinfected
      S[t] <<- S[t] - toRemove
      C[t] <<- C[t] - toRemove
    },
    add_susceptibles = function(toAdd) {
      # called when vaccinate for a strain/vaccine interaction
      S[t] <<- S[t] + toAdd
      C[t] <<- C[t] + toAdd
    },
    update_mitigation = function(mit) {
      mitig <<- mit
    }
  )
}

add_compartment = function(compartments, baseName, interventionName,
                           interactionName, pars_population, simLength, day, 
                           type=NULL) {
  stopifnot(type %in% c("vaccine", "variant", "wane"))
  if (type == "variant") {
    pars = compartments[[interventionName]]$pars()
    pars$I_import = NULL; pars$I_import_start = NULL  # no imports
    pars$I0 = 0*pars$I0; pars$S0 = 0  # new groups lack previous cases
  } else {
    pars = compartments[[baseName]]$pars()
    pars$S0 = 0; pars$Vschedule = rep(0, simLength)  # in case lacking base
  }
  
  pars$res = switch(
    type,
    vaccine = map2(compartments[[baseName]]$pars()$res, 
                   compartments[[interventionName]]$pars()$res, ~ max(.x, .y)),
    variant = map2(compartments[[baseName]]$pars()$res, 
                   compartments[[interventionName]]$pars()$res, ~ 1 - (1-.x)*(1-.y)),
    wane = compartments[[interventionName]]$pars()$omega
    # waneVaccine = compartments[[interventionName]]$pars()$res
    # waneVaccine = map2(compartments[[interventionName]]$pars()$res,
    #                    compartments[[baseName]]$pars()$res,  # flipped as base/intervention flipped 
    #                    ~1-(1-.x)*(1-compartments[[interventionName]]$pars()$omega*.y)/.y)
  )
  
  if (type == "variant") {
    compartments[[interactionName]] = 
      make_strain(interactionName, pars, pars_population, simLength, day)
    compartments[[interventionName]]$state()$mitig %>%  # mitig begins as NA 
      compartments[[interactionName]]$update_mitigation()
  } else {
    compartments[[interactionName]] = 
      make_vaccine(interactionName, pars, simLength, day)
  }
  
  compartments
}

make_vaccine = function(name, vaccinePars, simLength = 10, day=1) {
  switch(ifelse(is.null(vaccinePars$type), "other", vaccinePars$type),
         booster = make_booster(name, vaccinePars, simLength = 10, day=1),
         make_vaccine_base(name, vaccinePars, simLength = 10, day=1))
}

make_vaccine_base = function(name, vaccinePars, simLength = 10, day=1) {
  stopifnot(length(vaccinePars$Vschedule) >= simLength)
  force(name)
  S = c(vaccinePars$S0, rep(0, simLength))
  #C = S = c(vaccinePars$S0, rep(0, simLength))
  t = 2  # t=1 is initialization; implies Vschedule[t-1]
  S[t] = S[t - 1]
  
  list(
    state = function() {
      list(
        name = name, day = t, S = S[t], 
        Vschedule = vaccinePars$Vschedule[t - 1], I = NA
      )
    },
    history = function() {
      list(S = S, Vschedule = vaccinePars$Vschedule)
    },
    pars = function() {
      vaccinePars
    },
    get_vaccGroups = function(compartments) {
      str_subset(names(compartments), "[:digit:]", negate=T)
    },
    vaccinate = function(nVacc) {
      S[t] <<- S[t] + nVacc
      if (nVacc < vaccinePars$Vschedule[t - 1]) {
        warning(paste("More vaccine", name, "than candidates."))
      }
    },
    next_day = function() {
      t <<- t + 1
      S[t] <<- S[t - 1]
    },
    remove_susceptibles = function(toRemove) {  # called when reinfected
      S[t] <<- S[t] - toRemove
    },
    add_susceptibles = function(toAdd) {  # called when vaccinate
      S[t] <<- S[t] + toAdd
    }
  )
}

make_booster = function(...) {
  vaccine = make_vaccine_base(...)
  env = environment(vaccine$state)
  vaccine$get_vaccGroups = function(compartments) {
    names(compartments) %>% 
      str_subset("[:digit:]") %>%  # have been vaccinated
      str_subset(paste0(env$name,"$"), negate=T)  # last intervention wasn't boost
    #   str_subset("[-[:lower:]]$")  # waned
  }
  
  vaccine
}

# Control Functions -------------------------------------------------------
# currently the code only functions correctly for a single control
# for multiple controls, need to decide how to "stack" changes, or take max etc
infection_control = function(pars, simLength = 10, name = "newInf_mitig") {
  if (!is.null(pars$delta)) stopifnot(pars$type %in% c("mult", "add"))
  
  delta_up = pars$delta_up
  delta_down = pars$delta_down
  mitig = cases_restrict = cases_relax =
    R_restrict = R_relax = rep(NA, simLength + 1)
  mitig[1] = pars$mitig
  t = 1
  tPossible = 0  # can intervene immediately if need be
  
  list(
    state = function() {
      list(name = name, t = t, mitig = mitig[t])
    },
    history = function() {
      list(
        mitig = mitig, cases_restrict = cases_restrict,
        cases_relax = cases_relax, R_restrict = R_restrict, R_relax = R_relax
      )
    },
    pars = function() {
      pars
    },
    control = function(stats, compartments) {
      t2 = t - pars$delay
      mit = mitig[t]
      if (!(t2 < 1 | t < tPossible)) { # able to intervene
        # compute case number averages
        cases = pars$detRatio *
          mean(stats[max(1, t2 - pars$average + 1):t2, "newInfs"])
        if ((stats[t2, "R_e_t"] > pars$upper_R_lockdown) |
            ((cases > pars$upper) & (stats[t2, "R_e_t"] > pars$lower_R))) {
          # Increase Mitigation
          cases_restrict[t] <<- cases
          R_restrict[t] <<- stats[t2, "R_e_t"]
          mit = ifelse(pars$type == "mult", 
                       delta_up + mit * (1 - delta_up),
                       mit + delta_up)
        } else if ((cases < pars$lower) & (stats[t2, "R_e_t"] < pars$upper_R)) {
          # Decrease Mitigation
          cases_relax[t] <<- cases
          R_relax[t] <<- stats[t2, "R_e_t"]
          mit = ifelse(pars$type == "mult",
                       (mit - delta_down) / (1 - delta_down),
                       mit - delta_down)
        }
        mit = min(max(mit, 0), 1)
      }
      mitig[t + 1] <<- mit  # change takes effect tomorrow, based on "old" info
      strainNames = get_strainNames(compartments)
      for (g in compartments[strainNames]) { g$update_mitigation(mit) }
      if (mit != mitig[t]) { # a change was made, cannot intervene for pars$gap
        tPossible <<- t + pars$gap
      }
      compartments
    },
    next_day = function() {
      t <<- t + 1
    }
  )
}

optimal_control = function(...) {  # only par is delay
  control = infection_control(...)
  e = environment(control$state)
  delay = e$pars$delay
  e$name = paste0("optimal_d", delay)
  
  control$control = function (stats, compartments) {
    t2 = e$t-delay
    mit = ifelse(t2<=1, e$mitig, stats[e$t-delay, "M_star"])
    e$mitig[e$t+1] = min(max(mit, 0), 1)  # change takes effect tomorrow
    strainNames = get_strainNames(compartments)
    for (g in compartments[strainNames]) { g$update_mitigation(mit) }
    
    compartments
  }
  
  control
}

no_control = function(...) {  # only par is delay
  control = infection_control(...)
  e = environment(control$state)
  e$mitig = rep(e$mitig[1], length(e$mitig))
  e$name = "no_control"
  
  control$control = function (stats, compartments) {
    if (e$t==1) {
      strainNames = get_strainNames(compartments)
      for (g in compartments[strainNames]) { g$update_mitigation(e$mitig[1]) }
    }
    
    compartments
  }
  
  control
}

fixed_control = function(mitig_vec, ...) {  # only par is delay
  control = infection_control(...)
  e = environment(control$state)
  stopifnot(length(e$mitig) != mitig_vec)
  e$mitig = mitig_vec
  e$name = "fixed_control"
  
  control$control = function (stats, compartments) {
    strainNames = get_strainNames(compartments)
    for (g in compartments[strainNames]) { g$update_mitigation(e$mitig[e$t]) }
    
    compartments
  }
  
  control
}


# Resistance Waning -------------------------------------------------------

wane_compartments = function(compartments, pars_population, simLength, day) {
  # restrict based on "non-waning" criteria: 3x intervention (first vaccine=2x)
  # updated to 4x intervention; because otherwise more booster vacc than waned
  restrict_waning = function(toWane) {
    toWane %>% 
      .[map_lgl(., ~ str_count(., "[:alpha:]") +
                  str_count(., "[:digit:]") + str_detect(., "[:digit:]") < 4)]
  }
  # first select "possible" sets
  vaccines = names(compartments) %>%  # unwaned vaccine
    .[map_lgl(., ~ str_count(., "[:digit:]") > str_count(., "-"))] %>% 
    restrict_waning()
  variants = names(compartments) %>% 
    str_subset("[:upper:][:digit:]?$") %>% 
    restrict_waning()
  compartments = wane(vaccines, compartments, pars_population, 
                      simLength, day, vaccines=T)
  compartments = wane(variants, compartments, pars_population, 
                      simLength, day, vaccines=F)
  
  compartments
}

wane = function(namesToWane, compartments, pars_population, simLength, day, 
                vaccines=T) {
  for (name in namesToWane) {
    # see if need to do anything
    nWane = compartments[[name]] %$% 
      rpois(1, state()$S/pars()$nDaysWane) %>% 
      min(compartments[[name]]$state()$S)  # can't remove more than is in group
    if (is.na(compartments[[name]]$state()$S) | nWane==0) { next }
    # get new name; create new compartment if needed
    waneName = ifelse(vaccines, paste0(name,"-"),
                      gsub(x = name, "([[:upper:]])(-?)$", "\\L\\1\\2", perl=T))
    if (!(waneName %in% sapply(compartments, function(com) com$state()$name))) {
      if (vaccines) {
        interventionName = stringi::stri_extract_last_regex(name, "[:digit:]")
      } else {
        interventionName = stringi::stri_extract_last_regex(name, "[:upper:]")
      } 
      compartments = add_compartment(compartments, name, interventionName, 
                                     waneName, pars_population, simLength, 
                                     day, type="wane")
    }
    # remove people from from previous compartment, add to this one
    compartments[[name]]$remove_susceptibles(nWane)
    compartments[[waneName]]$add_susceptibles(nWane)
  }
  
  compartments
}


# Functions for Daily Statistics ------------------------------------------
newInfs = function(compartments, strainNames, ...) {
  sum(map_dbl(compartments[strainNames], ~ .$state()$I))
}

#use Lambda for the computations instead of Theta
R_e_t = function(compartments, controls, strainNames, seasonality = 1, ...) {
  mitig = min(map_dbl(controls, ~ .$state()$mitig))
  strains = compartments[strainNames]
  Re = (1 - mitig) * seasonality *
    sum(map_dbl(strains, ~ .$state()$GammaTilde * .$state()$Lambda)) /
    sum(map_dbl(strains, ~ .$state()$W))
  ifelse(is.nan(Re), 0, Re)
}

rhoV = function(compartments, parsList, ...) {
  infNames = parsList$strains$names %>% 
    unique(str_sub(., -1)) %>% 
    setNames(., .)
  map(infNames, ~ parsList$population$N^{-1} * parsList$strain$pars[[.]]$R *
        sum(map_dbl(compartments, function(g) g$state()$S * (1-g$pars()$res[[.]]))))
}


M_star = function(compartments, strainNames, R_bar = 1, seasonality = 1, ...) {
  strains = compartments[strainNames]
  Ms = 1 - R_bar * seasonality^-1 * sum(map_dbl(strains, ~ .$state()$W)) /
    sum(map_dbl(strains, ~ .$state()$GammaTilde * .$state()$Lambda))
  ifelse(is.nan(Ms), 0, Ms)
}


# Simulation Functions ----------------------------------------------------

#  must provide all base strains first in parsList$strains$names
sim_init = function(parsList, simLength = 10) {
  compartments = list()
  # fully susceptible group
  compartments[["_"]] = make_vaccine(  
    "_",
    list(
      res = as.list(setNames(
        rep(0, length(parsList$strains$names)), parsList$strains$names
      )),
      S0 = with(parsList$population, floor(pSuscept * N)),
      Vschedule = rep(0, simLength)
    ),
    simLength
  )
  # variants
  for (i in seq_along(parsList$strains$names)) {
    compartments[[parsList$strains$names[i]]] =
      make_strain(
        parsList$strains$names[i],
        parsList$strains$pars[[i]],
        parsList$population,
        simLength
      )
  }
  # vaccines
  for (i in seq_along(parsList$vaccines$names)) {
    compartments[[parsList$vaccines$names[i]]] =
      make_vaccine(parsList$vaccines$names[i],
                   parsList$vaccines$pars[[i]],
                   simLength)
  }
  
  compartments
}

infect = function(compartments, strainName, susceptibleName, pars_population,
                  simLength, day) {
  infName = compartments[[strainName]]$state()$infName
  susceptible = compartments[[susceptibleName]]
  newInfs = compartments[[strainName]]$create_new_infs(susceptible)
  assign("env", environment(), envir = .GlobalEnv)
  if (newInfs$I == 0) { return(compartments) }  # no infections; return early
  susceptible$remove_susceptibles(newInfs$I)
  if (susceptibleName == "_") { # no interaction
    compartments[[infName]]$add_new_infs(newInfs)
  } else {
    interactionName = susceptibleName %>% 
      str_replace(infName, "") %>% 
      str_replace(tolower(infName), "") %>%
      paste0(infName)
    # if new name, create strain group; only res really matters
    if (!(interactionName %in% names(compartments))) {
      compartments = add_compartment(compartments, susceptible$state()$name,
                                     infName, interactionName, pars_population, 
                                     simLength, day, type="variant")
    }
    compartments[[interactionName]]$add_new_infs(newInfs)
  }
  compartments
}

vaccinate = function(compartments, vaccineName, pars_population, simLength, 
                     day) {
  vaccine = compartments[[vaccineName]]
  
  # who gets vaccinated
  vaccGroups = vaccine$get_vaccGroups(compartments)
  groupSizes = map_dbl(compartments[vaccGroups], ~ .$state()$S)
  if (sum(groupSizes)==0) return(compartments)
  nVacc = min(sum(groupSizes), vaccine$state()$Vschedule)
  toVacc = rep(seq_along(groupSizes), times=groupSizes) %>% 
    sample(size = nVacc) %>% 
    table()
  names(toVacc) = names(groupSizes)[as.integer(names(toVacc))]
  vaccine$vaccinate(nVacc)  # add vaccinated people to this group susceptibles
  
  # move vaccinated people to appropriate groups
  for (name in names(toVacc)) {
    if (name == "_") { # in this case they stay in vaccine class
      compartments[[name]]$remove_susceptibles(toVacc[name])
    } else { # move to interaction class of appropriate strain
      vaccine$remove_susceptibles(toVacc[name]) #  remove from vaccine group
      compartments[[name]]$remove_susceptibles(toVacc[name])  # came from here
      interactionName = paste0(name, vaccineName)
      # if new name, create strain group; only res really matters
      if (!(interactionName %in% names(compartments))) {
        compartments = add_compartment(compartments, name, vaccineName,
                                       interactionName, pars_population, 
                                       simLength, day, type="vaccine")
      }
      compartments[[interactionName]]$add_susceptibles(toVacc[name])
    }
  }
  compartments
}


sim_multistrain = function(parsList, simLength, wane=F,
                           pars_controls=list(type="none"), history=F) {
  # write a check for options if wane=T: omega, nDaysWane
  if (is.null(parsList$population$seasonality)) {  # initialize seasonality
    parsList$population$seasonality = seasonality(parsList$population)
  }
  seasonality = parsList$population$seasonality
  iTime = parsList$population$iTime
  names(parsList$strains$pars) = parsList$strains$names 
  statNames = c("newInfs", "R_e_t", "M_star", "rhoV") # needs to match function name
  statFncs = map(statNames, get)
  names(statFncs) = statNames
  
  
  compartments = sim_init(parsList, simLength)
  
  if (pars_controls$type == "optimal") {
    pars_controls$type = "mult"
    controls = list(optimal_control(pars_controls, simLength))
  } else if (pars_controls$type == "none") {
    controls = list(no_control(pars_controls, simLength))
  } else {
    controls = list(infection_control(pars_controls, simLength))
  }
  
  names(controls) = map_chr(controls, ~ .$state()$name)
  for (co in controls) { compartments = co$control(stats, compartments) }
  N = parsList$population$N
  
  strainNames = get_strainNames(compartments)
  # Initial Statistics
  for (g in compartments[strainNames]) {
    g$compute_summary(compartments, yesterday=T)
  }
  stats_init = parsList$population %$% 
    map(statFncs, exec, compartments, controls = controls, R_bar = R_bar,
        strainNames = strainNames, seasonality = seasonality[iTime],
        parsList=parsList) %>% 
    unlist()
  stats = matrix(NA_real_, nrow = simLength + 1, ncol = length(stats_init))
  colnames(stats) = names(stats_init)
  stats[1, ] = stats_init
  
  for (day in seq_len(simLength)) {
    # each strain causes new infections
    for (strainName in strainNames) {
      for (susceptibleName in names(compartments)) {
        compartments = compartments %>% 
          infect(strainName, susceptibleName, parsList$population, simLength, day)
      }
    }
    # each vaccine "vaccinates" more people, proportional to group size
    for (vaccineName in parsList$vaccines$names) {
      compartments = compartments %>% 
        vaccinate(vaccineName, parsList$population,simLength, day)
    }
    
    # daily statistics
    strainNames = get_strainNames(compartments)
    for (g in compartments[strainNames]) { g$compute_summary(compartments) }
    stats[day + 1, ] = parsList$population %$% 
      map(statFncs, exec, compartments, controls = controls, R_bar = R_bar,
          strainNames = strainNames, seasonality = seasonality[day+iTime],
          parsList=parsList) %>% 
      unlist()
    
    #delete day=0 line for use of stats in control
    for (co in controls) compartments = co$control(stats[-1, ], compartments)
    if (wane) {
      compartments = compartments %>% 
        wane_compartments(pars_population, simLength, day)
    }
    for (g in compartments) g$next_day()
    for (co in controls) co$next_day()  # here for only a single if() check
  }
  
  if (history) {
    out = list(compartments = map(compartments, ~.$history()), stats = stats)
    if (!is.null(pars_controls)) out$controls = map(controls, ~.$history())
  } else {
    out = list(compartments = compartments, stats = stats)
    if (!is.null(pars_controls)) out$controls = controls
  }
  
  out
}


runSim = function(nRuns = 20, pars_compartments, pars_controls, simLength=10) {
  nuCores = detectCores() - 1
  registerDoParallel(nuCores)  # set up cluster and register
  output = foreach(i = seq_len(nRuns)) %dopar% {
    sim_multistrain(
      parsList = pars_compartments,
      pars_controls = pars_controls,
      simLength = simLength
    )
  }
  stopImplicitCluster()
  output
}
