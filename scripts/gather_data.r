resistance_data = sample_variant_data_truncnormal(resistance_ranges, 
                                                  variant_names = strain_names)

wane_data = wane_ranges %>% 
  sample_variant_data_truncnormal(wane_categories) %>% 
  cbind(setNames(rep(list(.$V), length(strain_names)-1),
                 str_subset(strain_names, "O", negate=T))) %>% 
  select(-V)

vaccine_pars = list(
  names = vaccine_names,
  pars = foreach(vacc = vaccine_names_full) %do% {
    list(
      res = resistance_data %>% filter(name == vacc) %>% select(-name, -type) %>% as.list,
      ###
      S0   = vaccination_group_sizes[[vacc]],
      Vschedule = vaccination_plan[[vacc]],
      omega = wane_data %>% filter(name == vacc) %>% select(-name, -type) %>% as.list,
      nDaysWane = 180
    )
  }
)

if (toBoost) {
  vaccine_pars[[2]][[which(vaccine_names_full == "booster")]]$type = "booster"
}

### LAMBDA ###
lambda = filter(resistance_data, name == "lambda") %>% 
  select_if(is.double)

### LAMBDA ###
E_gT = filter(resistance_data, name == "E[gT]") %>% 
  select_if(is.double)


### STRAIN RESISTANCES ###
strain_prevalence =
  list(A = 1 - delta_prev,
       D = delta_prev)

missing_strains = setdiff(strain_names, names(strain_prevalence))
strain_prevalence = c(strain_prevalence,
                      setNames(as.list(rep(0, length(missing_strains))), 
                               missing_strains))[strain_names]

#strain_prevalence$W = 1 - sum(unlist(strain_prevalence))
### give initial infectious population as (constant) proportion of observed cases
### this is also wrong; these aren't the values subtracted from population size
I0 = lapply(strain_prevalence, function(a) {round(a * cases_init)})

strain_pars = foreach(strain = strain_names) %do% {
  list(
    name = strain,
    iTime = 13,
    w = generationInt(iTime = 13, E = 4.46, Var = 2.63 ^ 2),  # E_gT[[strain]]
    k = 0.1,
    R = R0 * lambda[[strain]],
    res = resistance_data %>% filter(name == strain) %>% select(-name, -type) %>% as.list,
    ###
    I0 = I0[[strain]],
    S0 = S0[[strain]],
    I_import = import[[strain]],
    omega = wane_data %>% filter(name == strain) %>% select(-name, -type) %>% as.list,
    nDaysWane = 180
  )
}

### COMBINE PARAMTERS TO ONE LIST ###
pars_compartments = list(
  population = list(
    pSuscept = pSus,
    N = N,
    iTime = 13,
    w = generationInt(iTime = 13, E = 4.46, Var = 2.63 ^ 2),  # E_gT[[strain]]
    I0 = cases_init,
    date = forecasting_date,
    R_bar = R_observed,
    e_seasonality = 0.4
  ),
  strains = list(names = strain_names,
                 pars = strain_pars),
  vaccines = vaccine_pars
)

if (geometric_growth_model) {
  L = seasonality_day(date = forecasting_date)
  
  N = pars_compartments$population$N
  #strain_names = pars_compartments$strains$names
  
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
  #print(groups)
  W = with(pars_compartments$population, sum(rev(w) * I0))
  group_sizes = sapply(groups, function(g) {g$S0})
  
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
  
  mit_tilde_approx = R_observed * W / (L * sum(Gamma_g * Lambda_g)) #preliminary mitigation
  
  ###compute R_e_t^g
  R_e_t_g = mit_tilde_approx * L * unlist(lambda) * Gamma_g * R0
  #R_e_t_g
  
  growth_rate_total = estimate_growth_rate(pars_compartments$population, R = R_observed)
  
  prevalence_model =
    foreach(i = 1:length(strain_names)) %do% {
      group = pars_compartments$strains$pars[[i]]
      growth_rate = estimate_growth_rate(group, R = R_e_t_g[i])
      prev = strain_prevalence[[i]]
      iTime = pars_compartments$strains$pars[[i]]$iTime
      rev(prev * (growth_rate_total / growth_rate) ^ (0:(iTime - 1)))
    } %>% set_names(strain_names)
  
  prevalence_model$A = 1 - prevalence_model$D
  # adjust initial incidence
  for (i in 1:length(strain_names)) {
    pars_compartments$strains$pars[[i]]$I0 =
      round(cases_init * prevalence_model[[i]])
  }
  
}

#readjust current mitigation
pars_controls$mitig = max(
  compute_initial_mit(
    pars_compartments = pars_compartments,
    R_observed = R_observed,
    L = seasonality_day(date = forecasting_date,
                    e = pars_compartments$population$e_seasonality) #initial seasonality
  ),
  0
)
