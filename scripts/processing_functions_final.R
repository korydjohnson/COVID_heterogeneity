# processing --------------------------------------------------------------

get_incidence = function(compartments_, simLength, nu=13) {
  enframe(compartments_, value="history") %>%
    mutate(variant = str_sub(name, -1)) %>%
    filter(str_detect(variant, "[:Upper:]")) %>%
    mutate(I = map(history, ~tail(.$I, simLength+nu))) %>%
    group_by(variant) %>%
    summarise(I = list(reduce(I, `+`))) %>%
    add_row(variant = "Total", I = list(reduce(.$I, `+`))) %>%
    mutate(prevalence = map(I, function(i) i/.$I[[nrow(.)]]),
           prevalence = map(prevalence, ~ifelse(is.nan(.), 0, .))) %>%
    ungroup()
}

get_quantiles = function(grTib, q) {
  colExpand = setdiff(colnames(grTib), group_vars(grTib))
  grTib %>% 
    summarise_all(function(col) list(quantile(col, probs = q, na.rm=T))) %>% 
    unnest(all_of(colExpand)) %>% 
    ungroup() %>% 
    mutate(q_name = rep(as.character(1:length(q)), nrow(.)/length(q))) %>% 
    pivot_wider(names_from = q_name, 
                values_from = all_of(colExpand)) %>% 
    set_names(c(group_vars(grTib), unlist(map(colExpand, ~paste0(., "_", q * 100, "%")))))
}

process_simData = function(simData, q = c(0.025, 0.25, 0.5, 0.75, 0.975), nu=13) {
  simLength = nrow(simData[[1]]$stats)-1
  
  incidenceTib = enframe(simData, name="sample") %>% 
    mutate(incidenceTib = map(value, ~get_incidence(.$compartments, simLength))) %>% 
    select(-sample, -value) %>% 
    unnest(incidenceTib) %>% 
    mutate(day = list(-(nu-1):simLength)) %>%  # plot data from initialization
    unnest(c(I, prevalence, day)) %>% 
    group_by(day, variant) %>%
    get_quantiles(q)
  
  statTib = enframe(simData, name="sample") %>%
    unnest_wider(value) %>% 
    mutate(mit = map(controls, ~ 1 - .[[1]]$mitig),  # we plot \tilde M
           stats = map(stats, ~cbind(., day=1:nrow(.)-1)),
           stats = map(stats, as_tibble)) %>%  
    select(-compartments, -controls, -sample) %>%
    unnest_wider(stats) %>% 
    select(day, everything()) %>% 
    unnest(everything()) %>% 
    group_by(day) %>% 
    get_quantiles(q)
  
  list(
    stats = statTib,
    incidence = incidenceTib
  )
}

get_comparisons = function(simOut, scale_by = 10^5/8932664, detRatio=1/1.4) {
  mtilde = 1-simOut$controls[[1]]$mitig
  Rhat = simOut$stats[, "R_e_t"]
  mstar = mtilde/Rhat
  
  M_mean = mean(mtilde)
  M_median = median(mtilde)
  M_min = min(mtilde)
  M_max = max(mtilde)
  
  population = sum(map_dbl(simOut$compartments, ~tail(.$S, 1)))
  infections = sum(simOut$stats[,"newInfs"])* scale_by
  quarantine = max(
    stats::filter(simOut$stats[,"newInfs"], rep(1, 10), sides = 1)*
      scale_by * detRatio, 
    na.rm=T
  )
  
  rm(simOut)
  as.list(environment())
}

compare_simData = function(simData_both) {
  simData_both %>% 
    map(get_comparisons) %>% 
    enframe(name="sample") %>%
    separate(sample, sep="_", 
             into = c(NA, NA, "control", "bounds", "delay", "sample")) %>% 
    unnest_wider(value)
}

write_simData = function(file, csv_folder = "") {
  data_name = str_sub(file, start = 9, end = -7)
  load(paste0(data_folder, file))
  data = process_simData(simData)
  write.csv(
    data$mitigation_data,
    file = paste0(csv_folder, "mitData-", data_name, ".csv"),
    row.names = F
  )
  write.csv(
    data$multi_variant_forecast,
    file = paste0(csv_folder, "forecast-", data_name, ".csv"),
    row.names = F
  )
}

