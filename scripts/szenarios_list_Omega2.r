latest_EpiNow_R_estimates = read.csv("data/EpiNow_R_estimates.csv") %>% 
  mutate(date = as.Date(date))

#########################################################
# SCENARIO LIST
#########################################################
szenarios_list = list(
  list(
    forecasting_date = as.Date("2021-08-08"),
    delta_prev = 1,
    geometric_growth_model = F,
    vaccination_plan_with_past_data = F,
    R_observed = latest_EpiNow_R_estimates %>%
      filter(date <= "2021-08-08") %>% pull(median) %>% tail(1),
    name = "2021-08-08",  # -medianR-import1
    import = list(
      W = 0,
      A = 0,
      B = 0,
      G = 0,
      D = 0,
      O = 1
    )
  )
)

#########################################################
# CONTROL LIST
#########################################################

controls_list = list(
  #########################################################
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1,
    upper_R_lockdown = 1.2,
    lower_R = 0,
    delay = 21,
    gap = 14,
    lower = 318 / 2,
    upper = 318 * 3,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "proactive_25-150_delay21"  # -Rlockdown12
  ),
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1.0,
    upper_R_lockdown = 10,
    lower_R = 0,
    delay = 21,
    gap = 14,
    lower = 318 / 2,
    upper = 318 * 3,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "reactive_25-150_delay21"
  ),
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1,
    upper_R_lockdown = 1.2,
    lower_R = 0,
    delay = 21,
    gap = 14,
    lower = 318 / 2,
    upper = 318,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "proactive_25-50_delay21"  # -Rlockdown12
  ),
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1.0,
    upper_R_lockdown = 10,
    lower_R = 0,
    delay = 21,
    gap = 14,
    lower = 318 / 2,
    upper = 318,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "reactive_25-50_delay21"
  )
)
