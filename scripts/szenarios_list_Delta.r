latest_EpiNow_R_estimates = read.csv("data/EpiNow_R_estimates.csv") %>% 
  mutate(date = as.Date(date))

#########################################################
# SCENARIO LIST
#########################################################
szenarios_list = list(
  #########################################################
  ### 20% Delta, 2021-06-12
  list(
    forecasting_date = as.Date("2021-06-12"),
    delta_prev = 0.20,
    geometric_growth_model = T,
    vaccination_plan_with_past_data = T,
    R_observed = latest_EpiNow_R_estimates %>%
      filter(date <= "2021-06-12") %>% pull(median) %>% tail(1),
    name = "2021-06-12",  # -medianR
    import = list(
      W = 0,
      A = 0,
      B = 1,
      G = 1,
      D = 1
    )
  )
)

controls_list = list(
  # 25-50
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1.0,
    upper_R_lockdown = 1.2,
    lower_R = 0,
    delay = 7,
    gap = 14,
    lower = 318 / 2,
    upper = 318,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "proactive_25-50_delay7"  # -Rlockdown12
  ),
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1.0,
    upper_R_lockdown = 10,
    lower_R = 0,
    delay = 7,
    gap = 14,
    lower = 318 / 2,
    upper = 318,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "reactive_25-50_delay7"
  ),
  # 25-150
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1.0,
    upper_R_lockdown = 1.2,
    lower_R = 0,
    delay = 7,
    gap = 14,
    lower = 318 / 2,
    upper = 318 * 3,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "proactive_25-150_delay7"  # -Rlockdown12
  ),
  list(
    delta_up = 0.2,
    delta_down = 0.2,
    upper_R = 1.0,
    upper_R_lockdown = 10,
    lower_R = 0,
    delay = 7,
    gap = 14,
    lower = 318 / 2,
    upper = 318 * 3,
    type = "mult",
    average = 14,
    detRatio = 1/1.4,
    name = "reactive_25-150_delay7"
  )
)

