# General Population Parameters -------------------------------------------

#N = 8902600 # total population
N = 8932664 # total population
R0 = 3.5

# Import Infection Data ---------------------------------------------------

### Import cases history from WHO
cases = read.csv("data/cases_who.csv") %>% 
  mutate(date = as.Date(date)) %>% 
  rename(cases_total = cases) %>% 
  filter(date <= forecasting_date)

### Account for unobserved infections
infectionrate_2020 = 0.12 # wild type infections in 2020

# dark figures for different time periods
dark_dates = as.Date(c("2020-01-01", "2021-01-01", "2021-03-01", "2021-04-01"))
dark_figures = c(NA, 2.3, 2, 1.4)

infections_2020 = round(N * infectionrate_2020)
cummulative_cases_2020 = cases %>%
  filter(date >= "2020-01-01", date <= "2020-12-31") %>%
  pull(cases_total) %>%
  sum
dark_figures[1] = infections_2020 / cummulative_cases_2020

cases = cases %>% 
  mutate(darkFig = dark_figures[map_int(date, ~sum(.>=dark_dates))],
         cases_corrected = round(cases_total * darkFig),
         # previously_infected_total = cumsum(replace_na(cases_corrected, 0)),
         cases_ra = round(stats::filter(cases_corrected, rep(1 / 7, 7), sides = 1))) 

cases_init = tail(cases$cases_ra, 13) #/ifelse(detRatio==1, 1.4, 1)

### Compute weekly cases (from corrected cases)
cases_weekly = cases %>% 
  mutate(year = lubridate::year(date),
         KW = lubridate::epiweek(date)) %>% 
  group_by(year, KW) %>% 
  summarise(cases_weekly = sum(cases_corrected))

# Import Variant Data -----------------------------------------------------

### Import variant data taken from the AGES website
prevalence = read.csv("data/ages_variant_data.csv") %>%
  filter(KW <= lubridate::epiweek(forecasting_date))
# delta_prev = tail(prevalence$D, 1)
weekly_variants =  prevalence %>% 
  fill(everything()) %>% 
  left_join(filter(cases_weekly, year==2021), by="KW") %>% 
  mutate_at(vars(W:D), ~round(.*cases_weekly))

### Compute initial compartment sizes
S0 = weekly_variants %>% 
  select(W:D) %>% 
  colSums

missing_strains = setdiff(strain_names, names(S0))
S0[missing_strains] = 0
S0 = as.list(S0)
S0$W = S0$W + infections_2020 #add wild type infections of 2020 to this list
S0_unvaccinated =  S0  # used in paper_plots.r

# end major updates
# Import Vaccination Data -------------------------------------------------

### PARAMETERS
vaccines = c("BioNTechPfizer", "Moderna", "AstraZeneca", "Janssen")
immunity_delay = c(2 * 7, 2 * 7, 4 * 7, 4 * 7)
immunity_delay_booster = 7

#maximum_vaccine_uptake = 0.85 # percentage of whole population which will get a vaccine
vaccination_data_aut = read.csv("data/vaccination_data_aut.csv") %>% 
    mutate(date = as.Date(date)) %>%
    select(date, all_of(vaccines))

vaccination_data_aut_trunc = filter(vaccination_data_aut, date <= forecasting_date)
#get latest data for group size
vaccination_group_sizes = 
  foreach(vaccine = vaccines, .combine = cbind) %do% {
    data.frame(vax = tail(vaccination_data_aut_trunc[[vaccine]],
                          immunity_delay[vaccine == vaccines])[1])
  } %>% set_names(vaccines)

total_vaccines = sum(unlist(vaccination_group_sizes))

vaccination_group_sizes %<>%
  mutate(
    mRNA   = BioNTechPfizer + Moderna,
    vector = AstraZeneca + Janssen,
    booster = 0
  )

# Future Vaccinations/Vaccination Plan
# last couple of first doses and paste this together with daily doses

plan_length = max(simLength, 300)

if (!vaccination_plan_with_past_data) {
  daily_doses = foreach(i = 1:length(vaccines), .combine = c) %do% {
    vaccination_data_aut_trunc[[vaccines[i]]] %>% diff %>% tail(7) %>% median %>% round
  }
  
  vaccination_plan = foreach(i = 1:length(vaccines), .combine = cbind) %do% {
    first_dose_past = vaccination_data_aut_trunc[[vaccines[i]]] %>% 
      diff %>% tail(immunity_delay[i])
    c(first_dose_past,
      rep(daily_doses[i], plan_length - immunity_delay[i]))
  } %>% as.data.frame %>% set_names(vaccines)
} else { #### TAKE DATA FOR THE VACCINATION PLAN
  daily_doses = foreach(i = 1:length(vaccines), .combine = c) %do% {
    vaccination_data_aut %>% 
      filter(date <= vaccination_data_cutoff) %>% 
      .[[vaccines[i]]] %>% 
      diff %>% tail(7) %>% median %>% round
  }
  
  vaccination_plan = foreach(i = 1:length(vaccines), .combine = cbind) %do% {
    first_dose_past = vaccination_data_aut_trunc[[vaccines[i]]] %>% 
      diff %>% tail(immunity_delay[i])
    first_dose_data = vaccination_data_aut %>% 
      filter(date > as.Date(forecasting_date) & 
               (date <= as.Date(vaccination_data_cutoff))) %>% 
      .[[vaccines[i]]] %>% diff
    c(first_dose_past, first_dose_data,
      rep(daily_doses[i],
          max(0, plan_length - length(first_dose_past) - length(first_dose_data)))
    ) %>% tail(plan_length)
  } %>% as.data.frame %>% set_names(vaccines)
  
}

if (!toBoost) {
  Vschedule_booster = 0
} else {
  booster_data_aut = read.csv("data/booster_data_aut.csv") %>% 
    mutate(date = as.Date(date),
           doses = c(NA, diff(cumulativeDoses))) %>% 
    filter(date >= forecasting_date - immunity_delay_booster)
  
  daily_doses = booster_data_aut %>% 
    # filter(date <= vaccination_data_cutoff) %>% 
    pull(cumulativeDoses) %>% 
    diff %>% tail(7) %>% median %>% round
  
  booster_data = diff(pull(booster_data_aut, cumulativeDoses))
  Vschedule_booster = 
    c(booster_data, 
      rep(daily_doses, max(0, plan_length - length(booster_data)))
  ) %>% tail(plan_length)
}


# stop vaccinations when maximum vaccination uptake is reached
vaccination_plan[cumsum(rowSums(vaccination_plan)) +
                   total_vaccines  >  N * maximum_vaccine_uptake, ] = 0

vaccination_plan %<>%
  mutate(
    mRNA   = BioNTechPfizer + Moderna,
    vector = AstraZeneca + Janssen,
    booster = Vschedule_booster
  )


# Compute Summary Parameters ----------------------------------------------

### Overlap between infections and vaccinations
total_infected = S0 %>% unlist %>% sum
total_vaccinated = total_vaccines
total_uninfected = N - total_infected

### distribute vaccines between infected and uninfected population proportionally
uninfected_vaccinated = round(total_vaccinated * total_uninfected / N)
infected_vaccinated = round(total_vaccinated * total_infected / N)

#Compute proportion of population that was neither infected nor vaccinated
pSus = (total_uninfected - uninfected_vaccinated) / N #needed for setting up the simulation later

#Distribute vaccines among previously infected population
#ASSUMPTION: ONLY VACCINATE WILD AND ALPHA STRAIN!
W_vaccinated = round(S0$W * infected_vaccinated / (S0$W + S0$A))
A_vaccinated = round(S0$A * infected_vaccinated / (S0$W + S0$A))

#correct the strain groups
S0$W = S0$W - W_vaccinated
S0$A = S0$A - A_vaccinated

