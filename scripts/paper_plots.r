library(tidyverse)
library(magrittr)
library(foreach)
library(doParallel)
options(dplyr.summarise.inform = FALSE)
library(rvest)
library(foreach)
library(truncnorm)

source("scripts/simulation_final.R")
nu = 13
plot_strain_plots = T
R0=3.5

toBoost=F
toWane=F
maximum_vaccine_uptake = 0.85 # percentage of whole population which will get a vaccine

######################################################
# PREVALENCE PLOT
######################################################
ages_variant_data = read.csv("data/ages_variant_data.csv")
KW_dates = as.Date("2021-01-10") + (0:51) * 7

colour_palette <-
  c("#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7")

gather(ages_variant_data %>%
         set_names(c(
           "KW", "WT", "Alpha", "Beta", "Gamma", "Delta"
         )) %>%
         filter(KW <= 34),
       Variant,
       value = value,
       -KW) %>%
  mutate(Variant = factor(Variant, levels = unique(Variant))) %>%
  mutate(date = KW_dates[KW] - 3) %>%
  ggplot(aes(x = date)) +
  geom_col(aes(y = value, fill = Variant), alpha = 0.9) +
  theme_bw() +
  theme(legend.position = "bottom",
        aspect.ratio = 0.6) +
  
  scale_fill_manual(values = colour_palette) +
  ylab("Variant Prevalence") +
  xlab("Date") +
  scale_x_date(date_breaks = "months", date_labels = "%b %y")
  
ggsave("plots/AGES-prevalence.png",
       height = 7 * 0.7,
       width = 7,  dpi=600)

######################################################
# POPULATION MAKEUP PLOT
######################################################
variant_names = strain_names = c("W", "A", "B", "G", "D", "O")
vaccination_plan_with_past_data = F

colour_palette <- c(
  "#999999",
  "#333333",
  "#666666",
  "#E69F00",
  "#56B4E9",
  "#009E73",
  "#F0E442",
  "#0072B2",
  NA,
  "#CC79A7",
  "#330066"
  # "#D55E00", "#CC79A7"
)
#do JUNE population
forecasting_date = as.Date("2021-06-12")
sim_end = as.Date("2022-05-01")
simLength_with_end = T
if (simLength_with_end) {
  simLength = as.integer(sim_end - forecasting_date)
}
source("scripts/multi_variant_population_parameters_imports.r")
population_makeup = c(total_uninfected - uninfected_vaccinated, unlist(S0)) /
  (N - total_vaccinated)

names(population_makeup) = c("_", variant_names)
names(S0) = c("WT", "Alpha", "Beta", "Gamma", "Delta", "")

population_composition_June = tibble(
  what = "Total Population Composition",
  p = c(
    pSus,
    sum(unlist(S0)) / N,
    vaccination_group_sizes$mRNA / N + vaccination_group_sizes$vector /
      N
  ),
  status = c(
    "No Previous Infection or Vaccination",
    'Previously Infected, not Vaccinated',
    'Vaccinated'
  )
) %>% rbind(tibble(
  what = "Population of Previously Infected",
  p = unlist(S0) / sum(unlist(S0)),
  status = names(S0)
)) %>% rbind(tibble(
  what = "Populations of Vaccinated",
  p = c(
    vaccination_group_sizes$mRNA / (
      vaccination_group_sizes$mRNA + vaccination_group_sizes$vector
    ),
    vaccination_group_sizes$vector / (
      vaccination_group_sizes$mRNA + vaccination_group_sizes$vector
    )
  ),
  status = c("mRNA", "vector")
)) %>%
  mutate(
    what = factor(what, levels = rev(unique(what))),
    status = factor(status, levels = unique(status)),
    #data_date = forecasting_date
    data_date = "June 12, 2021"
  )

#do AUGUST population
forecasting_date = as.Date("2021-08-08")
simLength_with_end = T
if (simLength_with_end) {
  simLength = as.integer(sim_end - forecasting_date)
}
source("scripts/multi_variant_population_parameters_imports.r")
population_makeup = c(total_uninfected - uninfected_vaccinated, unlist(S0)) /
  (N - total_vaccinated)

names(population_makeup) = c("_", variant_names)
names(S0) = c("WT", "Alpha", "Beta", "Gamma", "Delta", "")


population_composition_August = tibble(
  what = "Total Population Composition",
  p = c(
    pSus,
    sum(unlist(S0)) / N,
    vaccination_group_sizes$mRNA / N + vaccination_group_sizes$vector /
      N
  ),
  status = c(
    "No Previous Infection or Vaccination",
    'Previously Infected, not Vaccinated',
    'Vaccinated'
  )
) %>% rbind(tibble(
  what = "Population of Previously Infected",
  p = unlist(S0) / sum(unlist(S0)),
  status = names(S0)
)) %>% rbind(tibble(
  what = "Populations of Vaccinated",
  p = c(
    vaccination_group_sizes$mRNA / (
      vaccination_group_sizes$mRNA + vaccination_group_sizes$vector
    ),
    vaccination_group_sizes$vector / (
      vaccination_group_sizes$mRNA + vaccination_group_sizes$vector
    )
  ),
  status = c("mRNA", "vector")
)) %>%
  mutate(
    what = factor(what, levels = rev(unique(what))),
    status = factor(status, levels = unique(status)),
    #data_date = forecasting_date
    data_date = "August 8, 2021"
  )



rbind(population_composition_June, population_composition_August) %>%
  ggplot(aes(x = what)) +
  geom_col(aes(y = p, fill = status), position = position_fill(reverse = TRUE)) +
  ylim(c(0, 1)) +
  coord_flip() +
  theme_bw() +
  theme(
    legend.position = "bottom",
    aspect.ratio = 0.3,
    axis.title.y = element_blank(),
  ) +
  scale_fill_manual(values = c(colour_palette)) +
  labs(fill = "") +
  ylab("Proportion") +
  facet_wrap(vars(data_date), nrow = 2)
#facet_wrap(vars(data_date), nrow = 1)


#ggsave("population-composition-horizontal.png", height = 7*0.4, width = 7)
ggsave("plots/population-composition-vertical.png",
       height = 7 * 0.7,
       width = 7,  dpi=600)


if (plot_strain_plots == T) {
  ######################################################
  # STRAIN PLOT
  ######################################################
  forecasting_date = as.Date("2021-08-08")
  sim_end = as.Date("2022-05-01")
  simLength = as.integer(sim_end - forecasting_date)
  nu = 13
  R0 = 3.5
  
  ### import current population makeup of Austria
  new_samples = T
  n_samples = 500
  file_suffix = ""
  
  
  variant_names = strain_names = c("W", "A", "B", "G", "D", "O")
  variant_names_full = c("WT", "Alpha", "Beta", "Gamma", "Delta", "Omega")
  
  variant_names_selection = c("A", "B", "G", "D", "O")
  variant_names_selection_full = c("Alpha", "Beta", "Gamma", "Delta", "Omega")
  
  vaccine_names = c("mRNA", "vector")
  
  vaccination_plan_with_past_data = F
  source("scripts/multi_variant_population_parameters_imports.r")
  
  population = c(list("_" = total_uninfected),
                 S0_unvaccinated)
  
  population = unlist(population[c("_", variant_names)])
  
  colour_palette = c("#E69F00",
                     "#56B4E9",
                     "#009E73",
                     "#F0E442",
                     "#0072B2",
                     "#D55E00",
                     "#CC79A7")
  colour_names = c(variant_names_full, "total")
  colours = set_names(as.list(colour_palette), colour_names)
  colour_palette = unlist(colours[variant_names_selection_full])
  
  resistance_data = read.csv("data/resistance_parameters.csv", stringsAsFactors = F) %>%
    select(c("name", "type", "bound", all_of(variant_names))) %>%
    filter((type != "infection") | (name %in% variant_names))
  
  if (new_samples == T) {
    tic = Sys.time()
    variant_data_samples =
      foreach(i = 1:n_samples, .combine = rbind) %do% {
        resistance_data_sample = sample_variant_data_truncnormal(resistance_data, variant_names) %>% mutate(sample = i)
      }
    print(Sys.time() - tic)
  }
  
  ### PLOT STRAIN DATA
  variant_samples = foreach(vax = c("mRNA", "vector"), .combine = rbind) %do% {
    gather(
      variant_data_samples %>% filter(sample <= 500),
      variant,
      value = value,
      -name,
      -type,
      -sample
    ) %>%
      filter(name %in% c("lambda", vax)) %>%
      select(-type) %>%
      spread(name, value) %>%
      mutate(R0 = lambda * R0) %>%
      filter(variant %in% variant_names_selection) %>%
      mutate(vax = vax) %>% 
      rename(vaccine = !!sym(vax))
  }
  
  for (i in 1:length(variant_names)) {
    variant_samples$variant[variant_samples$variant == variant_names[i]] = variant_names_full[i]
  }
  
  variant_samples  %>%
    mutate(variant = factor(variant)) %>% #, levels = variant_names_full)) %>%
    ggplot(aes(x = vaccine, y = R0)) +
    geom_point(aes(colour = variant, shape = factor(vax)), alpha = 0.3) +
    theme_bw() +
    theme(legend.position = "bottom",
          aspect.ratio = 0.6) +
    scale_colour_manual(values = colour_palette) +
    xlab("Vaccine Effectiveness") +
    ylab(expression(paste(lambda ^ V, R[0]))) +
    labs(colour = "Variant", shape = "Vaccine") +
    xlim(c(0, 1))
  
  ggsave(
    paste0("plots/variant-samples-truncnormal", file_suffix, ".png"),
    width = 7,
    height = 7 * 0.7, dpi=600
  )
  
  ######################################################
  # RHO PLOT
  ######################################################
  group_names = c("_",
                  variant_names,
                  foreach(vaccine = vaccine_names, .combine = c) %do% {
                    c(vaccine, paste0(variant_names, vaccine))
                  })
  r = seq(0, 1, 0.01)
  ratio_mRNA = vaccination_group_sizes$mRNA / 
    (vaccination_group_sizes$mRNA + vaccination_group_sizes$vector)
  ratio_vector = 1 - ratio_mRNA
  data =
    cbind(
      (1 - r) %>% map_dfr( ~ . * population),
      r %>% map_dfr( ~ . * population * ratio_mRNA),
      r %>% map_dfr( ~ . * population * ratio_vector)
    ) %>% set_names(group_names)
  
  
  variant_data =
    foreach(i = 1:max(variant_data_samples$sample), .combine = rbind) %do% {
      resistance_data_sample = variant_data_samples %>% filter(sample == i)
      resistances = get_resistances(
        resistance_data = resistance_data_sample,
        group_names = c(variant_names, "mRNA", "vector")
      )
      #resistances of interaction with vaccine
      for (vaccine in vaccine_names) {
        for (variant in variant_names) {
          res = set_names(list(map2(
            resistances[[variant]], resistances[[vaccine]],
            ~ max(.x, .y)
          )),
          paste0(variant, vaccine))
          resistances = c(resistances, res)
        }
      }
      resistances = resistances[group_names]
      
      foreach(variant = variant_names_selection, .combine = rbind) %do% {
        lambda = resistance_data_sample %>% filter(name == "lambda") %>% pull(variant)
        
        Gamma_tilde = data %>%
          select(names(resistances)) %>%
          apply(1, function(row) {
            sum(row * (1 - sapply(resistances, function(g) {
              g[[variant]]
            })))
          })
        
        rho_t = Gamma_tilde * lambda * R0 / N
        
        tibble(
          r = r,
          Gamma_tilde = Gamma_tilde,
          rho_t = rho_t,
          variant = variant
        )
      } %>% mutate(sample = i)
    }
  
  for (i in 1:length(variant_names)) {
    variant_data$variant[variant_data$variant == variant_names[i]] = variant_names_full[i]
  }
  ### PLOT QUANTILES
  variant_data %>%
    group_by(r, variant) %>%
    summarise(quant = list(quantile(
      rho_t, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
    ))) %>%
    unnest_wider(quant) %>%
    mutate(variant = factor(variant, levels =  variant_names_full)) %>%
    ggplot(aes(x = r)) +
    # geom_vline(aes(xintercept = 0.5), linetype = "dashed", alpha = 0.5) +
    geom_line(aes(y = `50%`, colour = variant))  +
    geom_ribbon(aes(
      ymin = `25%`,
      ymax = `75%`,
      fill = variant
    ), alpha = .4) +
    geom_ribbon(aes(
      ymin = `2.5%`,
      ymax = `97.5%`,
      fill = variant
    ), alpha = .25) +
    theme_bw() +
    theme(legend.position = "bottom",
          aspect.ratio = 0.5) +
    scale_colour_manual(values = colour_palette) +
    scale_fill_manual(values = colour_palette) +
    xlab("Fully Vaccinated") +
    ylab(expression(paste(rho[t] ^ V))) +
    labs(colour = "Variant", fill = "Variant") +
    xlim(c(0, 1))
  
  ggsave(
    paste0(
      "plots/variant-rho-quantiles-aut-truncnormal",
      file_suffix,
      ".png"
    ),
    width = 7,
    height = 0.5 * 7,  dpi=600
  )
}


# Resistance Summary ------------------------------------------------------

library(reshape2)
parameters = read_csv("data/resistance_parameters.csv") %>% 
  filter(bound=="middle") %>% 
  select(name, W, A, B, G, D, O) %>%
  slice(-c(1:2)) 
ordering = c("mRNA", "vector", "booster", "W", "A", "B", "G", "D", "O")
names(ordering) = c("mRNA", "vector", "booster", "WT", "Alpha", "Beta", 
                    "Gamma", "Delta", "Omega")
molten = melt(parameters) %>% 
  mutate(name = factor(name, levels = ordering, ordered=T) %>% 
           fct_recode(!!!ordering) %>% fct_rev(),
         variable = factor(variable, levels = ordering[-c(1,2)], ordered=T) %>% 
           fct_recode(!!!ordering))
ggplot(molten, aes(x=variable, y=name, fill=value)) + 
  geom_tile() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  scale_x_discrete(position = 'top') +
  scale_fill_gradient(limits = c(0,1), name="Resistance")
ggsave("plots/resistance_parameters.png", width = 7, height = 3, 
       dpi = 600)


# average resistance due to waning ----------------------------------------

n = 1e5
nDays = 200
daysWane = 180

getWaned = function(n, nDays, daysWane) {
  full = rep(n, nDays+1)
  for (day in 2:200) {
    x = rpois(1, full[day]/daysWane)
    full[day] = full[day] - x
    full[day+1] = full[day]
  }
  full
}

wanedList = replicate(1, getWaned(n, nDays, daysWane)) %>% 
  as_tibble()

omega1 = .84
omega2 = .5
omega3 = .6
omega4 = .05

aveRes = wanedList %>% 
  mutate(Delta = (omega1*V1 + omega2*(n-V1))/n,
         Omega = (omega3*V1 + omega4*(n-V1))/n) %>% 
  mutate(Day=1:n())
aveRes


aveRes %>%
  pivot_longer(2:3, names_to="Variant") %>% 
  ggplot(aes(x = Day)) +
  geom_line(aes(y = value, color=Variant)) +
  ylab(expression(paste("Average Resistance"))) +
  xlab("Days Since Max. Vaccine Effectiveness")+
  theme_bw() +
  scale_color_colorblind()
ggsave("plots/waning-resistance.png", width = 7, height = 3, 
        dpi = 600)
