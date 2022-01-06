library(tidyverse)
library(magrittr)
options(dplyr.summarise.inform = FALSE)
library(rvest)
library(foreach)
library(truncnorm)

source("scripts/simulation_final.R")

forecasting_date = as.Date("2021-08-08") # Date from which
nu = 13
R0 = 3.5
giantPlot=F

sim_end = as.Date("2022-05-01")
simLength = as.integer(sim_end - forecasting_date)
toBoost = F
toWane = F
maximum_vaccine_uptake = 0.85 # percentage of whole population which will get a vaccine


### import current population makeup of Austria
new_samples = T
n_samples = 5
file_suffix = ""

# Colours
colour_palette = c("#E69F00",
                   "#56B4E9",
                   "#009E73",
                   "#F0E442",
                   "#0072B2",
                   "#D55E00",
                   "#CC79A7")


# variants to be considered in the model
variant_names = strain_names = c("W", "A", "B", "G", "D", "O")
variant_names_full = c("WT", "Alpha", "Beta", "Gamma", "Delta", "Omega")
# varinats to plot
variant_names_selection = c("A", "B", "G", "D", "O")
variant_names_selection_full = c("Alpha", "Beta", "Gamma", "Delta", "Omega")

# colour selection
colour_names = c(variant_names_full, "total")
colours = set_names(as.list(colour_palette), colour_names)
colour_palette = unlist(colours[variant_names_selection_full])

vaccine_names = c("mRNA", "vector")

vaccination_plan_with_past_data = F # =T -> 'future' vaccination data will be used for the model.
# (irrelevant for these plots but needed for calling 'multi_variant_population_parameters_imports.r')
source("scripts/multi_variant_population_parameters_imports.r")
# Import and process case data, variant data and vaccination data

# Import resistance data (formerly from googledoc)
resistance_data = read.csv("data/resistance_parameters.csv", 
                           stringsAsFactors = F) %>%
  select(c("name", "type", "bound", all_of(variant_names))) %>%
  filter((type != "infection") | (name %in% variant_names))

# Draw variant samples
if (new_samples == T) {
  tic = Sys.time()
  variant_data_samples =
    foreach(i = 1:n_samples, .combine = rbind) %do% {
      resistance_data_sample = sample_variant_data_truncnormal(resistance_data, variant_names) %>% mutate(sample = i)
    }
  print(Sys.time() - tic)
}

group_names = c("_",
                variant_names,
                foreach(vaccine = vaccine_names, .combine = c) %do% {
                  c(vaccine, paste0(variant_names, vaccine))
                })
r = seq(0, 1, 0.01)
#ratio_mRNA = vaccination_group_sizes$mRNA / (vaccination_group_sizes$mRNA + vaccination_group_sizes$vector)
#ratio_vector = 1 - ratio_mRNA
ratio_mRNA = 0.8
ratio_vector = 1 - ratio_mRNA

#
#
#
#
#
### vary only over ratio of previously infected, keep mRNA ratio fixed:
if (giantPlot) {
  
  variant_dataset = foreach(infection_ratio = c(0, 0.2, 0.4), .combine = rbind) %do% {
    infection_ratio = infection_ratio
    
    population = c(list("_" = total_uninfected),
                   S0_unvaccinated) ## population before vaccination
    infection_ratio_aut = sum(unlist(population)[2:length(population)]) / N
    
    population = unlist(population[c("_", variant_names)]) / N
    population[1] = population[1] * (1 - infection_ratio) / (1 - infection_ratio_aut)
    population[2:length(population)] = population[2:length(population)] * infection_ratio /
      infection_ratio_aut
    
    data =
      cbind(
        (1 - r) %>% map_dfr(~ . * population),
        r %>% map_dfr(~ . * population * ratio_mRNA),
        r %>% map_dfr(~ . * population * ratio_vector)
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
          
          rho_t = Gamma_tilde * lambda * R0 #/ N
          
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
    
    variant_data %>%
      group_by(r, variant) %>%
      summarise(quant = list(quantile(
        rho_t, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
      ))) %>%
      unnest_wider(quant) %>%
      mutate(variant = factor(variant, levels =  variant_names_full)) %>%
      mutate(previous_infections = paste0(infection_ratio * 100, "% previous infections"))
  }
  
  ### PLOT DATA
  variant_dataset %>%
    ggplot(aes(x = r)) +
    geom_vline(aes(xintercept = 0.5), linetype = "dashed", alpha = 0.5) +
    geom_line(aes(y = `50%`, colour = variant))  +
    geom_ribbon(aes(ymin = `25%`,
                    ymax = `75%`,
                    fill = variant), alpha = .4) +
    geom_ribbon(aes(ymin = `2.5%`,
                    ymax = `97.5%`,
                    fill = variant), alpha = .25) +
    theme_bw() +
    theme(legend.position = "bottom",
          aspect.ratio = 0.5) +
    scale_colour_manual(values = colour_palette) +
    scale_fill_manual(values = colour_palette) +
    xlab("Fully Vaccinated") +
    ylab(expression(paste(rho[t] ^ V))) +
    labs(colour = "Variant", fill = "Variant") +
    xlim(c(0, 1)) +
    facet_wrap(vars(previous_infections), nrow = 3)
  
  ### SAVE PLOT
  ggsave(
    paste0("plots/variant-rho-quantiles-truncnormal-test",
           file_suffix,
           ".eps"),
    width = 7,
    height = 7, device=cairo_ps, dpi=600
  )
}

#
#
#
#
#
### vary over ratio of previously infected AND mRNA ratio:
# mRNA_ratio = c(0, 0.25, 0.5, 0.75, 1)
# blah = c(0, 0.2, 0.4)

variant_dataset =
  foreach(mRNA_ratio = c(0, 1), .combine = rbind) %do% {
    ratio_mRNA = mRNA_ratio
    ratio_vector = 1 - ratio_mRNA
    
    foreach(blah = c(0, 0.4), .combine = rbind) %do% {
      infection_ratio = blah
      
      population = c(list("_" = total_uninfected),
                     S0_unvaccinated) ## population before vaccination
      infection_ratio_aut = sum(unlist(population)[2:length(population)]) /
        N
      
      population = unlist(population[c("_", variant_names)]) / N
      population[1] = population[1] * (1 - infection_ratio) / (1 - infection_ratio_aut)
      population[2:length(population)] = population[2:length(population)] *
        infection_ratio / infection_ratio_aut
      
      data =
        cbind(
          (1 - r) %>% map_dfr( ~ . * population),
          r %>% map_dfr( ~ . * population * ratio_mRNA),
          r %>% map_dfr( ~ . * population * ratio_vector)
        ) %>% set_names(group_names)
      
      variant_data =
        foreach(i = 1:max(variant_data_samples$sample),
                .combine = rbind) %do% {
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
                    
                    rho_t = Gamma_tilde * lambda * R0 #/ N
                    
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
      
      variant_data %>%
        group_by(r, variant) %>%
        summarise(quant = list(quantile(
          rho_t, probs = c(0.025, 0.25, 0.5, 0.75, 0.975)
        ))) %>%
        unnest_wider(quant) %>%
        mutate(variant = factor(variant, levels =  variant_names_full)) %>%
        mutate(previous_infections = paste0(blah * 100, "% prev. infections"))
      
    } %>% mutate(
      vaccine_ratio = paste0(
        ratio_mRNA * 100,
        "% mRNA vaccines, ",
        ratio_vector * 100,
        "% vector vaccines"
      )
    )
  }

### PLOT DATA
variant_dataset %>%
  mutate(
    previous_infections = factor(previous_infections, levels = unique(previous_infections)),
    vaccine_ratio = factor(vaccine_ratio, levels = unique(vaccine_ratio))
  ) %>%
  ggplot(aes(x = r)) +
  # geom_vline(aes(xintercept = 0.5), linetype = "dashed", alpha = 0.5) +
  geom_line(aes(y = `50%`, colour = variant))  +
  geom_ribbon(aes(ymin = `25%`,
                  ymax = `75%`,
                  fill = variant), alpha = .4) +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`,
                  fill = variant), alpha = .25) +
  theme_bw() +
  theme(legend.position = "bottom",
        aspect.ratio = 0.5) +
  scale_colour_manual(values = colour_palette) +
  scale_fill_manual(values = colour_palette) +
  xlab("Fully Vaccinated") +
  ylab(expression(paste(rho[t] ^ V))) +
  labs(colour = "Variant", fill = "Variant") +
  xlim(c(0, 1)) +
  #facet_grid(rows = vars(vaccine_ratio), cols = vars(previous_infections))
  facet_grid(cols = vars(vaccine_ratio),
             rows = vars(previous_infections))
#facet_wrap(vaccine_ratio ~ previous_infections, ncol = 3)

### SAVE PLOT
ggsave(
  paste0("plots/variant-rho-quantiles-truncnormal-big-test",
         file_suffix,
         ".eps"),
  width = 10,
  height = 7, device=cairo_ps, dpi=600
)
