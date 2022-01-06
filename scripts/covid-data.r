library(tidyverse)
library(magrittr)
library(foreach)
options(dplyr.summarise.inform = FALSE)
library(rvest)

data_cutoff = as.Date("2021-12-15")


# R_EFF ESTIMATES ----------------------------------
### possibly not the R_eff estimates that will be used fot the paper
from_website = F
if(from_website){
    #Data from website
    EpiNow_R_estimates =
      readRDS(url("https://nowcasting.manulari.eu/np374bw/national/Austria/latest/bigr_estimates.rds", "rb")) %>% 
        select(date, median) %>%
        filter(date <= data_cutoff)
} else{
    #Data from file
    EpiNow_R_estimates = readRDS("bigr_estimates.rds", "rb") %>% 
        select(date, median) %>%
        filter(date <= data_cutoff)
}
write.csv(EpiNow_R_estimates, file = ("EpiNow_R_estimates.csv"),row.names = F)


# AGES VOC ----------------------------------
VOC_names = c("Alpha", "Beta", "Gamma", "Delta")
VOC_names_short = c("A", "B", "G", "D")

#Week 1-9
url = "https://web.archive.org/web/20210506123009/https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/"
mutation_data_aut_hist = (read_html(url) %>% html_table(fill = TRUE))[[2]][,-1] %>% 
    apply(.,2, function(chr){parse_number(gsub(",", ".",chr))/100}) %>% 
    t %>% as.data.frame
names(mutation_data_aut_hist) = c("Wild", "Alpha", "AlphaE","Beta", "Gamma", "VOC_not_specified")
mutation_data_aut_hist[is.na(mutation_data_aut_hist)] = 0
mutation_data_aut_hist %<>% 
    mutate(Alpha = Alpha+AlphaE) %>% 
    select(-AlphaE) %>%
    mutate(VOC_specified = Wild + Alpha + Beta + Gamma) %>%
    mutate(W = Wild + Wild/VOC_specified*VOC_not_specified,
           A = Alpha + Alpha/VOC_specified*VOC_not_specified,
           B = Beta + Beta/VOC_specified*VOC_not_specified,
           G = Gamma + Gamma/VOC_specified*VOC_not_specified,
          ) %>%
    select(W, A, B, G) %>%
    mutate(KW = 1:14, .before = 1)

ages_variant_data_hist = mutation_data_aut_hist

#week 10-12
url1 = "https://web.archive.org/web/20210506123009/https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/"
VOC_ratio1 = ((read_html(url1) %>% html_table(fill = TRUE))[[1]][10,-1] %>% 
    apply(.,2, function(chr){parse_number(gsub(",", ".",chr))/100}))[paste0("KW", 10:12)]
#week 13-16(15?) (VOC ratio) and week 10-16 for VOC differentiation
url2 = "https://web.archive.org/web/20210609161726/https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/"
VOC_ratio2 = 
((read_html(url2) %>% html_table(fill = TRUE))[[1]][10,-1] %>% 
    apply(.,2, function(chr){parse_number(gsub(",", ".",chr))/100}))[paste0("KW", 13:16)]

VOC_ratio = c(VOC_ratio1, VOC_ratio2)

mutation_data_aut = (read_html(url2) %>% html_table(fill = TRUE))[[2]] %>% 
    select(Varianten, paste0("KW", 10:16)) %>%
    .[c(-3, -5),]  %>% 
    mutate(variant = c("VOC_total", rep(VOC_names, each = 2))) %>%
    select(-Varianten) %>%
    group_by(variant) %>% 
    summarise_at(vars(paste0("KW", 10:16)),sum)

VOC_total = mutation_data_aut %>% filter(variant == "VOC_total") %>% select(-variant) %>% unlist

VOC_specified = mutation_data_aut %>% 
    filter(variant %in% VOC_names) %>% 
    select(-variant) %>% colSums

VOC_non_specified = VOC_total - VOC_specified

ages_variant_data_mid = 
foreach(VOC = VOC_names, .combine = cbind) %do% {
    ratios = 
    (mutation_data_aut  %>% 
         filter(variant == VOC) %>% 
         select(-variant) %>% 
         unlist)*(1+VOC_non_specified/VOC_specified)/VOC_total* VOC_ratio 

    tibble(ratios)
} %>% set_names(VOC_names_short) %>% 
    mutate(W = 1-VOC_ratio, .before = 1) %>%
    mutate(KW = 10:16, .before = 1) #rows sum up to 1, this is nice :)

## DATA KW 17 to KW 45
ages_url = "https://web.archive.org/web/20211208215839/https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/"
ages_data = read_html(ages_url) %>% html_table(fill = TRUE)
ages_data = ages_data[[1]][1:(nrow(ages_data[[1]])-1),]

ages_variant_data_new  = ages_data %>%
    mutate(VOC_total = `B.1.1.7 (Alpha)` + `B.1.351 (Beta)` + `P.1 (Gamma)` + `B.1.617.2 (Delta)`) %>%
    rename(cases_weekly = `Fälle gesamt`,
           KW = Kalenderwoche
          ) %>%
    mutate(Alpha_lower = `B.1.1.7 (Alpha)`/cases_weekly,
           Beta_lower = `B.1.351 (Beta)`/cases_weekly,
           Gamma_lower = `P.1 (Gamma)`/cases_weekly,
           Delta_lower = `B.1.617.2 (Delta)`/cases_weekly
          ) %>%
    mutate(Alpha_upper = `B.1.1.7 (Alpha)`/VOC_total,
           Beta_upper = `B.1.351 (Beta)`/VOC_total,
           Gamma_upper = `P.1 (Gamma)`/VOC_total,
           Delta_upper = `B.1.617.2 (Delta)`/VOC_total
          ) %>%
    mutate(VOC_upper = 1) %>%
    #mutate(VOC_lower = VOC_total/cases_weekly) %>% mutate(non_Delta = `B.1.1.7 (Alpha)` + `B.1.351 (Beta)` + `P.1 (Gamma)`) %>%
    #mutate(new_upper_Delta = 1-non_Delta/cases_weekly)
    rename(A = Alpha_upper,
           B = Beta_upper,
           G = Gamma_upper,
           D = Delta_upper
    ) %>% select(KW, A, B, G, D) %>% mutate(W = 0, .before = 2)

ages_variant_data = 
    rbind(ages_variant_data_hist %>% mutate(D = 0) %>% filter(KW <= 9),
          ages_variant_data_mid,
          ages_variant_data_new %>% filter(KW > 16)
)

# WRITE CSV
write.csv(ages_variant_data, file = "ages_variant_data.csv", row.names = F)

#ages_url = "https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/"


# VOC take 2 --------------------------------------------------------------

# read.csv("data/ages_variant_data_new_base.csv") %>%
#   mutate(VOCt = A+B+G+D+O) %>%
#   mutate_at(vars(2:6), ~./`T`)

# WHO CASE NUMBERS ----------------------------------
cases_who = read.csv("https://covid19.who.int/WHO-COVID-19-global-data.csv") %>%
  filter(Country == "Austria") %>%
  mutate(date = Date_reported %>% as.Date,
         cases = New_cases) %>%
  select(date, cases) %>%
  filter(date <= data_cutoff)

write.csv(cases_who, file = ("cases_who.csv"),row.names = F)

# VACCINATION DATA ----------------------------------
# former dataset used for initial submission (dataset will be retired mid December 2021):
#vaccination_data_aut = read.csv("https://info.gesundheitsministerium.gv.at/data/timeline-eimpfpass.csv",
#                                sep = ";") %>%
#    filter(Name == "Österreich") %>%
#    mutate(date = as.Date(substr(Datum, 1, 10))) %>%
#    rename(
#        BionNtechPfizer = EingetrageneImpfungenBioNTechPfizer_1,
#        Moderna = EingetrageneImpfungenModerna_1,
#        AstraZeneca = EingetrageneImpfungenAstraZeneca_1,
#        Janssen = EingetrageneImpfungenJanssen
#    ) %>%
#    select(date, BionNtechPfizer, Moderna, AstraZeneca, Janssen) %>%
#    filter(date <= "2021-11-30")

vaccination_data_aut = 
read.csv("https://info.gesundheitsministerium.gv.at/data/COVID19_vaccination_doses_timeline.csv", sep = ";") %>%
    filter(dose_number == 1, state_name == "Österreich") %>%
    mutate(date = as.Date(substr(date, 1, 10))) %>%
    select(date, vaccine, doses_administered_cumulative) %>%
    pivot_wider(names_from = vaccine, values_from = doses_administered_cumulative) %>%
    filter(date <= data_cutoff)
#vaccination_data_aut
write.csv(vaccination_data_aut, file = ("data/vaccination_data_aut.csv"),row.names = F)

# boosters
booster_data_aut = read.csv("data/COVID19_vaccination_doses_timeline.csv", sep = ",") %>%
  filter(dose_number == 3, state_name == "Österreich") %>%
  mutate(date = as.Date(substr(date, 1, 10))) %>%
  select(date, vaccine, doses_administered_cumulative) %>%
  group_by(date) %>% 
  summarise(cumulativeDoses = sum(doses_administered_cumulative))
write.csv(booster_data_aut, file = ("data/booster_data_aut.csv"),row.names = F)


