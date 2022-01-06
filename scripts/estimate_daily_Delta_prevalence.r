library(tidyverse)
library(magrittr)

#Probability generating function of the serial intervall
G_w = function(x, w = generationInt()){
    nu = length(w)
    return(sum(w*x^(1:nu)))
}

#(Numerical) inverse of probability generating function of the serial intervall
G_w_inv = function(x, w = generationInt()){
    uniroot(f = function(b){x - G_w(x = b, w)}, interval = c(0.001,10))$root
}

#Estimation of growth rate given R(_eff)
estimate_beta = function(R = 1, w = generationInt()){
    return(1/G_w_inv(x = 1/R, w = w))
}

estimate_growth_rate = function(w, R){
#estimate_growth_rate = function(g, R){
    #w = generationInt(g$iTime, g$E, g$Var)
    return(1/G_w_inv(x = 1/R, w = w))
}

geometric_growth_model = function(I_t = 1000, b = 1, s = 1) {
    #I_(t-s) = b^(-s)I_t
    round(I_t*b^(-s))
}

### estimate exponential model given weekly data
weekly_lsp = function(b, I_t, J) {
    t_tilde = length(J)*7
    sum((colSums(matrix(round(b^(1:(t_tilde)-t_tilde) * I_t), nrow = 7, ncol = length(J))) - J)^2)
    }

KW_dates = as.Date("2021-01-10")+(0:51)*7
ages_variant_data = read.csv("data/ages_variant_data.csv")

ages_data_total = read.csv("data/ages_variant_data_new_base.csv") %>%
    select(KW, T) %>%
    filter(KW <= 45)

ages_data = 
ages_variant_data %>% 
mutate(W = round(W*ages_data_total$T),
       A = round(A*ages_data_total$T),
       B = round(B*ages_data_total$T),
       G = round(G*ages_data_total$T),
       D = round(D*ages_data_total$T),
       cases_weekly = ages_data_total$T
      )

KW_start = 22
KW_end = 27

KW_start_model = 20
KW_end_model = 30
dates_model = seq(KW_dates[KW_start_model]-6, KW_dates[KW_end_model], by = 1)

ages_data_ = ages_data %>% filter(KW >= KW_start, KW <= KW_end)

J_D = ages_data_ %>% pull(D)

pars_Delta = optim(par = c(1.1, tail(J_D,1)/7),
  fn = function(pars){
      weekly_lsp(b = pars[1], I_t = pars[2], J = J_D)
  }
 )$par
#print(pars_Delta)

J_A = ages_data_ %>% pull(A)
pars_Alpha = optim(par = c(1.1, tail(J_A,1)/7),
  fn = function(pars){
      weekly_lsp(b = pars[1], I_t = pars[2], J = J_A)
  }
 )$par

s = -as.integer(dates_model - KW_dates[KW_end])

model_data =     
tibble(
    date = dates_model,
    model_Delta = geometric_growth_model(I_t = pars_Delta[2], 
                   b = pars_Delta[1],
                   s = s),
    model_Alpha = geometric_growth_model(I_t = pars_Alpha[2], 
                   b = pars_Alpha[1],
                   s = s)
) %>% select(model_Delta, model_Alpha) %>%
apply(2, function(col){colSums(matrix(col, nrow = 7, ncol = KW_end_model - KW_start_model +1
                                     ))}) %>% 
    as.data.frame %>% 
mutate(cases_total = model_Delta + model_Alpha) %>%
mutate(ratio_Delta = model_Delta/cases_total,
       ratio_Alpha = model_Alpha/cases_total) %>%
mutate(date = KW_dates[KW_start_model:KW_end_model])

model_data_daily =  
tibble(
        date = dates_model,
        model_Delta = geometric_growth_model(I_t = pars_Delta[2], 
                       b = pars_Delta[1],
                       s = s),
        model_Alpha = geometric_growth_model(I_t = pars_Alpha[2], 
                       b = pars_Alpha[1],
                       s = s)
    ) %>% select(model_Delta, model_Alpha) %>% 
    mutate(cases_total = model_Delta + model_Alpha) %>%
    mutate(ratio_Delta = model_Delta/cases_total,
           ratio_Alpha = model_Alpha/cases_total,
           date = dates_model)

ages_variant_data %>% mutate(date = KW_dates[KW]) %>%
    filter(KW >= 20, KW <= 30) %>%
    select(date, D)



ages_variant_data %>% mutate(date = KW_dates[KW]) %>%
    filter(KW >= 20, KW <= 30) %>%
              pull(date, D) %>%
            mutate(what = "Observed Proportion (Weekly Average)")

prevalence_data = 
    rbind(ages_variant_data %>% mutate(date = KW_dates[KW]) %>%
    filter(KW >= 20, KW <= 30) %>%
              select(date, D) %>%
              rename(ratio_Delta = D) %>%
            mutate(what = "Observed Proportion (Weekly Average)"),
        model_data %>%
              select(date, ratio_Delta) %>% 
            mutate(what = "Model (Weekly Average)")
    ) %>% mutate(what = factor(what, 
                               levels = c("Observed Proportion (Weekly Average)", "Model (Weekly Average)" )))

model_data_daily %>%
    ggplot() +
    geom_line(
              aes(x = date, y = ratio_Delta, linetype = "Model (Daily)"), 
              colour = "#0072B2"#, alpha = 0.7
    ) +
    geom_point( data = prevalence_data, aes(x = date-3, y = ratio_Delta, 
               colour = what), alpha = 0.6) +
    geom_hline(yintercept = 0.2, alpha = 0.5) + 
    geom_vline(xintercept = as.Date("2021-06-12"), linetype = "dashed") +
    theme_bw() +
    theme(legend.position = "bottom",
          aspect.ratio = 0.3) +
    scale_x_date(date_breaks = "months", date_labels = "%b") +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2), lim = c(0, 1)) +
    scale_colour_manual(values = c("#CC79A7","#0072B2")) +
    ylab("Delta Prevalence") +
    xlab("Date") +
    labs(linetype = "",colour = "") +
    guides(shape = guide_legend(order = 2), 
           col = guide_legend(order = 1))

ggsave("Delta-prevalence-estimate.eps", device=cairo_ps,  width = 7, height = 3.5)


