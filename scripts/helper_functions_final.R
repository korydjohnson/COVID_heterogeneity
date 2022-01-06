###############################
# Helper Functions 
###############################


seasonality_day = function(date,
                           e = 0.4,
                           t_peak = as.Date("2021-01-01")) {
  time_to_peak = -as.numeric(date - t_peak)
  (1 - e * (cos(2 * pi * (
    time_to_peak / 365.25 + 1 / 2
  )) + 1) / 2)
}

compute_rolling_averages <- function(I, ra = 7) {
  #Given a time series "I" this function computes the rolling averages over "ra" many days.
  #The first (ra-1) days the rolling averages will be computed over resp. shorter periods.
  #OUTPUT: RA is a vector of same length as I.
  RA <- rep(NA, length(I))
  for (t in 1:(ra - 1)) { RA[t] <- sum(I[1:t]) / t }
  for (t in ra:length(I)) { RA[t] <- sum(I[(t - ra + 1):t]) / ra }
  
  RA
}



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
  
  group_sizes = sapply(groups, function(g) {
    g$S0
  })
  
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