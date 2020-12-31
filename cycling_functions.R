library(tidyverse)
library(deSolve)

#############
# Data

## VEHICLE

s_mono_fit <- read_csv('results/parameters/summary_statistics/LCC1_vehicle.csv')
r_mono_fit <- read_csv('results/parameters/summary_statistics/LCC9_vehicle.csv')


LV_fit <- read_csv('results/parameters/summary_statistics/LV_vehicle.csv')

K_s <- s_mono_fit %>% 
  filter(par == 'theta[2]') %>% 
  select(mean, sd)


r_s <- s_mono_fit %>% 
  filter(par == 'theta[1]') %>% 
  select(mean, sd)


K_r <- r_mono_fit %>% 
  filter(par == 'theta[2]') %>% 
  select(mean, sd)


r_r <- r_mono_fit %>% 
  filter(par == 'theta[1]') %>% 
  select(mean, sd)

S0 <- LV_fit %>% 
  filter(par == 'y0_hat[1]') %>% 
  select(mean, sd)

R0 <- LV_fit %>% 
  filter(par == 'y0_hat[2]') %>% 
  select(mean, sd)

alpha <- LV_fit %>% 
  filter(par == 'theta[1]') %>% 
  select(mean, sd)

beta <- LV_fit %>% 
  filter(par == 'theta[2]') %>% 
  select(mean, sd)

sigma <- LV_fit %>% 
  filter(par == 'sigma') %>% 
  select(mean, sd)

pars_veh <- c(r_s = r_s,
              K_s = K_s,
              r_r = r_r,
              K_r = K_r,
              alpha = alpha,
              beta = beta)

### TREAT

s_mono_fit <- read_csv('results/parameters/summary_statistics/LCC1_treatment.csv')
r_mono_fit <- read_csv('results/parameters/summary_statistics/LCC9_treatment.csv')


LV_fit <- read_csv('results/parameters/summary_statistics/LV_treatment.csv')

K_s <- s_mono_fit %>% 
  filter(par == 'theta[2]') %>% 
  select(mean, sd)


r_s <- s_mono_fit %>% 
  filter(par == 'theta[1]') %>% 
  select(mean, sd)


K_r <- r_mono_fit %>% 
  filter(par == 'theta[2]') %>% 
  select(mean, sd)


r_r <- r_mono_fit %>% 
  filter(par == 'theta[1]') %>% 
  select(mean, sd)

S0 <- LV_fit %>% 
  filter(par == 'y0_hat[1]') %>% 
  select(mean, sd)

R0 <- LV_fit %>% 
  filter(par == 'y0_hat[2]') %>% 
  select(mean, sd)

alpha <- LV_fit %>% 
  filter(par == 'theta[1]') %>% 
  select(mean, sd)

beta <- LV_fit %>% 
  filter(par == 'theta[2]') %>% 
  select(mean, sd)

sigma <- LV_fit %>% 
  filter(par == 'sigma') %>% 
  select(mean, sd)

pars_treat <- c(r_s = r_s,
              K_s = K_s,
              r_r = r_r,
              K_r = K_r,
              alpha = alpha,
              beta = beta)

##############
# Define function

LV <- function(Time, State, parameters){
  
    r_s <- parameters[1]
    K_s <- parameters[2]
    
    r_r <- parameters[3]
    K_r <- parameters[4]
    
    alpha <- parameters[5]
    beta <- parameters[6]
    
    S = State[1]
    R = State[2]
    
    
    dSdt <- r_s * S * (1 - (S / K_s) - ((alpha * R) / K_s))
    
    dRdt <- r_r * R * (1 - (R / K_r) - ((beta * S) / K_r))
    
    list(c(dSdt = dSdt,
                dRdt = dRdt))
}

############

treat_cycle <- function(t_treat, t_veh, cycles, y0, pars_treat, pars_veh){
  # Function cycles treatment with n cycles with treatment length t_treat, vehicle length t_vehicle,
  # initial condition y0, vehicle parameterization pars_veh, and treatment parameterization pars_tre
  
  ####
  # instantiate vectors to hold pops
  
  S_pop <- vector(mode = 'double', length = (t_treat + t_veh) * cycles)
  R_pop <- vector(mode = 'double', length = (t_treat + t_veh) * cycles)
  
  ####
  
  # Move through the cycles
  for (i in 1:cycles){
      
      # increment population vectors
      incr <- (t_treat + t_veh) * (i-1) 
      
      # run 
      treat_out <- ode(y0, 1:t_treat, LV, pars_treat)
      veh_out <- ode(tail(treat_out, n = 1)[c(2, 3)], 1:t_veh, LV, pars_veh) # use last value of treat as y0 for veh
      
      S_pop[(incr+1):(incr + t_treat + t_veh)] <- c(treat_out[,2], veh_out[,2])
      R_pop[(incr+1):(incr + t_treat + t_veh)] <- c(treat_out[,3], veh_out[,3])
      
      y0 <- tail(veh_out, n = 1)[c(2, 3)]
      
  }
  
  list(S = S_pop,
       R = R_pop)
}

############

treat_cycle.rep <- function(nrep, t_treat, t_veh, cycles, y0, pars_treat, pars_veh){
  # Function runs treatment_cycle nrep times and stores results in a nrep x pop x time array
  
  ####
  # pre-allocate array
  pop <- tibble(S = numeric(),
                R = numeric(),
                t = numeric(), 
                rep = numeric(),
                condition = character())
  
  ####
  # generate random deviates
  pars_veh_full <- matrix(data = c(rnorm(nrep, pars_veh[[1]], pars_veh[[2]]),
                                     rnorm(nrep, pars_veh[[3]], pars_veh[[4]]),
                                     rnorm(nrep, pars_veh[[5]], pars_veh[[6]]),
                                     rnorm(nrep, pars_veh[[7]], pars_veh[[8]]),
                                     rnorm(nrep, pars_veh[[9]], pars_veh[[10]]),
                                     rnorm(nrep, pars_veh[[11]], pars_veh[[12]])),
                            nrow = nrep,
                            ncol = 6,
                            byrow = F)
  
  
  pars_treat_full <- matrix(data = c(rnorm(nrep, pars_treat[[1]], pars_treat[[2]]),
                                     rnorm(nrep, pars_treat[[3]], pars_treat[[4]]),
                                     rnorm(nrep, pars_treat[[5]], pars_treat[[6]]),
                                     rnorm(nrep, pars_treat[[7]], pars_treat[[8]]),
                                     rnorm(nrep, pars_treat[[9]], pars_treat[[10]]),
                                     rnorm(nrep, pars_treat[[11]], pars_treat[[12]])),
                            nrow = nrep,
                            ncol = 6,
                            byrow = F)
  
  ####
  # simulate treatment cycles
  
  for (rep in 1:nrep){
    
    p <- treat_cycle(t_treat, t_veh, cycles, y0, pars_treat_full[rep,], pars_veh_full[rep,])
    
    
    p <- p %>% 
      as_tibble() %>% 
      mutate(rep = rep,
             t = row_number(),
             treatment = rep(c(rep('Treatment', t_treat), rep('Vehicle', t_veh)), cycles))
    
    # infefficient hack but easier to use later w/ ggplot
    pop <- suppressMessages(full_join(pop, p))
  }
  
  # return array in long format
  fit <- pop %>% 
    pivot_longer(c(S, R), names_to = 'pop', values_to = 'N') %>% 
    mutate(pop = if_else(pop == 'S', 'LCC1', 'LCC9')) %>% 
    mutate(condition = paste(pop, treatment, sep = ', '))
  
  
  ####
  return(fit)
}

############


