# This file fits a Bayesian LV regression model

############
# Libraries

library(tidyverse)
library(rstan)
library(lubridate)
library(tidybayes)

############
# Stan setup

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
stan_file_path <- 'stan_model/LV_impute_log_full.stan'

###########
# Load data

sus <- read_csv('data/LCC1_coculture_clean.csv')
res <- read_csv('data/LCC9_coculture_clean.csv')


s_mono_fit <- read_csv('results/parameters/summary_statistics/LCC1_treatment.csv')
r_mono_fit <- read_csv('results/parameters/summary_statistics/LCC9_treatment.csv')


df <- full_join(res, sus)
###########
# Process data

df <- df %>% 
  dplyr::select(well, field, count, date, cell_type) %>% 
  filter(well == 'C3' | well == 'C4') %>% 
  group_by(well, date, cell_type) %>% 
  mutate(agg_count = sum(count)) %>% 
  dplyr::select(-count)

# Find date at which decline in well C4 begins
date.change <- df %>% 
  filter(well == 'C4') %>% 
  group_by(well) %>% 
  top_n(n=1) %>% 
  select(-field) %>% 
  unique() %>% 
  pull(date)

# Remove cell counts after decline
df <- df %>% 
  filter(date < date.change)

ggplot(df, aes(x = date, y = agg_count, color = well))+
  geom_point()

#########
# Add missingness

df.full <- df %>% 
  group_by(well, field, cell_type) %>% 
  complete(date = seq(from = ymd_hms(min(df$date)), # adds datetimes every two hours, "2019-12-06 16:00:00 UTC"
                      to = ymd_hms(date.change),
                      by = '2 hours')) %>% 
  mutate(rank = rank(date), # recalculate ranks
         interpolate = if_else(is.na(agg_count), T, F)) %>% 
  filter(field == 1) %>% 
  ungroup() %>% 
  select(well, cell_type, rank, agg_count) %>% 
  mutate(well = as.numeric(as.factor(well))) %>% 
  unique()

ggplot(df.full, aes(x = rank, y = agg_count, group = well))+
  geom_point()

########
# Set up data for Stan

t_obs <- df.full %>%
  filter(!is.na(agg_count),
         well == 1,
         cell_type == 'LCC1') %>% 
  pull(rank)

t_mis <- df.full %>% 
  filter(is.na(agg_count),
         well == 1,
         cell_type == 'LCC1') %>% 
  pull(rank)

ts = df.full %>%
  filter(well == 1,
         cell_type == 'LCC1') %>% 
  pull(rank)

y_obs = array(data = c(matrix(data = c(df.full %>% 
                                         filter(!is.na(agg_count),
                                                well == 1) %>% 
                                         pull(agg_count)),
                              ncol = 2),
                       matrix(data = c(df.full %>% 
                                         filter(!is.na(agg_count),
                                                well == 2) %>% 
                                         pull(agg_count)),
                              ncol = 2)), 
              dim = c(length(t_obs), 2, 2),
              dimnames = NULL)

y0 <- matrix(y_obs[1,1:2,],
             nrow = 2,
             byrow = T)


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

logistic_pars = c(r_s, r_r, K_s, K_r)

#######
# Specify inits

init <- list(r_s = r_s[[1]], 
             r_r = r_r[[1]], 
             K_s = K_s[[1]], 
             K_r = K_r[[1]],
             alpha = 0,
             beta = 0,
             sigma = array(.05),
             y0_hat  = rowMeans(y0))


inits <- list(init, init, init, init)


#######

compute_likelihood = 1


dat <- list(N_obs = length(t_obs),
            N_mis = length(t_mis),
            ts = ts,
            y0 = y0,
            R = 2,
            y_obs = y_obs,
            t_obs = t_obs,
            compute_likelihood = compute_likelihood,
            logistic_pars = logistic_pars)

#######
# Compile model

model <- stan_model(file = stan_file_path,
                    model_name = 'LV_impute_log')

#######
# Sample from model

fit <- sampling(object = model,
                pars = c('theta', 'sigma',  'y0_hat', 'y_rep'),
                data = dat,
                chains = 4,
                iter = 10000,
                sample_file = 'stan_samples/LV_treatment.csv',
                diagnostic_file = 'stan_diagnostics/LV_treatment.csv',
                seed = 42,
                control = list(max_treedepth = 15),
                init = inits)

trace <- traceplot(fit, pars = c('theta', 'sigma', 'y0_hat')) 
ggsave('results/images/traceplot/LV_treatment.jpeg', plot = trace)


##########

# Get 95% curvewise CIs
t <- fit %>% 
  spread_draws(y_rep[i, j, k]) %>% 
  ungroup() %>% 
  group_by(i, j) %>% 
  curve_interval(y_rep, .width = .95) %>% 
  mutate(cell_type = if_else(j == 1, 'LCC1', 'LCC9')) %>% 
  rename(rank = i) %>% 
  full_join(df.full) 

# But replace median curve with mean curve for smoothnness
s <- fit %>% 
  spread_draws(y_rep[i, j]) %>% 
  ungroup() %>% 
  group_by(i, j) %>% 
  mean_qi() %>% 
  select(i, j, y_rep) %>% 
  rename(rank = i,
         y_rep_mean = y_rep)

t <- full_join(s, t)


write_csv(t, 'results/fitted_draws/LV_treatment.csv')

ggplot(t, aes(x = rank, group = cell_type))+
  geom_point(aes(y = agg_count, color = cell_type))+
  geom_line(aes(y = y_rep_mean))+
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = interaction(.width, cell_type)), alpha = .1)+
  theme_bw()+
  labs(x = 'Time steps',
       y = 'Mean cell count',
       title = 'Coculture fit w/ treatment (with 95% CI)',
       color = 'Cell type')

ggsave('results/images/ppd/LV_treatment.jpeg')


# Save posterior output
write_csv(summary(fit)$summary %>% as_tibble(rownames = 'par'), 'results/parameters/summary_statistics/LV_treatment.csv')

# save interaction raw data
fit %>% 
  spread_draws(theta[i]) %>% 
  write_csv('results/parameters/LV_treatment.csv')

# Save phase space draws
fit %>% 
  spread_draws(y_rep[i, j]) %>% 
  pivot_wider(id_cols = c(i, .draw),
              names_from = j,
              values_from = y_rep) %>% 
  filter(.draw %in% seq(1, 150, 1)) %>% 
  write_csv('results/phase_space_draws/LV_treatment.csv')


