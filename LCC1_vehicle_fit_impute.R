
library(tidyverse)
library(rstan)
library(lubridate)
library(tidybayes)
library(ggdist)

############
# Stan setup

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
stan_file_path <- 'stan_model/logistic_impute_log.stan'

############
# Load and set up data

df <- read_csv("data/green_count_clean.csv")
#########
# Process data

df <- df %>% 
  select(well, field, rank, count, date) %>% 
  filter(well == 'A1' | well == 'A2') %>% 
  group_by(well, rank) %>% 
  mutate(agg_count = sum(count)) %>% 
  select(-count)

#########
# Add missingness

df.full <- df %>% 
  group_by(well, field) %>% 
  complete(date = seq(from = ymd_hms("2019-12-06 16:00:00 UTC"), # adds datetimes every two hours
                      to = ymd_hms("2019-12-19 14:00:00 UTC"),
                      by = '2 hours')) %>% 
  mutate(rank = rank(date), # recalculate ranks
         interpolate = if_else(is.na(agg_count), T, F))

########
# Take average of replicates

df.process <- df.full %>% 
  ungroup() %>% 
  group_by(rank) %>% 
  filter(rank > 24) %>% 
  mutate(rank = rank - 24) %>% 
  select(rank, agg_count, well) %>% 
  mutate(well = as.numeric(as.factor(well))) %>% 
  unique()

#########
# Preview data

ggplot(df.process, aes(x = rank, y = agg_count, color = as.factor(well)))+
  geom_point()+
  theme_bw()

########
# Set up data for Stan

t_obs <- df.process %>% 
  filter(well == 1) %>% 
  filter(!is.na(agg_count)) %>% 
  pull(rank)

t_mis <- df.process %>% 
  filter(well == 1) %>% 
  filter(is.na(agg_count)) %>% 
  pull(rank)

ts = df.process %>% 
  filter(well == 1) %>% 
  pull(rank)

y_obs <- matrix(c(df.process %>% 
  filter(!is.na(agg_count)) %>% 
  unique() %>% 
  pivot_wider(id_cols = rank,
              names_from = well,
              values_from = agg_count) %>% 
  pull(`1`), 
  df.process %>% 
    filter(!is.na(agg_count)) %>% 
    unique() %>% 
    pivot_wider(id_cols = rank,
                names_from = well,
                values_from = agg_count) %>% 
    pull(`2`)),
  ncol = 2,
  byrow = F)
  

dat <- list(t_obs = t_obs,
            t_mis = t_mis,
            N_obs = length(t_obs),
            N_mis = length(t_mis),
            R = 2,
            ts = ts,
            y0 = y_obs[1,],
            y_obs = y_obs,
            compute_likelihood = 1)
##########
# Compile model
model <- stan_model(file = stan_file_path,
                    model_name = 'LCC1_logistic_growth_impute')
##########
# Fit the model

# specify parameters of interest
pars <- c("theta", "sigma", 'y0_hat', 'y_rep')

# fit the mode
fit <- sampling(object = model,
                pars = pars,
                data = dat,
                chains = 4,
                iter = 10000,
                sample_file = 'stan_samples/LCC1_vehicle.csv',
                diagnostic_file = 'stan_diagnostics/LCC1_treatment.csv',
                seed = 12345)  

trace <- traceplot(fit, pars = c('theta', 'sigma', 'y0_hat')) 
ggsave('results/images/traceplot/LCC1_cehicle.jpeg', plot = trace)

##########
# Visualize PPC

t <- fit %>% 
  spread_draws(y_rep[i, j]) %>% 
  curve_interval(y_rep, .width = .95) %>% 
  rename(rank = i) %>% 
  left_join(df.process) 

s <- fit %>% 
  spread_draws(y_rep[i, j]) %>% 
  mean_qi(y_rep, na.rm = T) %>% 
  rename(rank = i,
         y_rep_mean = y_rep) %>% 
  select(rank, y_rep_mean)
t <- full_join(s, t)

write_csv(t, 'results/fitted_draws/LCC1_vehicle.csv')


ggplot(t, aes(x = rank))+
  geom_line(aes(x = rank, y = y_rep_mean))+
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = .15)+
  geom_point(aes(rank, agg_count))+
  theme_bw()+
  labs(x = 'Sampling time point', y = 'Cells per well')
ggsave('results/images/ppd/LCC1_vehicle_fit.jpeg', dpi = 300)

# Print summary stats
print(fit, pars = c("theta", "sigma", 'y0_hat'))

# Save summary stats
write_csv(summary(fit)$summary %>% as_tibble(rownames = 'par'), 'results/parameters/summary_statistics/LCC1_vehicle.csv')

# Save raw parameter ppd draws
fit %>% 
  spread_draws(theta[i]) %>% 
  write_csv('results/parameters/LCC1_vehicle.csv')

