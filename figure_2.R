## This file creates figure 2

#################
# Libraries
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(patchwork)
library(ggdist)

#################
# Read in data

LCC1_V <- read_csv('results/fitted_draws/LCC1_vehicle.csv') %>% 
  mutate(condition = 'LCC1, Vehicle') %>% 
  mutate(treatment = 'Vehicle')

LCC1_T <- read_csv('results/fitted_draws/LCC1_treatment.csv') %>% 
  mutate(condition = 'LCC1, Treatment') %>% 
  mutate(treatment = 'Treatment')

LCC9_V <- read_csv('results/fitted_draws/LCC9_vehicle.csv') %>% 
  mutate(condition = 'LCC9, Vehicle') %>% 
  mutate(treatment = 'Vehicle')

LCC9_T <- read_csv('results/fitted_draws/LCC9_treatment.csv') %>% 
  mutate(condition = 'LCC9, Treatment') %>% 
  mutate(treatment = 'Treatment')

df <- full_join(LCC1_V, LCC1_T) %>% 
  full_join(LCC9_V) %>% 
  full_join(LCC9_T) %>% 
  filter(.width == .95) %>% 
  mutate(treatment = as.factor(treatment)) %>% 
  mutate(treatment = relevel(treatment, ref = 'Vehicle'),
         condition = as.factor(condition))

##

LCC1_V <- read_csv('results/parameters/LCC1_vehicle.csv') %>%
  mutate(condition = 'LCC1, Vehicle') %>% 
  mutate(treatment = 'Vehicle')


LCC1_T <- read_csv('results/parameters/LCC1_treatment.csv') %>%
  mutate(condition = 'LCC1, Treatment') %>% 
  mutate(treatment = 'Treatment')

LCC9_V <- read_csv('results/parameters/LCC9_vehicle.csv') %>%
  mutate(condition = 'LCC9, Vehicle') %>% 
  mutate(treatment = 'Vehicle')

LCC9_T <- read_csv('results/parameters/LCC9_treatment.csv') %>%
  mutate(condition = 'LCC9, Treatment') %>% 
  mutate(treatment = 'Treatment')

df.2 <- full_join(LCC1_V, LCC1_T) %>% 
  full_join(LCC9_V) %>% 
  full_join(LCC9_T) %>% 
  mutate(i = if_else(i == 1, 'Intrinsic growth rate (r)',
                 'Carrying capacity (K)'))

#################
# Panel 1

pal <- c(brewer.pal(9, 'Blues')[c(4, 7)], brewer.pal(9, 'Greens')[c(4, 7)])

df$condition <- ordered(df$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))

p1 <- ggplot(df, aes(x = rank, group = condition))+
  geom_line(aes(y = y_rep_mean, color = condition), size = 1.25, alpha=0.75)+
  geom_point(aes(y = agg_count, color = condition), size = 1.75)+
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = interaction(condition, .width), fill = condition), alpha = .35, show.legend = F)+
  theme_minimal()+
  scale_color_manual(values = pal[c(1, 3, 2, 4)])+
  scale_fill_manual(values = pal[c(1, 3, 2, 4)])+
  labs(x = 'Time steps',
       y = 'Cell count',
       color = 'Condition',
       tag = 'A')+
  theme(legend.position = c(0.75, 0.85),
        legend.text=element_text(size=18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.title=element_text(size = 20),
        strip.text.x = element_text(size = 20))+
  facet_wrap(~condition)

#################
# Panel 2

df.2$i <- ordered(df.2$i, c('Intrinsic growth rate (r)', 'Carrying capacity (K)'))

p2 <- ggplot(df.2, aes(x = theta, y = condition, fill = condition))+
  stat_halfeye(orientation = 'horizontal',
               show.legend = F,
               normalize = 'xy',
               .width = c(.95, .4),
               aes(fill_ramp = stat(cut_cdf_qi(cdf, 
                                               .width = c(.5, .95, 1), 
                                               labels = scales::percent_format(accuracy = 1)))))+
  scale_fill_ramp_discrete(range = c(.9, .3), na.translate = F)+
  facet_wrap(~i, scales = 'free_x')+
  scale_y_discrete(limits = c('LCC9, Treatment','LCC9, Vehicle', 'LCC1, Treatment', 'LCC1, Vehicle'))+
  theme_minimal()+
  scale_fill_manual(values = pal[c(3, 1, 4, 2)])+
  labs(x = '', y = '', tag = 'B')+
  theme(strip.text.x = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15))


#################
# Save figure

ggsave('results/images/figures/fig2.jpeg', p1 | p2, scale = 2)


