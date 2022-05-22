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
  mutate(condition = 'LCC1 (SENS), VEH') %>% 
  mutate(treatment = 'Vehicle')

LCC1_T <- read_csv('results/fitted_draws/LCC1_treatment.csv') %>% 
  mutate(condition = 'LCC1 (SENS), TRT') %>% 
  mutate(treatment = 'Treatment')

LCC9_V <- read_csv('results/fitted_draws/LCC9_vehicle.csv') %>% 
  mutate(condition = 'LCC9 (RES), VEH') %>% 
  mutate(treatment = 'Vehicle')

LCC9_T <- read_csv('results/fitted_draws/LCC9_treatment.csv') %>% 
  mutate(condition = 'LCC9 (RES), TRT') %>% 
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
  mutate(condition = 'LCC1 (SENS), VEH') %>% 
  mutate(treatment = 'Vehicle')


LCC1_T <- read_csv('results/parameters/LCC1_treatment.csv') %>%
  mutate(condition = 'LCC1 (SENS), TRT') %>% 
  mutate(treatment = 'Treatment')

LCC9_V <- read_csv('results/parameters/LCC9_vehicle.csv') %>%
  mutate(condition = 'LCC9 (RES), VEH') %>% 
  mutate(treatment = 'Vehicle')

LCC9_T <- read_csv('results/parameters/LCC9_treatment.csv') %>%
  mutate(condition = 'LCC9 (RES), TRT') %>% 
  mutate(treatment = 'Treatment')

df.2 <- full_join(LCC1_V, LCC1_T) %>% 
  full_join(LCC9_V) %>% 
  full_join(LCC9_T) %>% 
  mutate(i = if_else(i == 1, 'Intrinsic growth rate (r)',
                 'Carrying capacity (K)'))

#################
# Panel 1

pal <- c(brewer.pal(9, 'Blues')[c(4, 7)], brewer.pal(9, 'Greens')[c(4, 7)])

df$condition <- ordered(df$condition, c('LCC1 (SENS), VEH','LCC1 (SENS), TRT', 'LCC9 (RES), VEH', 'LCC9 (RES), TRT'))

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
        strip.text.x = element_text(size = 18),
        plot.tag = element_text(size = 25))+
  facet_wrap(~condition)

#################
# Panel 2

df.2$i <- ordered(df.2$i, c('Intrinsic growth rate (r)', 'Carrying capacity (K)'))
df.2$condition <- ordered(df.2$condition, c('LCC1 (SENS), VEH', 'LCC1 (SENS), TRT', 'LCC9 (RES), VEH', 'LCC9 (RES), TRT'))

p2 <- df.2 %>% 
  filter(i == 'Intrinsic growth rate (r)') %>% 
  ggplot(aes(x = condition, y = theta, fill = condition))+
  stat_halfeye(orientation = 'vertical',
               show.legend = F,
               normalize = 'xy',
               .width = c(.95, .5),
               aes(fill_ramp = stat(cut_cdf_qi(cdf, 
                                               .width = c(.5, .95, 1), 
                                               labels = scales::percent_format(accuracy = 1)))))+
  scale_fill_ramp_discrete(range = c(.9, .3), na.translate = F)+
  #facet_wrap(~i, scales = 'free_y')+
  scale_x_discrete(labels = c('LCC1,\nVEH', 'LCC1,\nTRT', 'LCC9,\nVEH','LCC9,\nTRT'))+
  theme_minimal()+
  labs(title = 'Intrinsic growth rate (r)', x = '', y = '')+
  scale_fill_manual(values = pal[c(1, 3, 2, 4)])+
  labs(x = '', y = '', tag = 'B')+
  theme(strip.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 13),
        plot.tag = element_text(size = 25),
        plot.title = element_text(hjust = 0.5, size = 20))+
  ylim(0.01, 0.05)


p3 <- df.2 %>% 
  filter(i == 'Carrying capacity (K)') %>% 
  ggplot(aes(x = condition, y = theta, fill = condition))+
  stat_halfeye(orientation = 'vertical',
               show.legend = F,
               normalize = 'xy',
               .width = c(.95, .5),
               aes(fill_ramp = stat(cut_cdf_qi(cdf, 
                                               .width = c(.5, .95, 1), 
                                               labels = scales::percent_format(accuracy = 1)))))+
  scale_fill_ramp_discrete(range = c(.9, .3), na.translate = F)+
  #facet_wrap(~i, scales = 'free_y')+
  scale_x_discrete(labels = c('LCC1,\nVEH', 'LCC1,\nTRT', 'LCC9,\nVEH','LCC9,\nTRT'))+
  labs(title = 'Carrying capacity (K)', x = '', y = '')+
  theme_minimal()+
  scale_fill_manual(values = pal[c(1, 3, 2, 4)])+
  theme(strip.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 13),
        plot.tag = element_text(size = 25),
        plot.title = element_text(hjust = 0.5, size = 20))+
  ylim(3500, 16000)

#################
# Save figure

design <- "
  1123
  1123
"

ggsave('results/images/figures/fig2.jpeg', p1 + p2 + p3 + plot_layout(design = design), scale = 1, dpi = 300, width = 15, height = 8)


