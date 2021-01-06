## This file creates figure 3

#################
# Libraries
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(patchwork)
library(ggdist)

#################
# Read in data


LV_V <- read_csv("results/parameters/LV_vehicle.csv") %>% 
  mutate(treatment = 'Vehicle')

LV_T <- read_csv("results/parameters/LV_treatment.csv") %>% 
  mutate(treatment = 'Treatment')

df <- full_join(LV_V, LV_T) %>% 
  filter(i < 3) %>% 
  mutate(i = if_else(i == 1, 'Net effect of LCC9 \non LCC1 (beta)', 'Net effect of LCC1 \non LCC9 (alpha)')) %>% 
  mutate(treatment = as.factor(treatment)) %>% 
  mutate(treatment = relevel(treatment, ref = 'Vehicle')) 


LV_V <- read_csv("results/phase_space_draws/LV_vehicle.csv") %>% 
  mutate(treatment = 'Vehicle')

LV_T <- read_csv("results/phase_space_draws/LV_treatment.csv") %>% 
  mutate(treatment = 'Treatment')

df.2 <- full_join(LV_V, LV_T) %>% 
  rename(LCC1 = `1`,
         LCC9 = '2',
         rank = i)

#################
# Make figure

pal <- c(brewer.pal(9, 'Greens')[6], brewer.pal(9, 'Blues')[6])

ggplot(df, aes(x = theta, y = as.factor(i), fill = treatment, group = treatment))+
  stat_halfeye(orientation = 'horizontal', .width = c(.95, .5),
               aes(fill = treatment,
                   fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .95, 1), labels = scales::percent_format(accuracy = 1)))))+
  scale_fill_ramp_discrete(range = c(.9, .3), na.translate = F)+
  scale_fill_manual(values = pal[c(2, 1)])+
  theme_minimal()+
  labs(x = 'Per-capita interspecific decrease in available carrying capacity',
       y = '',
       fill = 'Condition',
       fill_ramp = 'Interval')+
  theme(legend.position = c(0.75, 0.65),
        legend.text=element_text(size=10),
        axis.text=element_text(size=15),
        axis.title=element_text(size=15),
        legend.title=element_text(size = 15))

#################
# Save figure

ggsave('results/images/figures/fig4.jpeg', scale = 1)  

#################

ggplot(df.2, aes(x = LCC1, y = LCC9, group = interaction(.draw, treatment), color = treatment))+
  geom_line(alpha = .025)+
  geom_abline(linetype = 2, color = 'grey')+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme_minimal()



