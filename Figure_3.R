
#################
# Libraries
library(tidyverse)
library(tidybayes)
library(RColorBrewer)
library(patchwork)
library(ggdist)

#################
# Read in data

LV_V <- read_csv("results/fitted_draws/LV_vehicle.csv") %>% 
  mutate(treatment = 'Vehicle')

LV_T <- read_csv("results/fitted_draws/LV_treatment.csv") %>% 
  mutate(treatment = 'Treatment')

df <- full_join(LV_V, LV_T) %>% 
  mutate(condition = if_else(cell_type == 'LCC1', 
                             if_else(treatment == 'Vehicle', 'LCC1, Vehicle', 'LCC1, Treatment'),
                             if_else(treatment == 'Vehicle', 'LCC9, Vehicle', 'LCC9, Treatment'))) %>% 
  filter(.width == .95) %>% 
  mutate(condition = as.factor(condition)) %>% 
  mutate(condition = relevel(condition, ref = 'LCC1, Vehicle'))


###

LV_V <- read_csv("results/parameters/LV_vehicle.csv") %>% 
  mutate(condition = if_else((i == 3) | (i == 5), 'LCC1, Vehicle', 'LCC9, Vehicle'))

LV_T <- read_csv("results/parameters/LV_treatment.csv") %>% 
  mutate(condition = if_else((i == 3) | (i == 5), 'LCC1, Treatment', 'LCC9, Treatment'))

df.2 <- full_join(LV_V, LV_T) %>% 
  mutate(condition = as.factor(condition)) %>% 
  mutate(condition = relevel(condition, ref = 'LCC1, Vehicle')) %>% 
  filter(i > 2) %>% 
  mutate(i = if_else((i == 3) | (i == 4), 'Intrinsic growth rate (r)',
                     'Carrying capacity (K)'))


#################

pal <- c(brewer.pal(9, 'Blues')[c(4, 7)], brewer.pal(9, 'Greens')[c(4, 7)])

df$treatment <- ordered(df$treatment, c('Vehicle','Treatment'))

df$condition <- ordered(df$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))


p1 <- ggplot(df, aes(x = rank, group = condition))+
  geom_line(aes(y = y_rep_mean, color = condition), alpha= 0.75, size = 1.25)+
  geom_point(aes(y = agg_count, color = condition), size = 1.75)+
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = interaction(condition, .width), fill = condition), alpha = .35, show.legend = F)+
  theme_minimal()+
  scale_color_manual(values = pal[c(1, 3, 2, 4)])+
  scale_fill_manual(values = pal[c(1, 3, 2, 4)])+
  labs(y = 'Cell count',
       x = 'Time steps',
       color = 'Condition',
       tag = 'A')+
  theme(legend.position = c(0.2, 0.9),
        legend.text=element_text(size=18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.title=element_text(size = 20),
        strip.text.x = element_text(size = 20))+
  facet_wrap(~treatment)

#################

df.2$i <- ordered(df.2$i, c('Intrinsic growth rate (r)', 'Carrying capacity (K)'))

p2 <- ggplot(df.2, aes(x = theta, y = condition, fill = condition))+
  stat_halfeye(orientation = 'horizontal',
               show.legend = F,
               normalize = 'xy',
               .width = c(.95, .5),
               aes(fill_ramp = stat(cut_cdf_qi(cdf, 
                                               .width = c(.5, .95, 1), 
                                               labels = scales::percent_format(accuracy = 1)))))+
  scale_fill_ramp_discrete(range = c(.9, .3), na.translate = F)+
  facet_wrap(~i, scales = 'free_x')+
  scale_y_discrete(limits = c('LCC9, Treatment','LCC9, Vehicle', 'LCC1, Treatment', 'LCC1, Vehicle'))+
  theme_minimal()+
  scale_fill_manual(values = pal[c(1, 2, 4, 3)])+
  labs(x = '', y = '', tag = 'B')+
  theme(strip.text.x = element_text(size = 20),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 20))



#################
# Save figure

ggsave('results/images/figures/fig3.jpeg', p1 | p2, scale = 2)
