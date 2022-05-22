###############
# Libraries
library(tidyverse)
source('cycling_functions.R') # loads data and functions
library(RColorBrewer)
library(patchwork)

###############
# Simulate assuming that Resistant is 50% less fit than sensitive
# in the presence of vehicle. All else is unchanged.

source('cycling_functions.R') # loads data and functions

#pars_treat$r_s.mean <- pars_veh$r_s.mean * (pars_treat$r_s.mean / pars_veh$r_s.mean)
pars_veh$r_r.mean <- pars_veh$r_s.mean * .5

#pars_treat$K_s.mean <- pars_veh$K_s.mean * (pars_treat$K_s.mean / pars_veh$K_s.mean)
pars_veh$K_r.mean <- pars_veh$K_s.mean * .5

###############
# Simulate

fit <- treat_cycle.rep(200, 200, 100, 5, c(1000, 1000), pars_treat, pars_veh)

###############
# Plot

pal <- c(brewer.pal(9, 'Blues')[c(4, 8)], brewer.pal(9, 'Greens')[c(4, 8)])
fit$condition <- ordered(fit$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))

p1 <- ggplot(fit, aes(x = t, y = N, group = interaction(rep, pop), color = condition))+
  geom_line(alpha = .1, show.legend = T)+
  scale_color_manual(values = pal[c(1, 3, 2, 4)])+
  theme_minimal()+
  labs(y = 'Cell count', 
       x = 'Time steps',
       color = 'Condition')+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.position = c(0.75, 0.15),
        legend.text=element_text(size=18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.title=element_text(size = 20),
        strip.text.x = element_text(size = 15),
        plot.tag = element_text(size = 25))+
  ylim(0, 16500)

###############
# Simulate assuming that Resistant is 33% less fit than sensitive
# in the presence of vehicle. All else is unchanged.

source('cycling_functions.R') # loads data and functions

#pars_treat$r_s.mean <- pars_veh$r_s.mean * (pars_treat$r_s.mean / pars_veh$r_s.mean)
pars_veh$r_r.mean <- pars_veh$r_s.mean * .66

#pars_treat$K_s.mean <- pars_veh$K_s.mean * (pars_treat$K_s.mean / pars_veh$K_s.mean)
pars_veh$K_r.mean <- pars_veh$K_s.mean * .66

###############
# Simulate

fit <- treat_cycle.rep(200, 200, 100, 5, c(1000, 1000), pars_treat, pars_veh)

###############
# Plot

pal <- c(brewer.pal(9, 'Blues')[c(4, 8)], brewer.pal(9, 'Greens')[c(4, 8)])
fit$condition <- ordered(fit$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))

p2 <- ggplot(fit, aes(x = t, y = N, group = interaction(rep, pop), color = condition))+
  geom_line(alpha = .1, show.legend = T)+
  scale_color_manual(values = pal[c(1, 3, 2, 4)])+
  theme_minimal()+
  labs(y = 'Cell count', 
       x = 'Time steps',
       color = 'Condition')+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.position = c(0.75, 0.15),
        legend.text=element_text(size=18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.title=element_text(size = 20),
        strip.text.x = element_text(size = 15),
        plot.tag = element_text(size = 25))+
  ylim(0, 16500)

###############
# Save fig

###############
# Simulate assuming that Resistant is 33% less fit than sensitive
# in the presence of vehicle. All else is unchanged.

source('cycling_functions.R') # loads data and functions

#pars_treat$r_s.mean <- pars_veh$r_s.mean * (pars_treat$r_s.mean / pars_veh$r_s.mean)
pars_veh$r_r.mean <- pars_veh$r_s.mean

#pars_treat$K_s.mean <- pars_veh$K_s.mean * (pars_treat$K_s.mean / pars_veh$K_s.mean)
pars_veh$K_r.mean <- pars_veh$K_s.mean

###############
# Simulate

fit <- treat_cycle.rep(200, 200, 100, 5, c(1000, 1000), pars_treat, pars_veh)

###############
# Plot

pal <- c(brewer.pal(9, 'Blues')[c(4, 8)], brewer.pal(9, 'Greens')[c(4, 8)])
fit$condition <- ordered(fit$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))

p3 <- ggplot(fit, aes(x = t, y = N, group = interaction(rep, pop), color = condition))+
  geom_line(alpha = .1, show.legend = T)+
  scale_color_manual(values = pal[c(1, 3, 2, 4)])+
  theme_minimal()+
  labs(y = 'Cell count', 
       x = 'Time steps',
       color = 'Condition')+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.position = c(0.75, 0.15),
        legend.text=element_text(size=18),
        axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        legend.title=element_text(size = 20),
        strip.text.x = element_text(size = 15),
        plot.tag = element_text(size = 25))+
  ylim(0, 16500)

###############
# Save fig

ggsave('results/images/supplementary figures/figS4.jpeg', p1 + p2 + p3, scale = 1, width = 15, height =8)





