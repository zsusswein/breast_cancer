
###############
# Libraries
library(tidyverse)
source('cycling_functions.R') # loads data and functions
library(RColorBrewer)

###############

fit <- treat_cycle.rep(100, 20, 20, 20, c(1000, 1000), pars_treat, pars_veh)

###############

pal <- c(brewer.pal(9, 'Blues')[c(4, 8)], brewer.pal(9, 'Greens')[c(4, 8)])
fit$condition <- ordered(fit$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))

p1 <- ggplot(fit, aes(x = t, y = N, group = interaction(rep, pop), color = condition))+
  geom_line(alpha = .1)+
  scale_color_manual(values = pal[c(1, 2, 3, 4)])+
  theme_minimal()+
  labs(y = 'Cell count', 
       x = 'Time steps',
       color = 'Condition',
       tag = 'A')+
  guides(colour = guide_legend(override.aes = list(alpha = 1)))+
  theme(legend.position = c(0.75, 0.15),
        legend.text=element_text(size=18),
        axis.text=element_text(size=10),
        axis.title=element_text(size=20),
        legend.title=element_text(size = 20),
        strip.text.x = element_text(size = 15))+
  ylim(0, 15000)

###############
# Refit assuming strong constant competition
pars_treat$beta.mean = 1
pars_treat$beta.sd = .0001
pars_treat$alpha.mean = 1
pars_treat$alpha.sd = .0001

pars_veh$beta.mean = 1
pars_veh$beta.sd = .0001
pars_veh$alpha.mean = 1
pars_veh$alpha.sd = .0001

fit <- treat_cycle.rep(100, 20, 20, 20, c(1000, 1000), pars_treat, pars_veh)

###########
pal <- c(brewer.pal(9, 'Blues')[c(4, 8)], brewer.pal(9, 'Greens')[c(4, 8)])
fit$condition <- ordered(fit$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))

p2 <- ggplot(fit, aes(x = t, y = N, group = interaction(rep, pop), color = condition))+
  geom_line(alpha = .1, show.legend = F)+
  scale_color_manual(values = pal[c(1, 2, 3, 4)])+
  theme_minimal()+
  labs(y = '', 
       x = 'Time steps',
       color = 'Condition',
       tag = 'B')+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 20))+
  ylim(0, 15000)

###########

ggsave('results/images/figures/fig5.jpeg', p1 | p2, scale = 2)

  
