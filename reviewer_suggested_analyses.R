# This file implements simulations of treatment cycling with different treatment 
# frequencies, ranging from no treatment to continuous treatment with intermediate
# frequencies in between.

# Libraries ---------------------------------------------------------------
library(tidyverse)
source('cycling_functions.R') # loads data and functions
library(RColorBrewer)
library(patchwork)


# Handler function --------------------------------------------------------
# Set up function to run treatment sim and generate the associated figure

fit_analysis <- function(nrep, t_treat, t_veh, cycles, y0, pars_treat, pars_veh){
  
  fit <- treat_cycle.rep(nrep, t_treat, t_veh, cycles, y0, pars_treat, pars_veh)
  
  pal <- c(brewer.pal(9, 'Blues')[c(4, 8)], brewer.pal(9, 'Greens')[c(4, 8)])
  fit$condition <- ordered(fit$condition, c('LCC1, Vehicle','LCC1, Treatment', 'LCC9, Vehicle', 'LCC9, Treatment'))
  
  p1 <- ggplot(fit, aes(x = t, y = N, group = interaction(rep, pop), color = condition))+
    geom_line(alpha = .1, show.legend = F)+
    scale_color_manual(values = pal[c(1, 3, 2, 4)])+
    theme_minimal()+
    labs(y = 'Cell count', 
         x = 'Time steps',
         color = 'Condition')+
    guides(colour = guide_legend(override.aes = list(alpha = 1)))+
    theme(legend.position = c(0.75, 0.15),
          legend.text=element_text(size=18),
          axis.text=element_text(size=20),
          axis.title=element_text(size=20),
          legend.title=element_text(size = 20),
          strip.text.x = element_text(size = 15),
          plot.tag = element_text(size = 25))+
    ylim(0, 16500)
  
  return(p1)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
# Set up grid -------------------------------------------------------------

times <- expand_grid(t_treat = seq(0, 150, 50),
                     t_veh = seq(0, 150, 50)) %>% 
  filter(!(t_treat == 0 & t_veh == 0))


myplots <- list(ggplot())

for(i in 1:nrow(times)){
  
  fig <- fit_analysis(100, times[[i,1]], times[[i,2]], 10, c(1000, 1000), pars_treat, pars_veh)
    
  print(i)
  print(fig)
  myplots[[i+1]] <- fig  # add each plot into plot list}
  
}

library(gridExtra)
n <- length(myplots)
nCol <- floor(sqrt(n))
do.call("grid.arrange", c(myplots, ncol=4))

ggsave('results/images/cycling_grid.png')
