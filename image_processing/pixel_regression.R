library(tidyverse)

data <- read_csv('monoculture_estimate.csv')

ggplot(data, aes(count, estimate, label=X1))+
  geom_point()+
  geom_label(aes(color = X1))

# the log-log has slight better properties but he estimate is about the same
fit <- lm(data=data, count~estimate)
summary(fit)
  
ggplot(data, aes(estimate, count, label=X1))+
  geom_point()+
  geom_line(aes(estimate, fit$fitted.values))+
  geom_abline(slope=1, alpha=.2, linetype=2)+
  labs(x='Algorithm count',
       y='Manual count',
       tag = 'B')+
  theme_minimal()+
  theme(axis.text=element_text(size=10),
  axis.title=element_text(size=15))
ggsave('monoculture_corr_plot.jpg')

cor.test(data$estimate, data$count, alternative='two.sided')

fit <- lm(data=data, scale(count)~scale(estimate))
summary(fit)

ggplot(data)+
  geom_histogram(aes(difference/count), fill='black', bins=20)+
  theme_classic()+
  labs(title = 'Histogram of absolute deviation', x='Absolute deviation (algorithm count to manual)', y='Count')
8