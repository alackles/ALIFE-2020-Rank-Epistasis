rm(list = ls())
library(ggplot2)
setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/MABE')

data_fit = read.csv('snapshots/3/fit/pop.csv')
data_fit$treatment = 'Fit'
data_flat = read.csv('snapshots/3/flat/pop.csv')
data_flat$treatment = 'Flat'

data_fit_flat = read.csv('snapshots/3/fit_flat/pop.csv')
data_fit_flat$treatment = 'Fit'
data_fit_flat$update = data_fit_flat$update + 501
data_flat_fit = read.csv('snapshots/3/flat_fit/pop.csv')
data_flat_fit$treatment = 'Flat'
data_flat_fit$update = data_flat_fit$update + 501

data = rbind(data_fit, data_flat, data_fit_flat, data_flat_fit)

ggplot(data[data$update <= 500,], aes(update, score_AVE, color = treatment)) + 
  geom_hline(aes(yintercept = 1, color = 'Flat'), linetype = 2, alpha = 0.75) + 
  geom_hline(aes(yintercept = 2, color = 'Fit'), linetype = 2, alpha = 0.75) + 
  geom_vline(aes(xintercept = 500), linetype = 2, alpha = 0.75)+
  geom_line() +
  xlab('Update') + 
  ylab('Average score') +
  labs(color = 'Environment') + 
  theme(legend.position = 'bottom') + 
  scale_x_continuous(limits = c(0, 1501)) +
  ggsave('transition_prior.pdf', units = 'in', width = 8, height = 6) +
  ggsave('transition_prior.png', units = 'in', width = 8, height = 6)

ggplot(data, aes(update, score_AVE, color = treatment)) + 
  geom_hline(aes(yintercept = 1, color = 'Flat'), linetype = 2, alpha = 0.75) + 
  geom_hline(aes(yintercept = 2, color = 'Fit'), linetype = 2, alpha = 0.75) + 
  geom_vline(aes(xintercept = 500), linetype = 2, alpha = 0.75)+
  geom_line() +
  xlab('Update') + 
  ylab('Average score') +
  labs(color = 'Environment') + 
  theme(legend.position = 'bottom') + 
  scale_x_continuous(limits = c(0, 1501)) +
  ggsave('transition.pdf', units = 'in', width = 8, height = 6) +
  ggsave('transition.png', units = 'in', width = 8, height = 6)
