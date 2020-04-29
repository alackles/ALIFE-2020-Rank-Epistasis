rm(list = ls())
library(ggplot2)

data = read.csv('./data/data_combined_TRD_0.csv')
data = data[!is.na(data$beta),]

min_x = min(data$edit_distance_mean)
max_x = max(data$edit_distance_mean)
min_y = min(data$beta)
max_y = max(data$beta)
m = (max_y - min_y) / (max_x - min_x)
b = min_y - m * min_x

x = seq(min_x, max_x, length.out = 1000)
y = m * x + b
data_line = data.frame(x = x, y = y)

ggplot(data, aes(x = edit_distance_mean, y = beta)) + 
  geom_line(data_line, mapping = aes(x = x, y = y)) +
  geom_point(aes(color = update), alpha = 0.7) + 
  xlab('Mean edit distance across bits') + 
  ylab('Mean beta across orgs') + 
  ggsave('./plots/correlation.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/correlation.png', units = 'in', width = 8, height = 6)


data_plot = data[data$FIX == 0,]
ggplot(data_plot, aes(x = edit_distance_mean, y = beta)) + 
  geom_point(aes(color = update), alpha = 0.7) + 
  xlab('Mean edit distance across bits') + 
  ylab('Mean beta across orgs') + 
  facet_grid(cols = vars(MUT_factor), rows = vars(factor(K_factor, levels = paste0('K = ', c(10, 5, 3))))) + 
  ggsave('./plots/correlation_random.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/correlation_random.png', units = 'in', width = 8, height = 6)

data_plot = data[data$FIX == 1,]
ggplot(data_plot, aes(x = edit_distance_mean, y = beta)) + 
  geom_point(aes(color = update), alpha = 0.7) + 
  xlab('Mean edit distance across bits') + 
  ylab('Mean beta across orgs') + 
  facet_grid(cols = vars(MUT_factor), rows = vars(factor(K_factor, levels = paste0('K = ', c(10, 5, 3))))) + 
  ggsave('./plots/correlation_fixed.pdf', units = 'in', width = 8, height = 6) + 
  ggsave('./plots/correlation_fixed.png', units = 'in', width = 8, height = 6)

data_plot = data[data$FIX == 0 & data$update == 2500,]
ggplot(data_plot, aes(x = edit_distance_mean, y = beta)) + 
  geom_point(alpha = 0.7) + 
  xlab('Mean edit distance across bits') + 
  ylab('Mean beta across orgs') + 
  facet_grid(cols = vars(MUT_factor), rows = vars(factor(K_factor, levels = paste0('K = ', c(10, 5, 3))))) + 
  ggsave('./plots/correlation_random_2500.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/correlation_random_2500.png', units = 'in', width = 8, height = 6)

data_stats = data.frame(data = matrix(nrow = 0, ncol = 9))
colnames(data_stats) = c('fix', 'mut', 'k', 'correlation', 'p_value', 't', 'conf_int_level', 'conf_int_lower', 'conf_int_upper')
data_test = data[data$FIX == 0 & data$update == 2500,]
for(mut in unique(data_test$MUT)){
  data_test_mut = data_test[data_test$MUT == mut,]
  for(k in unique(data_test_mut$K)){
    data_test_k = data_test_mut[data_test_mut$K == k,]
    results = cor.test(data_test_k$edit_distance_mean, data_test_k$beta)
    data_stats[nrow(data_stats) + 1,] = c(
      0, 
      mut, 
      k,
      results[['estimate']],
      results[['p.value']],
      results[['statistic']],
      0.95,
      results[['conf.int']][1],
      results[['conf.int']][2]
      )
  }  
}
write.csv(data_stats, './stats/random_2500_correlation.csv')

