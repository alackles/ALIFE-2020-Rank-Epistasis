rm(list = ls())
library(ggplot2)

setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/analysis/data/2')
data_fit_L = read.csv('edit_distance_fit_L.csv')
data_fit_L$locus = as.numeric(data_fit_L$locus)
data_fit_L$gen = as.numeric(data_fit_L$gen)
data_fit_L$edit_distance = as.numeric(data_fit_L$edit_distance)
data_fit_L$environment = 'Fit'
data_fit_L$metric = 'Levenshtein'

data_fit_DL = read.csv('edit_distance_fit_DL.csv')
data_fit_DL$locus = as.numeric(data_fit_DL$locus)
data_fit_DL$gen = as.numeric(data_fit_DL$gen)
data_fit_DL$edit_distance = as.numeric(data_fit_DL$edit_distance)
data_fit_DL$environment = 'Fit'
data_fit_DL$metric = 'Damerau-Levenshtein'

data_flat_L = read.csv('edit_distance_flat_L.csv')
data_flat_L$locus = as.numeric(data_flat_L$locus)
data_flat_L$gen = as.numeric(data_flat_L$gen)
data_flat_L$edit_distance = as.numeric(data_flat_L$edit_distance)
data_flat_L$environment = 'Flat'
data_flat_L$metric = 'Levenshtein'

data_flat_DL = read.csv('edit_distance_flat_DL.csv')
data_flat_DL$locus = as.numeric(data_flat_DL$locus)
data_flat_DL$gen = as.numeric(data_flat_DL$gen)
data_flat_DL$edit_distance = as.numeric(data_flat_DL$edit_distance)
data_flat_DL$environment = 'Flat'
data_flat_DL$metric = 'Damerau-Levenshtein'

data = rbind(data_fit_L, data_flat_L, data_fit_DL, data_flat_DL)
data$gen_str = paste0('Generation ', data$gen)

ggplot(data[data$gen %in% c(0,250, 500),], aes(x = locus, y = edit_distance, color = metric)) + 
  geom_line() + 
  facet_grid(rows = vars(gen_str), cols = vars(environment))

gens = c(500)
ggplot(data[data$gen %in% gens & data$environment == 'Fit' & data$locus <= 50,], aes(x = locus, y = edit_distance, fill = metric)) + 
  geom_col(data[data$gen %in% gens & data$environment == 'Fit' & data$metric == 'Levenshtein' & data$locus <= 50,], mapping = aes(), alpha = 0.5) + 
  geom_col(data[data$gen %in% gens & data$environment == 'Fit' & data$metric == 'Damerau-Levenshtein' & data$locus <= 50,], mapping = aes(), alpha = 0.5) + 
  facet_grid(rows = vars(gen_str), cols = vars(environment)) + 
  theme(legend.position = 'bottom')

# ggplot(data[data$gen %in% c(0,250,499),], aes(x = edit_distance)) +
#   geom_histogram(binwidth = 5) +
#   facet_grid(rows = vars(treatment), cols = vars(gen_str)) +
#   xlab('Edit Distance')

y = data[data$environment == 'Fit' & data$metric == 'Levenshtein',]$edit_distance - data[data$environment == 'Fit' & data$metric == 'Damerau-Levenshtein',]$edit_distance 
x = 1:length(y)
y = sort(y)
data_diff = data.frame(x = x, y = y)
ggplot(data_diff, aes(x, y)) + 
  geom_line() +
  xlab('x') + 
  ylab('Levenshtein Distance - Damerau-Levenshtein Distance') + 
  ggsave('distance_metric_difference_line.pdf', units = 'in', width = 8, height = 6) +
  ggsave('distance_metric_difference_line.png', units = 'in', width = 8, height = 6)

ggplot(data_diff, aes(x = y)) + 
  geom_histogram(binwidth = 1, fill = '#66ee99', color = '#000000') +
  xlab('Levenshtein Distance - Damerau-Levenshtein Distance') + 
  ylab('Count') + 
  ggsave('distance_metric_difference_histogram.pdf', units = 'in', width = 8, height = 6) +
  ggsave('distance_metric_difference_histofram.png', units = 'in', width = 8, height = 6)

ggplot(data[data$gen %in% c(0,250, 500),], aes(x = locus, color = metric)) + 
  geom_line(mapping = aes(y = edit_distance, color = 'Edit Distance')) + 
  geom_line(mapping = aes(y = weighted_edit_distance, color = 'Weighted Edit Distance')) + 
  facet_grid(cols = vars(gen_str), rows = vars(environment)) +
  theme(legend.position = 'bottom') +
  ylab('Edit Distance') +
  xlab('Locus') +
  labs(color = 'Value') + 
  ggsave('weight_effect.pdf', units = 'in', width = 8, height = 6) + 
  ggsave('weight_effect.png', units = 'in', width = 8, height = 6) 

pop_size = 200
data_plot = data[data$locus == 1 & data$edit_distance != 0,]
data_plot$unique_genotypes = (data_plot$weighted_edit_distance * pop_size) / data_plot$edit_distance

ggplot(data_plot, aes(x = gen, color = metric)) + 
  geom_line(mapping = aes(y = unique_genotypes, color = 'Edit Distance')) + 
  facet_grid(cols = vars(environment)) +
  theme(legend.position = 'bottom') +
  ylab('Unique Genotypes') +
  xlab('Generation') +
  labs(color = 'Value') + 
  ggsave('unique_genomes.pdf', units = 'in', width = 8, height = 6) +
  ggsave('unique_genomes.png', units = 'in', width = 8, height = 6)


ggplot(data_plot, aes(x = gen, color = metric)) + 
  geom_line(mapping = aes(y = weighted_edit_distance/edit_distance, color = 'Edit Distance')) + 
  facet_grid(cols = vars(environment)) +
  theme(legend.position = 'bottom') +
  ylab('Edit Distance') +
  xlab('Locus') +
  labs(color = 'Value') + 
  scale_y_continuous(limits = c(0,1))
