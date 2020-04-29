###################################
####            FIT            ####
###################################
rm(list = ls())
library(ggplot2)

setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/scripts/data/1')
data = read.csv('edit_distance_fit.csv')
data$locus = as.numeric(data$locus)
data$gen = as.numeric(data$gen)
data$edit_distance = as.numeric(data$edit_distance)

ggplot(data[data$gen < 10,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen))

ggplot(data[data$locus < 10,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus))

ggplot(data[data$gen > 495,], aes(x = edit_distance)) + 
  geom_histogram(binwidth = 1) + 
  facet_grid(rows = vars(gen))

ggplot(data[data$locus < 10,], aes(x = gen, y = edit_distance, color = locus)) + 
  geom_point()

ggplot(data, aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggsave('locus_faceted_by_gen_fit.pdf', units = 'in', width = 20, height = 500, limitsize = F)

ggplot(data, aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_fit.pdf', units = 'in', width = 20, height = 200, limitsize = F)

ggplot(data[data$locus == 1,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_1_fit.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$locus == 116,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_116_fit.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$locus == 155,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_155_fit.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$gen == 1,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggtitle('Generation 1') +
  scale_y_continuous(limits = c(0,20)) +
  ggsave('locus_faceted_by_gen_1_fit.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$gen == 250,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggtitle('Generation 250') +
  scale_y_continuous(limits = c(0,20)) +
  ggsave('locus_faceted_by_gen_250_fit.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$gen == 500,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggtitle('Generation 500') +
  scale_y_continuous(limits = c(0,20)) +
  ggsave('locus_faceted_by_gen_500_fit.png', units = 'in', width = 10, height = 2, limitsize = F)

###################################
####            FLAT           ####
###################################
rm(list = ls())
library(ggplot2)

setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/scripts/data/1')
data = read.csv('edit_distance_flat.csv')
data$locus = as.numeric(data$locus)
data$gen = as.numeric(data$gen)
data$edit_distance = as.numeric(data$edit_distance)

ggplot(data[data$gen < 10,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen))

ggplot(data[data$locus > 195,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus))

ggplot(data[data$gen > 490,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen))

ggplot(data[data$gen > 495,], aes(x = edit_distance)) + 
  geom_histogram(binwidth = 1) + 
  facet_grid(rows = vars(gen))

ggplot(data[data$locus < 10,], aes(x = gen, y = edit_distance, color = locus)) + 
  geom_point()

ggplot(data, aes(x = locus, y = edit_distance)) + 
  geom_point() + 
  facet_grid(rows = vars(gen)) + 
  ggsave('locus_faceted_by_gen_flat.pdf', units = 'in', width = 20, height = 500, limitsize = F)

ggplot(data[data$gen == 1,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggtitle('Generation 1') +
  scale_y_continuous(limits = c(0,20)) +
  ggsave('locus_faceted_by_gen_1_flat.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$gen == 250,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggtitle('Generation 250') +
  scale_y_continuous(limits = c(0,20)) +
  ggsave('locus_faceted_by_gen_250_flat.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$gen == 500,], aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggtitle('Generation 500') +
  scale_y_continuous(limits = c(0,20)) +
  ggsave('locus_faceted_by_gen_500_flat.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data, aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_flat.pdf', units = 'in', width = 20, height = 200, limitsize = F)


ggplot(data, aes(x = locus, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(gen)) + 
  ggsave('locus_faceted_by_gen_flat_2.pdf', units = 'in', width = 20, height = 500, limitsize = F)

ggplot(data[data$locus == 1,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_1_flat.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$locus == 116,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_116_flat.png', units = 'in', width = 10, height = 2, limitsize = F)

ggplot(data[data$locus == 155,], aes(x = gen, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(locus)) + 
  ggsave('gen_faceted_by_locus_155_flat.png', units = 'in', width = 10, height = 2, limitsize = F)


############################################
rm(list = ls())
library(ggplot2)

setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/scripts/data/1')
data_fit = read.csv('edit_distance_fit.csv')
data_fit$locus = as.numeric(data_fit$locus)
data_fit$gen = as.numeric(data_fit$gen)
data_fit$edit_distance = as.numeric(data_fit$edit_distance)

data_flat = read.csv('edit_distance_flat.csv')
data_flat$locus = as.numeric(data_flat$locus)
data_flat$gen = as.numeric(data_flat$gen)
data_flat$edit_distance = as.numeric(data_flat$edit_distance)


num_gens = length(unique(data_fit$gen))
data_integral = data.frame(data = matrix(nrow = num_gens * 2, ncol = 4))
for(gen in 1:num_gens){
  integral_fit = sum(data_fit[data_fit$gen == gen,]$edit_distance)
  integral_flat = sum(data_flat[data_flat$gen == gen,]$edit_distance)
  var_fit = var(data_fit[data_fit$gen == gen,]$edit_distance)
  var_flat = var(data_flat[data_flat$gen == gen,]$edit_distance)
  data_integral[gen,] = c(gen, integral_fit, 'fit', var_fit)
  data_integral[gen + num_gens,] = c(gen, integral_flat, 'flat', var_flat)
}

colnames(data_integral) = c('gen', 'integral', 'landscape', 'var')
data_integral$gen = as.numeric(data_integral$gen)
data_integral$integral = as.numeric(data_integral$integral)
data_integral$var = as.numeric(data_integral$var)

ggplot(data_integral, aes(gen, integral, color = landscape)) + 
  geom_line() + 
  ylab('Area under the curve (sum of edit distance / loci)') +
  ggsave('integral.png', units = 'in', width = 10, height = 4, limitsize = F)

ggplot(data_integral, aes(gen, var, color = landscape)) + 
  geom_line()
