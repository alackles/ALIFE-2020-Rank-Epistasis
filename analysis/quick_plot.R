rm(list = ls())

library(ggplot2)
setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/cse845/')

data = read.csv('mutant_data.csv')

updates_to_plot = c(0, 100, 250, 500, 1000, 1500, 2000)
data_plot = data[data$update %in% updates_to_plot,]

ggplot(data_plot, aes(x = bit_idx, y = edit_distance)) + 
  geom_line() + 
  facet_grid(rows = vars(update))