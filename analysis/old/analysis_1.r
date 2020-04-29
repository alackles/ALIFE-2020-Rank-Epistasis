rm(list = ls())

library(ggplot2)

setwd('~/documents/school/cur_sp20/multi_disciplinary_research/project/MABE')

data = read.csv('org_data.csv')
num_traits = ncol(data) - 3
pop_size = nrow(data[data$gen == 0,])
  
data_per_gen = data.frame(data = matrix(nrow = 0, ncol = 3))


for(cur_gen in 1:max(unique(data$gen))){
  gen_data = data[data$gen == cur_gen,]
  for(int_val in 1:(2^(num_traits) - 1)){
    data_per_gen[nrow(data_per_gen) + 1, ] = c(cur_gen, int_val, nrow(gen_data[gen_data$int_val == int_val,]))
  }
  print(cur_gen)
}
colnames(data_per_gen) = c('gen', 'int_val', 'count')

data_per_gen$pct = data_per_gen$count / pop_size


ggplot(data_per_gen, aes(x = gen, y = int_val, fill = pct)) + 
  geom_tile()

ggplot(data_per_gen, aes(x = gen)) + 
  geom_line(data = data_per_gen[data_per_gen$int_val == 365,], aes(y = count, color = '101101101')) +
  geom_line(data = data_per_gen[data_per_gen$int_val == 146,], aes(y = count, color = '010010010'))

ggplot(data_per_gen[data_per_gen$count > 2,], aes(x = int_val)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = 365)) + 
  geom_vline(aes(xintercept = 146))

data_50 = data_per_gen[data_per_gen$gen == 50,]
data_50[order(data_50$count, decreasing = T),]
