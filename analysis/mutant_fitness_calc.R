rm(list = ls())
library(ggplot2)

# Load in data, make column names more legible
data = read.csv('./data/raw/mutant_fitness_TRD_0_final.csv')
colnames(data) = c('X', 'X2', 'update', 'org_idx', 'num_mutations', 'fitness_avg', 'fitness_max', 'fitness_min', 
                   'condition', 'seed', 'FIX', 'K', 'ED', 'MUT', 'TRD', 'VEL')
# Give each row an identifier
data$row_id = 1:nrow(data)


print(paste0('All rows: ', nrow(data)))
# Calculate additioinal columns for original organisms
for(seed in c(502:550)){
  # Identify original organisms (no mutations)
  data_originals = data[data$num_mutations == 0 & data$seed == seed,]
  # Setup additional columns that will be filled out
  data_originals$alpha = NA
  data_originals$beta = NA
  data_originals$one_step_fitness_avg = NA
  data_originals$two_step_fitness_avg = NA
  data_originals$alpha_el = NA
  data_originals$beta_el = NA
  data_seed = data[data$seed == seed & data$ED == 0,]
  print(seed)
  print(nrow(data_seed))
  for(fix in unique(data_seed$FIX)){
    print(paste0('  Fix: ', fix))
    data_fix = data_seed[data_seed$FIX == fix,]
    print(paste0('  ', nrow(data_fix)))
    for(mut in unique(data_fix$MUT)){
      print(paste0('    Mut: ', mut))
      data_mut = data_fix[data_fix$MUT == mut,]
      print(paste0('    ', nrow(data_mut)))
      for(update in unique(data_mut$update)){
        print(paste0('      Update: ', update))
        data_update = data_mut[data_mut$update == update,]
        print(paste0('      ', nrow(data_update)))
        for(ed in unique(data_update$ED)){
          print(paste0('  ED: ', ed))
          data_ed = data_update[data_update$ED == ed,]
          print(paste0('  ED: ', nrow(data_ed)))
          for(k in unique(data_ed$K)){
            data_k = data_ed[data_ed$K == k,]
            print(paste0('        K: ', k))
            print(paste0('        ', nrow(data_k)))
            for(org_idx in unique(data_k$org_idx)){
              data_org = data_k[data_k$org_idx == org_idx,]
              row = data_org[data_org$num_mutations == 0,]
              row_one = data_org[data_org$num_mutations == 1,]
              row_two = data_org[data_org$num_mutations == 2,]
              #alpha = -1 * log(row_one$fitness_avg/row$fitness_avg)
              #beta = log2((-1 * log(row_two$fitness_avg/row$fitness_avg)) / alpha)
              x_1 = 1
              y_1 = log(row_one$fitness_avg / row$fitness_avg)
              x_2 = 2
              y_2 = log(row_two$fitness_avg / row$fitness_avg)
              beta_el = (y_2 - (x_2 / x_1) * y_1) / (x_2 * x_1 - x_2^2)
              alpha_el = (-1 * beta_el * x_1^2 - y_1) / x_1
              #data_originals[data_originals$row_id == row$row_id,]$alpha = alpha
              #data_originals[data_originals$row_id == row$row_id,]$beta = beta
              data_originals[data_originals$row_id == row$row_id,]$one_step_fitness_avg = row_one$fitness_avg
              data_originals[data_originals$row_id == row$row_id,]$two_step_fitness_avg = row_two$fitness_avg
              data_originals[data_originals$row_id == row$row_id,]$alpha_el = alpha_el
              data_originals[data_originals$row_id == row$row_id,]$beta_el = beta_el
            }
          }
        }
      }
    }
  }
  data_write = data_originals[!is.na(data_originals$one_step_fitness_avg),]
  write.csv(data_write[data_write$seed == seed,], paste0('./data/mutant_fitness_parts/data_mutant_fitness__seed_', seed, '.csv'))
}

data_all = data.frame(data = matrix(nrow = 0, ncol = ncol(data_write)))
for(seed in 502:550){
  print(seed)
  data_seed = read.csv(paste0('./data/mutant_fitness_parts/data_mutant_fitness__seed_', seed, '.csv'))
  data_all = rbind(data_all, data_seed)
}
write.csv(data_all[3:ncol(data_all)], paste0('./data/data_mutant_fitness__seed_all.csv'))

