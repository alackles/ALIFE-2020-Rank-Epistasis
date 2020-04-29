rm(list = ls())
library(dplyr)

#data = read.csv('./data/raw/edit_distance_scraped.csv')
data = read.csv('./data/raw/edit_distance_TRD_1_final.csv')
#data_fit = read.csv('./data/raw/edit_distance_scraped_fit.csv')
#data = rbind(data, data_fit)
data = data[data$seed != 501,]
colnames(data) = c('X', 'X2', 'update', 'bit_idx', 'edit_distance',
                   'unique_genomes', 'condition' ,'seed', 'FIX', 'K', 'ED', 'MUT', 'TRD', 'VEL')
data = data[data$ED == 1 & data$update <= 2500,]

data_grouped = dplyr::group_by(data, update, seed, FIX, K, ED, MUT, TRD, VEL)
data_summary = dplyr::summarize(data_grouped, 
                                edit_distance_mean = mean(edit_distance), 
                                edit_distance_max = max(edit_distance), 
                                edit_distance_min = min(edit_distance),
                                average_unique_genomes = mean(unique_genomes))

data_mutant = read.csv('data/data_mutant_fitness__seed_all_TRD_0.csv')
data_mutant = data_mutant[data_mutant$ED ==1, ]
data_mutant_grouped = dplyr::group_by(data_mutant, update, FIX,K, MUT, ED, seed)
data_mutant_summary = dplyr::summarize(data_mutant_grouped, 
                                alpha_mean = mean(alpha_el),
                                beta_mean = mean(beta_el))

data_summary$alpha = NA
data_summary$beta = NA
for(update in unique(data_summary$update)){
  print(update)
  data_update = data_summary[data_summary$update == update,]
  data_mutant_update = data_mutant_summary[data_mutant_summary$update == update,]
  for(fix in unique(data$FIX)){
    data_fix = data_update[data_update$FIX == fix,]
    data_mutant_fix = data_mutant_update[data_mutant_update$FIX == fix,]
    for(k in unique(data_fix$K)){
      data_k = data_fix[data_fix$K == k,]
      data_mutant_k = data_mutant_fix[data_mutant_fix$K == k,]
      for(mut in unique(data_k$MUT)){
        data_mut = data_k[data_k$MUT == mut,]
        data_mutant_mut = data_mutant_k[data_mutant_k$MUT == mut,]
        for(ed in unique(data_mut$ED)){
          data_ed = data_mut[data_mut$ED == ed,]
          data_mutant_ed = data_mutant_mut[data_mutant_mut$ED == ed,]
          for(vel in unique(data_ed$VEL)){
            data_vel = data_ed[data_ed$VEL == vel,]
            data_mutant_vel = data_mutant_ed[data_mutant_ed$ED == vel,]
            for(seed in unique(data_ed$seed)){
              data_mutant_seed = data_mutant_ed[data_mutant_ed$seed == seed,]
              if(nrow(data_mutant_seed) == 0){
                next
              }
              data_summary[
                data_summary$update == update &
                data_summary$FIX == fix &
                data_summary$K == k &
                data_summary$MUT == mut &
                data_summary$seed == seed,
                  ]$alpha = data_mutant_seed$alpha_mean
              data_summary[
                data_summary$update == update &
                data_summary$FIX == fix &
                data_summary$K == k &
                data_summary$MUT == mut &
                data_summary$seed == seed,
                  ]$beta = data_mutant_seed$beta_mean
            }
          }
        }
      }
    }
  }
}
data_summary$K_str = paste0('K = ', data_summary$K)
data_summary$K_factor = factor(data_summary$K_str, levels  = paste0('K = ', c(3,5,10)))
data_summary$FIX_str = 'Fixed NK'
data_summary[data_summary$FIX == 0,]$FIX_str = 'Random NK'
data_summary$ED_str = 'Levenshtein'
data_summary[data_summary$ED == 1,]$ED_str = 'Damerau-Levenshtein'
data_summary$MUT_str = paste0('Mut. Rate: ', data_summary$MUT)
data_summary$MUT_factor = factor(data_summary$MUT_str, levels  = paste0('Mut. Rate: ', c(0.001, 0.01, 0.1)))
write.csv(data_summary, 'data/data_combined.csv')

