######## Old data loading for edit_distance_plots.R #########

data = read.csv('edit_distance_scraped.csv')
data_fit = read.csv('edit_distance_scraped_fit.csv')
data = rbind(data, data_fit)
data = data[data$seed != 101,]
colnames(data) = c('X', 'update', 'bit_idx', 'edit_distance', 'seed', 'FIX', 'K', 'ED', 'MUT')
data = data[data$ED == 0,]
data_grouped = dplyr::group_by(data, update, FIX, ED, K, MUT, seed)
data_summary = dplyr::summarise(data_grouped, 
  edit_distance_mean = mean(edit_distance), 
  edit_distance_max = max(edit_distance),
  edit_distance_min = min(edit_distance))

data_summary$K_str = paste0('K = ', data_summary$K)
data_summary$K_factor = factor(data_summary$K_str, levels  = paste0('K = ', c(3,5,10)))
data_summary$FIX_str = 'Fixed NK'
data_summary[data_summary$FIX == 0,]$FIX_str = 'Random NK'
data_summary$ED_str = 'Levenshtein'
data_summary[data_summary$ED == 1,]$ED_str = 'Damerau-Levenshtein'
data_summary$MUT_str = paste0('Mut. Rate: ', data_summary$MUT)
data_summary$MUT_factor = factor(data_summary$MUT_str, levels  = paste0('Mut. Rate: ', c(0.001, 0.01, 0.1)))


######## Sigmoid plots ########
for(seed in c(103)){
  print(paste0('Seed ', seed))
  for(fix in unique(data$FIX)){
    print(paste0('\tFIX ', fix))
    for(k in unique(data[data$FIX == fix,]$K)){
      print(paste0('\t\tK ', k))
      for(ed in unique(data$ED)){
        print(paste0('\t\t\tED ', ed))
        for(mut in unique(data$MUT)){
          print(paste0('\t\t\t\tMUT ', mut))
          data_condition = data[
            data$update %in% seq(0,5000, by = 100) & 
            data$seed == seed & 
            data$FIX == fix &
            data$K == k &
            data$ED == ed &
            data$MUT == mut
          ,]
          
          data_condition$plot_order = NA
          for(update in unique(data_condition$update)){
            data_read_only = data_condition[data_condition$update == update,]
            order_vec = order(data_read_only$edit_distance)
            data_condition[data_condition$update ==update,][order_vec,]$plot_order = 1:nrow(data_read_only)
          }
          
          # ggplot(data_test, aes(x = plot_order, y = edit_distance)) + 
          #   geom_line() +
          #   facet_wrap(vars(update))
          
          ggplot(data_condition, aes(x = plot_order, y = edit_distance, color = update, group = update)) + 
            geom_line(alpha = 0.5) + 
            facet_grid(rows = vars(K), cols = vars(MUT)) + 
            ggtitle(paste0('Fixed: ', fix, ', Edit distance: ', ed)) +
            ggsave(
              paste0('plots/ordered__seed_', seed, '__k_', k, '__fix_', fix, '__ed_', ed, '__mut_', mut, '.png'), 
              units = 'in', width = 6, height = 6)
        }
      }
    }
  }
}


######## Misc. from edit_distance_plot.R ########
data_test = data[
  data$update %in% seq(0,2500, by = 100) & 
  data$seed == 102 & 
  data$FIX == 1 &
  data$K == 3 &
  data$ED == 0 &
  data$MUT == 0.001
  ,]

data_test$plot_order = NA
for(update in unique(data_test$update)){
  cat(update, nrow(data_test[data_test$update == update,]), '\n')
  data_read_only = data_test[data_test$update == update,]
  order_vec = order(data_read_only$edit_distance)
  data_test[data_test$update == update,][order_vec,]$plot_order = 1:nrow(data_read_only)
  #data_test[data_test$update == update,][order_vec,]$plot_order = 1:nrow(data_read_only)
  #print(data_test[data_test$plot_order,]$edit_distance)
}

# ggplot(data_test, aes(x = plot_order, y = edit_distance)) + 
#   geom_line() +
#   facet_wrap(vars(update))

ggplot(data_test, aes(x = plot_order, y = edit_distance, color = update, group = update)) + 
  geom_line(alpha = 0.5)
