rm(list = ls())
library(ggplot2)
library(ggridges)
library(dplyr)

data = read.csv('./data/data_combined.csv')
data$VEL_str = paste0('Velocity: ', data$VEL)
data$VEL_factor = factor(data$VEL_str, (paste0('Velocity: ', c(0.1, 0.05, 0.01))))

######## Random NK ########
for(vel in unique(data$VEL)){
  data_plot = data[data$FIX == 0 & data$VEL == vel,]
  # Without min/max ribbons
  ggplot(data_plot, aes(x = update, y = edit_distance_mean)) +
    geom_line(alpha = 0.1, aes(group = seed)) +
    xlab('Generation') +
    ylab('Mean Edit Distance') +
    ggtitle(paste0('Treadmill, Velocity: ', vel)) +
    facet_grid(cols = vars(factor(K_factor, levels = paste0('K = ', c(3,5,10)))), rows = vars(MUT_factor)) +
    labs(color = 'Edit Distance Metric') +
    theme(legend.position = 'bottom') + 
    ggsave(paste0('./plots/random_nk_vel_', vel, '.pdf'), units = 'in', width = 8, height = 6) +
    ggsave(paste0('./plots/random_nk_vel_', vel, '.png'), units = 'in', width = 8, height = 6)
  rm(data_plot)
}
for(k in unique(data$K)){
  data_plot = data[data$FIX == 0 & data$K == k,]
  # Without min/max ribbons
  ggplot(data_plot, aes(x = update, y = edit_distance_mean)) +
    geom_line(alpha = 0.1, aes(group = seed)) +
    xlab('Generation') +
    ylab('Mean Edit Distance') +
    ggtitle(paste0('Treadmill, K: ', k)) +
    facet_grid(cols = vars(factor(VEL_factor, levels = paste0('Velocity: ', c(0.1, 0.05, 0.01)))), rows = vars(MUT_factor)) +
    labs(color = 'Edit Distance Metric') +
    theme(legend.position = 'bottom') + 
    ggsave(paste0('./plots/random_nk_k_', k, '.pdf'), units = 'in', width = 8, height = 6) +
    ggsave(paste0('./plots/random_nk_k_', k, '.png'), units = 'in', width = 8, height = 6)
  rm(data_plot)
}

######## Effect of K (Random) ########
# Show difference in edit distance as K changes in random tables

#### Lines
data_plot = data[data$FIX == 0 & data$MUT == 0.01 & data$ED == 0,]
data_plot$group = paste(data_plot$seed, data_plot$K)
ggplot(data_plot, aes(x = update, y = edit_distance_mean, color = as.factor(K), group = group)) + 
  geom_line(alpha = 0.3) +
  xlab('Generation') + 
  ylab('Edit Distance (Mean of all bits)') + 
  labs(color = 'K') + 
  guides(alpha = F) +
  scale_color_manual(values = c('#66c2a5', '#fc8d62','#8da0cb')) + 
  ggsave('./plots/K_lines.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/K_lines.png', units = 'in', width = 8, height = 6)

data_plot_grouped = dplyr::group_by(data_plot, update, FIX, K, ED , MUT, K_factor, MUT_str) 
data_plot_summary = dplyr::summarize(data_plot_grouped, edit_distance_grand_mean = mean(edit_distance_mean))
ggplot(data_plot, aes(x = update, y = edit_distance_mean, color = as.factor(K))) + 
  geom_line(aes( group = group), alpha = 0.3) +
  xlab('Generation') + 
  ylab('Edit Distance (Mean of all bits)') + 
  labs(color = 'K') + 
  guides(alpha = F) +
  geom_line(data = data_plot_summary, aes(y = edit_distance_grand_mean), size = 1.5) +
  scale_color_manual(values = c('#66c2a5', '#fc8d62','#8da0cb')) + 
  ggsave('./plots/K_lines_with_grand_mean.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/K_lines_with_grand_mean.png', units = 'in', width = 8, height = 6)

rm(data_plot)
rm(data_plot_grouped)
rm(data_plot_summary)

#### Boxes
gen_vec = c(500, 1000, 1500, 2000, 2500)
data_plot = data[data$update %in% gen_vec & data$FIX == 0 & data$MUT == 0.001 & data$ED == 1 & data$K == 3,]

ggplot(data_plot, aes(x = as.factor(update), y = edit_distance_mean, fill = as.factor(VEL))) + 
  geom_boxplot(position = position_dodge(0.5), width = 0.5) + 
  scale_fill_manual(values = c('#66c2a5', '#fc8d62','#8da0cb')) + 
  xlab('Generation') +
  ylab('Edit Distance (Mean of all bits)') +
  labs(fill = 'Velocity') +
  ggsave('./plots/vel_boxes_k_3.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/vel_boxes_k_3.png', units = 'in', width = 8, height = 6)

#### Stats
data_stats = data.frame(data = matrix(nrow = 0, ncol = 5))
colnames(data_stats) = c('update', 'a', 'b', 'p_value', 'p_adjust')
for(update in gen_vec){
  data_update = data_plot[data_plot$update == update,]
  data_stats[nrow(data_stats) + 1,] = c(
    update, 
    0.01, 
    0.05, 
    wilcox.test(
      data_update[data_update$VEL == 0.01, ]$edit_distance_mean,
      data_update[data_update$VEL == 0.05,]$edit_distance_mean,
      paired = T
    )$p.value,
    1
  )
  data_stats[nrow(data_stats) + 1,] = c(
    update, 
    0.01, 
    0.1, 
    wilcox.test(
      data_update[data_update$VEL == 0.01, ]$edit_distance_mean,
      data_update[data_update$VEL == 0.1,]$edit_distance_mean,
      paired = T
    )$p.value,
    1
  )
  data_stats[nrow(data_stats) + 1,] = c(
    update, 
    0.05, 
    0.1, 
    wilcox.test(
      data_update[data_update$VEL == 0.05, ]$edit_distance_mean,
      data_update[data_update$VEL == 0.1,]$edit_distance_mean,
      paired = T
    )$p.value,
    1
  )
  data_stats[data_stats$update == update,]$p_adjust = p.adjust(data_stats[data_stats$update == update,]$p_value, 'bonferroni')
}
write.csv(data_stats, './stats/vel_effect_k_3.csv')
rm(data_stats)