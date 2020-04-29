rm(list = ls())
library(ggplot2)
library(ggridges)
library(dplyr)

data = read.csv('./data/data_combined.csv')


######## Fixed NK - Fit and flat #########

data_plot = data[data$FIX == 1,]
# Without min/max ribbons
ggplot(data_plot, aes(x = update, y = edit_distance_mean)) +
  geom_line(alpha = 0.1, aes(group = seed)) +
  xlab('Generation') +
  ylab('Mean Edit Distance') +
  facet_grid(cols = vars(K_factor), rows = vars(MUT_factor)) +
  labs(color = 'Edit Distance Metric') +
  theme(legend.position = 'bottom') + 
  ggsave('./plots/fixed_nk.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/fixed_nk.png', units = 'in', width = 8, height = 6)
rm(data_plot)

######## Random NK ########
data_plot = data[data$FIX == 0,]
# Without min/max ribbons
ggplot(data_plot, aes(x = update, y = edit_distance_mean)) +
  geom_line(alpha = 0.1, aes(group = seed)) +
  xlab('Generation') +
  ylab('Mean Edit Distance') +
  facet_grid(cols = vars(K_factor), rows = vars(MUT_factor)) +
  labs(color = 'Edit Distance Metric') +
  theme(legend.position = 'bottom') + 
  ggsave('./plots/random_nk.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/random_nk.png', units = 'in', width = 8, height = 6)
rm(data_plot)

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
data_plot = data[data$update %in% gen_vec & data$FIX == 0 & data$MUT == 0.01 & data$ED == 0,]

ggplot(data_plot, aes(x = as.factor(update), y = edit_distance_mean, fill = as.factor(K))) + 
  geom_boxplot(position = position_dodge(0.5), width = 0.5) + 
  scale_fill_manual(values = c('#66c2a5', '#fc8d62','#8da0cb')) + 
  xlab('Generation') +
  ylab('Edit Distance (Mean of all bits)') +
  labs(fill = 'K') +
  ggsave('./plots/K_boxes.pdf', units = 'in', width = 8, height = 6) +
  ggsave('./plots/K_boxes.png', units = 'in', width = 8, height = 6)

#### Stats
data_stats = data.frame(data = matrix(nrow = 0, ncol = 5))
colnames(data_stats) = c('update', 'a', 'b', 'p_value', 'p_adjust')
for(update in gen_vec){
  data_update = data_plot[data_plot$update == update,]
  data_stats[nrow(data_stats) + 1,] = c(
    update, 
    3, 
    5, 
    wilcox.test(
      data_update[data_update$K == 3, ]$edit_distance_mean,
      data_update[data_update$K == 5,]$edit_distance_mean,
      paired = T
    )$p.value,
    1
  )
  data_stats[nrow(data_stats) + 1,] = c(
    update, 
    3, 
    10, 
    wilcox.test(
      data_update[data_update$K == 3, ]$edit_distance_mean,
      data_update[data_update$K == 10,]$edit_distance_mean,
      paired = T
    )$p.value,
    1
  )
  data_stats[nrow(data_stats) + 1,] = c(
    update, 
    5, 
    10, 
    wilcox.test(
      data_update[data_update$K == 5, ]$edit_distance_mean,
      data_update[data_update$K == 10,]$edit_distance_mean,
      paired = T
    )$p.value,
    1
  )
  data_stats[data_stats$update == update,]$p_adjust = p.adjust(data_stats[data_stats$update == update,]$p_value, 'bonferroni')
}
write.csv(data_stats, './stats/k_effect.csv')
rm(data_stats)

########## Boxplots in the fixed space ########
data_plot = data[data$update == 2500 & data$FIX == 1,]
#ggplot(data, aes(x = MUT_factor, y = edit_distance_mean)) + 
#  geom_boxplot(alpha = 0.1)
ggplot(data_plot, aes(x = edit_distance_mean, y = factor(MUT, levels = c(0.1, 0.01, 0.001)))) + 
  geom_density_ridges(
    jittered_points = TRUE,
    position = 'raincloud',
    point_alpha = 0.6, alpha = 0.7,
    scale = 0.7
  ) + 
  xlab('Edit distance') +
  ylab('Mutation rate (per locus)') + 
  ggsave('./plots/fixed_nk_raincloud.pdf', units = 'in', width = 8, height = 6) + 
  ggsave('./plots/fixed_nk_raincloud.png', units = 'in', width = 8, height = 6) 

#### Stats
data_stats = data.frame(data = matrix(nrow = 0, ncol = 4))
colnames(data_stats) = c('mut_a', 'mut_b', 'p_value', 'p_adjust')
data_update = data_plot
data_stats[nrow(data_stats) + 1,] = c(
  0.1, 
  0.01, 
  wilcox.test(
    data_update[data_update$MUT == 0.1, ]$edit_distance_mean,
    data_update[data_update$MUT == 0.01,]$edit_distance_mean,
    paired = T
  )$p.value,
  1
)
data_stats[nrow(data_stats) + 1,] = c(
  0.1, 
  0.001, 
  wilcox.test(
    data_update[data_update$MUT == 0.1, ]$edit_distance_mean,
    data_update[data_update$MUT == 0.001,]$edit_distance_mean,
    paired = T
  )$p.value,
  1
)
data_stats[nrow(data_stats) + 1,] = c(
  0.01, 
  0.001, 
  wilcox.test(
    data_update[data_update$MUT == 0.01, ]$edit_distance_mean,
    data_update[data_update$MUT == 0.001,]$edit_distance_mean,
    paired = T
  )$p.value,
  1
)
data_stats$p_adjust = p.adjust(data_stats$p_value, 'bonferroni')
write.csv(data_stats, './stats/fixed_nk_mut_effect.csv')
