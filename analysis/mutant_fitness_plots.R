rm(list = ls())
library(ggplot2)
library(dplyr)

data = read.csv('./data/data_mutant_fitness__seed_all.csv')

data_plot = data[data$FIX == 0,]
data_plot = data_plot[data_plot$update %in% c(500,1000,1500,2000,2500),]

ggplot(data_plot, aes(x = update, y = alpha_el, color = as.factor(update))) + 
  geom_hline(aes(yintercept = 0), alpha = 0.2) +
  geom_point(alpha = 0.2) + 
  ggtitle('Alpha (Elena and Lenski)') +
  facet_grid(cols = vars(MUT), rows = vars(K)) 
 
ggplot(data_plot, aes(x = update, y = beta_el, color = seed)) + 
  geom_hline(aes(yintercept = 0), alpha  = 0.2) +
  geom_jitter(alpha = 0.2) + 
  ggtitle('Beta (Elena and Lenski)') +
  facet_grid(cols = vars(MUT), rows = vars(K))

ggplot(data_plot, aes(x = as.factor(update), y = beta_el, fill = as.factor(seed))) +
  geom_hline(aes(yintercept = 0), alpha  = 0.2) +
  geom_boxplot() +#position = position_dodge2(1000)) + 
  ggtitle('Beta (Elena and Lenski)') +
  facet_grid(cols = vars(MUT), rows = vars(K))



# Draw seeds separately

data_grouped = dplyr::group_by(data, 
                               update, FIX,K, MUT, ED, seed)
data_summary = dplyr::summarize(data_grouped, 
                                alpha_el_mean = mean(alpha_el),
                                beta_el_mean = mean(beta_el))

data_summary$group = paste0(data_summary$K, '_', data_summary$seed)
ggplot(data_summary, aes(x = update, y = alpha_el_mean, color = as.factor(K), group = group)) +
  geom_line(alpha = 0.4) + 
  facet_grid(rows = vars(FIX), cols = vars(MUT))

ggplot(data_summary, aes(x = update, y = beta_el_mean, color = as.factor(K), group = group)) +
  geom_line(alpha = 0.4) + 
  facet_grid(rows = vars(FIX), cols = vars(MUT))


# Average over seeds
data_grouped = dplyr::group_by(data, 
                               update, FIX,K, MUT, ED)
data_summary = dplyr::summarize(data_grouped, 
                                alpha_el_mean = mean(alpha_el),
                                beta_el_mean = mean(beta_el))

ggplot(data_summary, aes(x = update, y = alpha_el_mean, color = as.factor(K))) +
  geom_line(aes(linetype = as.factor(MUT))) + 
  facet_grid(rows = vars(FIX))

ggplot(data_summary, aes(x = update, y = beta_el_mean, color = as.factor(K))) +
  geom_line(aes(linetype = as.factor(MUT))) + 
  facet_grid(rows = vars(FIX))



# 
# 
# ggplot(data_test, aes(x = alpha_el), fill = 'blue') + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_histogram(binwidth = 0.005) + 
#   ggtitle('Alpha (Elena and Lenski)')
#   facet_grid(rows = vars(MUT), cols = vars(update)) 
# 
# ggplot(data_test[data_test$update == max(data_test$update),], aes(x = beta_el), fill = 'blue') + 
#   geom_vline(aes(xintercept = 0)) +
#   geom_histogram(binwidth = 0.005) + 
#   ggtitle('Beta (Elena and Lenski)') +
#   facet_grid(cols = vars(MUT), rows = vars(K))
# 
# data_test$one_step_fraction = data_test$one_step_fitness_avg / data_test$fitness_avg
# data_test$two_step_fraction = data_test$two_step_fitness_avg / data_test$fitness_avg
# min_val = min(min(data_test$one_step_fraction), min(data_test$two_step_fraction))
# max_val = max(max(data_test$one_step_fraction), max(data_test$two_step_fraction))
# x = seq(min_val, max_val, by = 0.001)
# data_line = data.frame(x = x, y = x)
# ggplot(data_test, aes(x = one_step_fraction, y = two_step_fraction, color = as.factor(update))) +
#   geom_line(data_line, mapping = aes(x = x, y = y)) +
#   geom_point(alpha = 0.4) + 
#   facet_grid(rows = vars(MUT), cols = vars(K))
# 
# data_test$diff_expected  = (1 - log(data_test$one_step_fitness_avg / data_test$fitness_avg)) - (1 - 2 * data_test$alpha_el)
# ggplot(data_test[data_test$update == 2000,], aes(x = diff_expected)) + 
#   geom_histogram() + 
#   facet_grid(rows = vars(MUT), cols = vars(K))
# 
# 
# 
# # OLD
# 
# #data_plot = data_originals[!is.na(data_originals$beta),]
# data_plot = data_originals[!is.na(data_originals$one_step_fitness_avg),]
# data_plot$one_step_fraction = data_plot$one_step_fitness_avg / data_plot$fitness_avg
# data_plot$two_step_fraction = data_plot$two_step_fitness_avg / data_plot$fitness_avg
# min_val = min(min(data_plot$one_step_fraction), min(data_plot$two_step_fraction))
# max_val = max(max(data_plot$one_step_fraction), max(data_plot$two_step_fraction))
# x = seq(min_val, max_val, by = 0.001)
# data_line = data.frame(x = x, y = x)
# ggplot(data_plot, aes(x = one_step_fraction, y = two_step_fraction, color = as.factor(update))) +
#   geom_line(data_line, mapping = aes(x = x, y = y)) +
#   geom_point(alpha = 0.4) + 
#   facet_grid(rows = vars(MUT))
# 
# data_plot = data_plot[data_plot$MUT == 0.001,]
# ggplot(data_plot, aes(x = update, y = beta, color = as.factor(update))) +
#   geom_hline(aes(yintercept = 1)) +
#   geom_point(alpha = 0.4)
# 
# ggplot(data_plot, aes(x = update, y = alpha, color = as.factor(update))) +
#   geom_point(alpha = 0.4)
# 
# 
# 
# # Beta > 1
# data_solo = data_plot[data_plot$update == 2000 & data_plot$beta > 1,]
# x = seq(0, 3, by = 0.01)
# y = -1 * data_solo$alpha * (x^data_solo$beta)
# data_line = data.frame(x = x, y = y)
# 
# data_points = data.frame(x = c(0, 1, 2), y = c(0, log(data_solo$one_step_fraction), log(data_solo$two_step_fraction)))
# ggplot(data_line, aes(x, y)) + 
#   geom_line() +
#   geom_point(data = data_points, aes(x, y))
# 
# # Beta < 1
# data_solo = data_plot[data_plot$update == 2000 & data_plot$beta == min(data_plot[data_plot$update == 2000,]$beta),][1,]
# x = seq(0, 3, by = 0.01)
# y = -1 * data_solo$alpha * (x^data_solo$beta)
# data_line = data.frame(x = x, y = y)
# 
# data_points = data.frame(x = c(0, 1, 2), y = c(0, log(data_solo$one_step_fraction), log(data_solo$two_step_fraction)))
# ggplot(data_line, aes(x, y)) + 
#   geom_line() +
#   geom_point(data = data_points, aes(x, y))
# 
# # Two-step avg. fitness > one-step avg. fitness
# data_solo = data_plot[(data_plot$one_step_fitness_avg < data_plot$two_step_fitness_avg),][1,]
# x = seq(0, 3, by = 0.01)
# y = -1 * data_solo$alpha * (x^data_solo$beta)
# data_line = data.frame(x = x, y = y)
# 
# data_points = data.frame(x = c(0, 1, 2), y = c(0, log(data_solo$one_step_fraction), log(data_solo$two_step_fraction)))
# ggplot(data_line, aes(x, y)) + 
#   geom_line() +
#   geom_point(data = data_points, aes(x, y))