rm(list = ls())
library(ggplot2)

#w_0 = 0.524493
w_0 = 0.503661
#w_1 = 0.499336
w_1 = 0.486666
w_2 = 0.487320
#w_2 =0.503585

alpha = -1 * log(w_1/w_0)
beta = log2((-1 * log(w_2/w_0)) / alpha)

point_vec = c(w_0, w_1, w_2)
mod_point_vec = log(point_vec / w_0)
data_points = data.frame(x = 0:2, y = point_vec, y_mod = mod_point_vec)
x = seq(0, 3, by = 0.01)
y = -1 * alpha * (x^beta)
data_line = data.frame(x = x, y = y)

ggplot(data_line, aes(x, y)) + 
  geom_line() +
  geom_point(data = data_points, aes(x, y_mod))
  