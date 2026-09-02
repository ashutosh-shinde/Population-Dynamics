logistic_growth = function(x0, r, delta_t, t_max) {
  # x0: initial population density, in [0, 1]
  # r: Max growth rate
  # delta_t: small time step
  # t_max: total time to simulate

  time_steps = seq(0, t_max, by = delta_t)
  x = numeric(length(time_steps))
  x[1] = x0

  for (i in 2:length(time_steps)) {
    x_prev = x[i - 1]
    r_t = r * (1 - x_prev) # growth rate density: b_t - d_t
    x[i] = x_prev + r_t * x_prev * delta_t
    # r_t * delta_t = per capita growth over the interval delta_t
  }

  data.frame(time = time_steps, x = x)
}


#trajectory = logistic_growth(x0 = 0.23, r = 0.11, delta_t = 0.01,
#                             t_max = 50)
#plot(trajectory$time, trajectory$x, type = "l", xlab = "Time", ylab =
#    "Population density (x)", main = "Continuous Time Logistic Growth")
