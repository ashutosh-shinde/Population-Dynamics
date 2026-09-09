# Set Harvester functional As one of following.
harvest_fn = function(x, type, c = NULL, h0 = NULL, x0 = NULL,
                      theta = NULL)
  {
  switch(as.character(type),
         "0" = rep(c, length(x)), # constant
         "1" = h0 * x, # linear (Holling I)
         "2" = (h0 * x) / (x0 + x), # Holling II
         "3" = (h0 * x^theta) / (x0^theta + x^theta), # Holling III
         stop("type must be 0, 1, 2, or 3")
  )
}


Popn_trajectory_holling = function(x_init, r_m, delta_t, t_max,
                            harvest_type = NULL, c = NULL, h0 = NULL,
                            x0 = NULL, theta = NULL) {
  # x_init: initial state population density, in [0, 1]
  # r_m: intrinsic growth rate
  # delta_t: very small time step
  # t_max: total time to simulate
  # harvest_type: NULL (no harvesting) or 0/1/2/3
  # c, h0, x0, theta: harvester parameters (as needed by harvest_type)

  # Note : in other codes I have taken x_0 as initial state of
  # population-size.
  # here we have taken harvester rate rate x0 is half saturation
  # density.


  time_steps = seq(0, t_max, by = delta_t)
  x = numeric(length(time_steps))
  x[1] = x_init

  for (i in 2:length(time_steps)) {
    x_prev = x[i - 1]
    r_t = r_m * (1 - x_prev)

    h_t = if (is.null(harvest_type)) {
      0
    } else {
      harvest_fn(x_prev, harvest_type, c = c, h0 = h0, x0 = x0,
                 theta = theta)
    }
    # h_t is total harvest rate density, it is not a percapita value
    # like r_t.

    x[i] = x_prev + (r_t * x_prev - h_t) * delta_t
    x[i] = max(x[i], 0)  # prevent x_0 < 0, NBS}
  }
  data.frame(time = time_steps, x = x)
}


#Traj = Popn_trajectory_holling(x_init = 0.1, r_m = 1.2, delta_t = 0.01,
#                               t_max = 50, harvest_type = 2,
#                               h0 = 0.3, x0 = 0.2)

# plot(Traj$time, Traj$x, type = "l", xlab = "Time", ylab =
#     "Population density (x)", main = "Logistic Growth with Holling
#     Type II Harvesting")


