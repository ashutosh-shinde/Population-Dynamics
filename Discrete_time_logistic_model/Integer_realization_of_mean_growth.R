#source('Discrete_time_logistic_growth_model/Logistic_growth_model.R')
#source('Discrete_time_logistic_growth_model/Predictability_of_model.R')


#Simulates one stochastic realization of discrete-time logistic growth.
popn_trajectory_int = function(N_0, R_m, K, time_steps, distribution = c(
                               "poisson", "negbinom", "normal", "none"), disp =1
                               , rtn = 'density' )
  {
  # N0    : initial population size
  # Rm    : Max growth rate popuation can attain.
  # K     : population size when R_t becomes 0 => extinction.
  # K = 1000 is default => N_0 = 100 ~ X_0 = 0.1
  # time_steps : number of steps to simulate
  # distribution : noise model for demographic stochasticity
  # disp  : dispersion parameter, used only for negbinom (size argument)
  # rtn   : put rtn = Raw to get Xt vector, otherwise fn will returen Xt.

  distribution = match.arg(distribution)
  # user-provided character string matches one of several pre-defined
  # valid distribution choices.
  Nt = numeric(time_steps)
  Xt = numeric(time_steps)
  Nt[1] = N_0
  Xt[1] = N_0/K


  for (i in 2:time_steps) {
    expd_Nt = R_m * Nt[i - 1] * (1 - Nt[i - 1] / K)
    expd_Nt = max(expd_Nt, 0)   # To make sure expd_Nt > 0.

    # Only different step from previous model. Here we actually choose N_t from
    # random sampling under perticular distribution with mean [N_t-1 * R_m].
    # Insted of taking it equal to expected value. Well here we are ignoring
    # many stocastic factors and just taking integer of N.
    # We will almost recover Previouse model when we take K <- too large.
    Nt[i] = switch(distribution,
                    poisson  = rpois(1, lambda = expd_Nt),
                    negbinom = rnbinom(1, mu = expd_Nt, size = disp),
                    normal   = max(0, round(rnorm(1, mean = expd_Nt, sd =
                                                    sqrt(expd_Nt)))),
                    none     = expd_Nt
    ) # switch : choose value of N_t[i] from appropriet distrubution.

    Xt[i] = Nt[i]/K
    # Store density for same trajectory
  }

  if(rtn == 'raw') return(Nt)
  else return(Xt)
}
#set.seed(1)
# Uses same starting point for random value generating function so you get
# same result each time.

#===============================================================================

# Code here onwards is experimented out of my own curiosity; and it has some
# serious limitations and mathematical errors.
#===============================================================================


#same code as befor just changes the trajetory function.
popn_plots_for_R_m_int = function(N_0_list, R_m, K, time_steps, distribution)
  {
  colors = rainbow(length(N_0_list))
  eqpt = rep((1 - 1/R_m), times = time_steps)
  # This is constant line of equilibrium population size, for respective R_m

  plot(1, type ="n", xlim = c(1, time_steps), ylim = c(0,1),
       xlab= "t", ylab = expression(X[t]), main = paste(distribution,
       "logistic Map with K =", K, " R_m =", R_m), cex.main = 0.7)
  # Set the empty base plot.

  # Each loop plots population size values for each x_0 in X_0.
  for (i in 1:length(N_0_list)) {
    Xt = popn_trajectory_int(N_0 = N_0_list[i] , R_m = R_m, K = K, time_steps =
                               time_steps, distribution = distribution)
    # Set Xt for each n_0 in N_0
    lines(1:time_steps, Xt, col= colors[i], lwd=1.2)
    # Draw population size at t for respective x_0
    lines(1:time_steps, eqpt, col = "black")
    # Draw a black line equlibrium value of non-stocastic(!) model
  }

  X_0 = N_0_list/K
  legend("topright", legend = paste("x0 =", X_0), col = colors, lty = 1,
         lwd = 2, cex = 0.4) # To see what color represent what x_0.
}


plot_each_R_m_int = function(N_0_list, R_m_list, K,
                         time_steps, distribution)
  {
  for (R_m in R_m_list){
    popn_plots_for_R_m_int(N_0_list = N_0_list, R_m = R_m, K = K, time_steps =
                            time_steps, distribution = distribution)
  }
}




N_0_list_1 =
R_m_list_1 =

# plot_each_R_m_int(N_0_list = c(100, 200, 500, 800), R_m_list = c(0.5, 1.5, 2.5, 3.2,
#              3.6, 3.9), K = 1000, time_steps = 100, distribution = "poisson")


#===============================================================================
#bifurcation: same code

Plot_Bifurcation_int = function(N_0, R_m_list, K, time_steps, distribution){

  plot(NULL, xlim = c( min(R_m_list), max(R_m_list)), ylim  = c(0,1),
       xlab = expression(R_m), ylab = "Attractor points",
       main = "Bifurcation Diagram")

  for (R_m in R_m_list) {
    Xt = popn_trajectory_int(N_0 = N_0, R_m = R_m, K = K, time_steps =
                              time_steps, distribution = distribution)
    last_200 = tail(Xt, 200)        # Keep last 200 values
    points(rep(R_m, 200), last_200, pch = ".", col = "black")
  }
}

# Plot_Bifurcation_int(N_0 = 100, R_m_list = seq(0, 3.9, 0.01), K = 1000, time_steps
# = 1000, distribution = "poisson")
