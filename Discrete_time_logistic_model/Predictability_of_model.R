source("Logistic_growth_model.R")

#What we mean by predictability of a model?
# It is about, for initial measurment with small error can we approximately
# predict the future states?.
# if for two very close initial states, if we get drastically different
# trajectories then we call the system unpredictable as error measuring initial
# condition is practically unavoidable.


# If system is not predictable that means it may show the butterfly effect:
# A flap of butterfly in Bombay may cause tornado in Tokio.
# small initial changes can cause huge, unpredictable outcomes.

#(Precise epsilon delta definations are in README)

#===============================================================================


R_m = c(0.5, 1.5, 2.5, 3.2, 3.6, 4)
t1 = 1:100 # better visualization.
t2 = 1: 1000 # for long run.
error = c(0.1, 0.01, 0.00001)


# Plot two graphs
# 1. {phi_t(x_0 + errors} against t. lets call them 'error trajectories'.
# This plots time trajectory of very nearby initial conditions,
# so we can see how they evolve.
# 2.  Et = {Xt_error} - Xt, Plot Et against t. lets call it 'divergence of E_t'
# This shows how initial error grows decays over time.

Error_trajectories = function(x_0, R_m, t, eps = error){

  if (length(eps) == 3){
    colors = c('blue','darkgreen','red','black')
  } else {colors = rainbow(length(eps))} # For general any error set.

  # Base plot for error trajactories
  plot(1, type = 'n', xlim = c(1, length(t)), ylim = c(0, 1) ,xlab = "t",
       ylab = expression(x[t]), main = bquote("Error trajectories for " * R[0]
                                              == .(R_m)~' at '~ x[0] == .(x_0)))


  Xt_base = popn_trajectory(x_0, R_m, t) # Trajectory of x_0
  Et_matrix = matrix(NA, nrow = length(t) , ncol = length(eps))
  # Matrix with i-th column consist of vectors Et.
  # Where Et majors error in two trajectories over time.
  #each row corresponds to different values of error form list eps.

  for (i in 1:length(eps)){
    Xt_error = popn_trajectory((x_0+eps[i]), R_m, t)
    # Calculate trajectory for initial population of x_0+error under function f.
    lines(t, Xt_error, col= colors[i], lwd=1.2)
    # Draws Xt_error.
    Et_matrix[ , i] = Xt_error - Xt_base # Records each column
  }
  lines(t, Xt_base, col= colors[length(eps)+1], lwd=1.2)
  # Draws trajectory of x_0 on the same base graph.

  legend("topright", legend = paste(' error = ', c(eps, 0)), col = colors,
         lty = 1, lwd = 2, cex = 0.5)


  # another base graph to plot Et vs t
  plot(1, type = 'n', xlim = c(1, length(t)), ylim = c(-1, 1), xlab = "t", ylab
       = expression(E[t]),
       main = bquote("Divergence of " * E[t] * ' for ' * R[0] == .(R_m) ~
                       "at"~ x[0] == .(x_0)))

  for (i in 1:length(eps)){
    lines(t, Et_matrix[, i], col = colors[i], lwd = 1.2)
    #draws Et for each error value on base graph.
  }

  legend("topright", legend = paste('error = ', eps), col = colors, lty = 1,
         lwd = 2, cex = 0.5)
}

# Plot all the graphs for list of x_0, list of R_m for given dynamic
# fnunction f.
Plot_error_trajectories =function(X_0_list = c(0.1, 0.7), R_m_list = R_m,t= t1,
                                  eps = error){
  for( i in X_0_list){
    for (j in R_m_list){
      Error_trajectories(i, j, t, eps)
    }
  }
}



#Plot_error_trajectories(c(0.01, 0.9, 0.1, 0.11))
#Plot_error_trajectories(X_0_list = c(0.1), R_m_list = R_m, t= t1, eps= error)

# We noticed that R_m = 0.5, 1.5, 2.5 are predictable but 3.6 and 3.4 are not
# At R_m = 3.2, system goes out of phase for at x_0= 0.1, error=0.1.
# But for smaller error it is predictable.


# Def chaos: (speaking loosely) small change is initial condition leads
# to exponential divergence in trajectories.


