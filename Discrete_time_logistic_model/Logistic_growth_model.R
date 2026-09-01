# Discrete logistic population growth model.

# N_t : population at time t., N_0 = initial population
# R_t: growth rate for the population N_t.(Composed of births and deaths.)
# Note that:  R_t = births-deaths + 1.

# we are interested in two scenarios: 
#1. how population changes over t for different R_m and different x_0.
#2. when population attains equilibrium x*, that is when x_(t+1) = x_t = x* 
#and how this equilibrium changes over different R_m and x_0.

# inputs:
t1 = 1:100 # Use for population size oscilations over different X_0.
t2 = 1:1000 # Use for bifurcation 
# As we are taking last 200 values. This is to get cleaner subsequential limits.

X_0 = c(0, 0.01, 0.2, 0.4, 0.5, 0.7, 0.9, 0.99, 1) # To check for diff X_0. 
X_0_1 = c(0, 0.1, 1) # When you don't want a chaotic graph.  
# The choice initial x_0 is immaterial over long run, except x_0 =  0, 1

R_m_list_b = seq(0, 4, 0.01)  # Use for final bifurcation diagram
R_m_list_1 = c(0.5, 1.5, 2.5, 3.2, 3.6, 4)  
# Use this to see overall how population size changes over different R_m 
R_m_list_2 = seq(2.95, 3.5, 0.05) # Use this to see how populations size starts
# to oscillates and how that oscillations changes with R_m.


popn_trajectory = function(x_0, R_m, time_steps){ 
  Xt = numeric(length(time_steps)) 
  # Vector to store all vlaues of population trajectory
  Xt[1]=x_0 # Set initial population
  for (i in 2:length(time_steps)){
    Xt[i] = R_m * Xt[i-1] * (1 - Xt[i-1]) # Calculating all X_t
  }
  Xt # Vector consisting of population size at all time steps,
}


# For given R_m, this function plot the population size over all t 
# For each x_0 in vector X_0. on the same sheet.
popn_plots_for_R_m = function(R_m, time_step){
  colors = rainbow(length(X_0))
  eqpt = rep((1 - 1/R_m), times = length(time_step))
  # This is constant line of equilibrium population size, for respective R_m
  
  plot(1, type ="n", xlim = c(1, length(time_step)), ylim = c(0,1),
       xlab= "t", ylab = expression(X[t]),  
       main = paste("Logistic Map, R_m =", R_m))  # Set the empty base plot.
  
  # Each loop plots population size values for each x_0 in X_0.
  for ( i in 1:length(X_0)) { 
    Xt = popn_trajectory(X_0[i], R_m, time_step) # Set Xt for each x_0 in X_0
    lines(time_step, Xt, col= colors[i], lwd=1.2) 
    # Draw population size at t for respective x_0
    lines(time_step, eqpt, col = "black") 
    # draw a black constant line of Equlibrium value.
  }
  
  legend("topright", legend = paste("x0 =", X_0), col = colors, lty = 1,
         lwd = 2, cex = 0.4) # To see what color represent what x_0.
}




# Different graphs of population size time series for each R_m in R_m_list1.
plot_each_R_m = function(R_m_list){
  for (R_m in R_m_list){
    popn_plots_for_R_m(R_m, t1)
  }
}
#plot_each_R_m(R_m_list_1)

#===============================================================================
# Part 2
#Eqilibrium

# You can see that for R_m in [0, 1] there is only 1 equilibrium, which is 0.
# and all x_0 converge there.
# Also for R_m in [1, 3) all x_0 except 0, 1 attains value (1-1/R_m).
# which is an equilibrium indeed.
# For x_0 = 0, 1 Xt converge to 0. For x_0 = 0, it's trivial.
# And for x_0 = 1 R_m becomes Zero.


# Plot x_eq over R.
Eqb = function(R_m_list_1){
  1 - 1/R_m_list_1
}
plot_eqb_curve = function(){
  curve(Eqb, from = 1, to= 4, xlab= expression(R[0]), ylab= 'x*', col= "blue")
}


# Stability of equilibrium: 


# Lets plot the subsequential limits of x_t over different R_m.
# Since it is not so easy to find limits analytically we will try to plot
# the tail of values of Xt, so if Xt oscillates it shall give us the 
# bifurcation graph.
# As we have observed in last plots:
# Till 3 there is one single limit point at which population time series
# converge.

Plot_Bifurcation = function(x_0){
  plot(NULL, xlim = c( min(R_m_list_b), max(R_m_list_b)), ylim  = c(0,1), 
       xlab = expression(R_m), ylab = "Attractor points",
       main = "Bifurcation Diagram")
  
  for (R_m in R_m_list_b) {
    Xt = popn_trajectory(x_0, R_m, t2)          # Run the logistic map
    last_200 = tail(Xt, 200)        # Keep last 200 values
    points(rep(R_m, 200), last_200, pch = ".", col = "black")
  }
}

#Plot_Bifurcation(0.1)
