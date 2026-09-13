library(deSolve)
library(ggplot2)
library(msm)

#=======================================================================
# Consumer-Resource (C-R) competition model
# Two consumers competing for one resource.
#
# Model:
# dR/dt  = I - b1*R/(R01+R)*C1 - b2*R/(R02+R)*C2
# dC1/dt = g1*b1*R/(R01+R)*C1 - d1*C1
# dC2/dt = g2*b2*R/(R02+R)*C2 - d2*C2
# State = (C1, C2, R)
# We plot the C1-C2 phase plane.

#=======================================================================
#=======================================================================

# Orbits

cr_orbits = function(state, params, t) {

  ode_rhs = function(t, state, params) {

    C1 = state[1]
    C2 = state[2]
    R  = state[3]

    I   = params$I
    b1  = params$b1
    b2  = params$b2
    R01 = params$R01
    R02 = params$R02
    g1  = params$g1
    g2  = params$g2
    d1  = params$d1
    d2  = params$d2

    dR = I - (b1*C1*R)/(R01 + R) - (b2*C2*R)/(R02 + R)
    dC1 = g1*b1*R/(R01 + R)*C1 - d1*C1
    dC2 = g2*b2*R/(R02 + R)*C2 - d2*C2
    list(c(C1 = dC1, C2 = dC2, R = dR))
  }


  # Transfer state into data.frame structure
  # irrespective of whether vector or matrix was supplied.
  if (is.null(dim(state))) {
    state = data.frame(C1 = state[1], C2 = state[2], R  = state[3])
  } else {
    state = as.data.frame(state)
  }


  orbit_structure = data.frame()

  n = nrow(state)

  for (i in 1:n) {

    current_state =
      c(C1 = state[i, "C1"],
        C2 = state[i, "C2"],
        R  = state[i, "R"])

    Traj = as.data.frame(ode( y = current_state, times = t,
                              func = ode_rhs, parms = params,
                              method = "radau"))

    Traj$orbit = i

    orbit_structure =
      rbind(orbit_structure, Traj)
  }


  return(orbit_structure)
}


#=======================================================================
# Consumer break-even resource levels

R1_star = function(params) {
  params$d1 * params$R01 /(params$g1 * params$b1 - params$d1)
}

R2_star = function(params) {
  params$d2 * params$R02 /(params$g2 * params$b2 - params$d2)
}


#=======================================================================
# C1-C2 nullclines
# We use the resource equilibrium condition:
# dR/dt = 0

# I = b1*R/(R01+R)*C1 + b2*R/(R02+R)*C2
# For C1 nullcline:
# C1 = 0 or R = R1_star
# Substituting R = R1_star into dR/dt = 0 gives:
# C1 = [I - b2*R1_star/(R02+R1_star)*C2] / [b1*R1_star/(R01+R1_star)]
# Similarly for C2.

nullcline_C1 = function(C1, params) {
  R1 = R1_star(params)
  (params$I-(params$b1*R1/(params$R01+R1))*C1)/
    ((params$b2*R1)/(params$R02+R1))

}

nullcline_C2 = function(C1, params) {
  R2 = R2_star(params)
  (params$I-(params$b1*R2/(params$R01+R2))*C1)/
    ((params$b2*R2)/(params$R02+R2))
}


#=======================================================================
# Phase Portrait

plot_phase_portrait_CR = function(states,params,
                               t = seq(0, 50, length.out = 200)) {

  orbits = cr_orbits(states, params, t)

  R1 = R1_star(params)
  R2 = R2_star(params)

  nc1_C1_intercept = params$I/(params$b1*R1/(params$R01+R1))
  nc2_C1_intercept = params$I/(params$b1*R2/(params$R01+R2))
  nc1_C2_intercept = params$I/(params$b2*R1/(params$R02+R1))
  nc2_C2_intercept = params$I/(params$b2*R2/(params$R02+R2))

  max_C1_needed = max(nc1_C1_intercept, nc2_C1_intercept)
  max_C2_needed = max(nc1_C2_intercept, nc2_C2_intercept)

  C1_max = 1.15 * max_C1_needed
  C2_max = 1.15 * max_C2_needed

  C1_seq = seq(0, C1_max, length.out = 200)

  nc1 = data.frame(C1 = C1_seq, C2 = nullcline_C1(C1_seq, params))
  nc2 = data.frame(C1 = C1_seq, C2 = nullcline_C2(C1_seq, params))

  nc1 = nc1[nc1$C1 >= 0,]
  nc2 = nc2[nc2$C2 >= 0,]


  p = ggplot() +
    geom_line(data = nc1, aes(C1, C2, colour = 'nC1'),
              linewidth = 0.7) +
    geom_line(data = nc2, aes(C1, C2, colour = 'nC2'),
              linewidth = 0.7) +
    geom_path(data = orbits, aes(C1, C2, group = orbit),
              colour = "steelblue", linewidth = 0.3) +
    scale_color_manual(name="Nullclines",
                       values=c('nC1'="firebrick", 'nC2'="darkgreen")) +
    coord_cartesian( xlim = c(0, C1_max), ylim = c(0, C2_max)) +
    labs(x = expression(C[1]),y = expression(C[2]),
      title = "Consumer-Resource competition",
      subtitle = paste("R1* = ", round(R1, 6),", R2* =", round(R2, 6))) +
    theme_minimal()

  return(p)
}


#=======================================================================

states_init_generate = function(params, n){
  R1 = R1_star(params)
  R2 = R2_star(params)

  nc1_C1_intercept = params$I/(params$b1*R1/(params$R01+R1))
  nc2_C1_intercept = params$I/(params$b1*R2/(params$R01+R2))
  nc1_C2_intercept = params$I/(params$b2*R1/(params$R02+R1))
  nc2_C2_intercept = params$I/(params$b2*R2/(params$R02+R2))

  max_C1_needed = max(nc1_C1_intercept, nc2_C1_intercept)
  max_C2_needed = max(nc1_C2_intercept, nc2_C2_intercept)

  C1_max = 1.15 * max_C1_needed
  C2_max = 1.15 * max_C2_needed

  C1_init = c(runif(n, 0, C1_max), 0.2)
  C2_init = c(runif(n, 0, C2_max), 0.3)
  R_init = rep(10, length(C1_init))

  C_init = cbind(C1_init, C2_init, R_init)
  colnames(C_init) = c("C1", "C2", "R")

  return(C_init)
}
#=======================================================================
#=======================================================================
# Trial


# params = list(
#   I = 12,
#   b1 = 4,
#   b2 = 3,
#   R01 = 5,
#   R02 = 5,
#   g1 = 1.4,
#   g2 = 1,
#   d1 = 0.5,
#   d2 = 0.8
# )
#
# states = states_init_generate(params, 10)
# plot_phase_portrait_CR(states = states, params = params)

