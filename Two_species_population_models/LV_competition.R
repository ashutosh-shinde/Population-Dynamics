library(deSolve)
library(ggplot2)
library(msm)
#=======================================================================

# Lotka-Volterra two-species competition

# Model:
#  dN1/dt = r1*N1*(1 - (N1 + a21*N2)/K1)
#  dN2/dt = r2*N2*(1 - (N2 + a12*N1)/K2)

# states can be matrix(multiple initial condition) or vector
# params is always a named list: list(r1=, r2=, K1=, K2=, a12=, a21=)
# t is time vector to evaluate ode. e.g. seq(0, 50, length.out= 200)
#=======================================================================

lv_orbits = function(state, params, t) {

  ode_rhs = function(t, state, params){
    with(as.list(c(state, params)), {
      dN1 = r1 * N1 * (1 - (N1 + a21 * N2) / K1)
      dN2 = r2 * N2 * (1 - (N2 + a12 * N1) / K2)
      list(c(N1 = dN1, N2 = dN2))
    })
  }

  # Transfer state into data.frame structure
  # Irrespective of initial state was passed as vector or matrix.
  if (is.null(dim(state))) {
    state = data.frame(N1 = state[1], N2 = state[2])
  } else {
    state = as.data.frame(state)
  }

  orbit_structure = data.frame()

  n = nrow(state)

  for (i in 1:n) {
    current_state = c(N1 = state[i, "N1"], N2 = state[i, "N2"])

    Traj = as.data.frame(ode(y = current_state, times = t,
                            func = ode_rhs, parms = params,
                            method = "radau"))
    # Returns matrix with coloumn (t, N1, N2) with each row storing
    # state at time t.

    Traj$orbit = i
    #set new coloum orbit for index

    orbit_structure = rbind(orbit_structure, Traj)
    # Multiple traj stacked
  }

  return(orbit_structure)
}


#=======================================================================
# Nullclines
nullcline_1 = function(N1, K1, a21) (K1 - N1) / a21
nullcline_2 = function(N1, K2, a12) K2 - a12 * N1
# Note : We are treating both nullclines as function of N1.
# i.e. It gives values of N2 for each input vale of N1. on phase plane


#=======================================================================
# stability
# Fixed points assumed directly: E0, E1, E2, E3.
# existance of E3 is cheaked.
# Jacobian and eigenvalues are computed at fixed points.
# Fixed points are classified accordingly.

lv_stability = function(params) {

  r1 = params$r1
  r2 = params$r2
  K1 = params$K1
  K2 = params$K2
  a12 = params$a12
  a21 = params$a21

  E0 = c(N1 = 0, N2 = 0)
  E1 = c(N1 = K1, N2 = 0)
  E2 = c(N1 = 0, N2 = K2)

  # check if E3 exist, i.e.
  # 1. nullclines do intersect determinanat !=0, so unique soln exist.
  # 2. intersection point is in positive quadrant.
  D = 1 - a12 * a21

  if (abs(D) > 1e-12) {
    E3 = c(N1 = (K1 - a21 * K2) / D, N2 = (K2 - a12 * K1) / D)
  } else {
    E3 = c(N1 = NA_real_, N2 = NA_real_)
  }

  E3_valid = !anyNA(E3) && E3["N1"] > 0 && E3["N2"] > 0

  jacobian = function(N1, N2) {
    J11 = r1 * (1 - (2 * N1 + a21 * N2) / K1)
    J12 = -r1 * a21 * N1 / K1
    J21 = -r2 * a12 * N2 / K2
    J22 = r2 * (1 - (2 * N2 + a12 * N1) / K2)
    matrix(c(J11, J12, J21, J22), nrow = 2, byrow = TRUE)
  }

  # lable the fixed point on basis of corresponding eigenvalues.
  classify = function(eig) {
    label = ""

    if (any(abs(Im(eig)) > 1e-8)) {
      re = Re(eig)[1]
      if (re < -1e-8) {
        label = "stable spiral"
      } else if (re > 1e-8) {
        label = "unstable spiral"
      } else {
        label = "center"
      }
    } else {
      eig = Re(eig)
      if (all(eig < -1e-8)) {
        label = "stable node"
      } else if (all(eig > 1e-8)) {
        label = "unstable node"
      } else if (any(eig > 1e-8) && any(eig < -1e-8)) {
        label = "saddle (unstable)"
      } else {
        label = "non-hyperbolic"
      }
    }

    return(label)
  }


  # Collect the fixed points to analyze
  eq_names = c("E0", "E1", "E2")
  eq_points = list(E0 = E0, E1 = E1, E2 = E2)

  if (E3_valid) {
    eq_names = c(eq_names, "E3")
    eq_points$E3 = E3
  }


  stability_table = data.frame()

  n = length(eq_names)

  for (i in 1:n) {
    name = eq_names[i]
    point = eq_points[[name]]

    J = jacobian(point["N1"], point["N2"])
    eig = eigen(J, only.values = TRUE)$values # Give eigenvalues
    # No need to compute eigenvectors here.
    label = classify(eig)


    # print(paste(name, ": N1 =", point["N1"], ", N2 =", point["N2"],
    #             "| Eig =", eig[1], ",", eig[2],
    #             "| Classification =", label))

    is_stable = grepl("^stable", label)
    # Searches for a text (stable) inside label-variable
    #  and returns TRUE or FALSE

    row = data.frame(equilibrium = name,
                     N1 = point["N1"], N2 = point["N2"],
                     stable = is_stable, classification = label)

    stability_table = rbind(stability_table, row)
  }

  return(stability_table)
}


#=======================================================================
#=======================================================================

# Phase Portrait

plot_phase_portrait_LV = function(states,params,
                               t = seq(0, 50, length.out = 200)) {
  # store orbits and fixed points
  orbits = lv_orbits(states, params, t)
  stab = lv_stability(params)

  # Ecological interpretation of the fixed points.

  # if ( length(stab$stable[stab$equilibrium == "E3"]) > 0 &&
  #      stab$stable[stab$equilibrium == "E3"])
  stable_E0 = stab$stable[stab$equilibrium == "E0"]
  stable_E1 = stab$stable[stab$equilibrium == "E1"]
  stable_E2 = stab$stable[stab$equilibrium == "E2"]
  stable_E3 = stab$stable[stab$equilibrium == "E3"]

  outcome =
    if (stable_E1 && stable_E2) {
      "Bistability"

    } else if (length(stable_E3) > 0 && stable_E3) {
      "Stable coexistence"

    } else if (stable_E1) {
      "SP1 excludes SP2"

    } else if (stable_E2) {
      "SP2 excludes SP1"

    } else if (stable_E0) {
      "Exctinction"
    } else {
      "No stable equilibrium / other"
    }

  stab$label = paste0( stab$equilibrium,": ", stab$classification)



  # Set up plot boundries w.r.t.  nullclines
  N1_max = 1.3 * max(params$K1, params$K2 / params$a12)
  N2_max = 1.3 * max(params$K2, params$K1 / params$a21)

  N1_seq = seq(0, N1_max, length.out = 200) # grid on x axis

  # Store nulclines
  nc1 = data.frame(N1 = N1_seq,
                   N2 = nullcline_1(N1_seq, params$K1, params$a21))
  nc2 = data.frame(N1 = N1_seq,
                   N2 = nullcline_2(N1_seq, params$K2, params$a12))


  # Plot


  ## group = orbit: Tells R to draw each simulation trajectory as its
  # own separate independent line
  p = ggplot() +
    geom_line(data = nc1, aes(N1, N2), colour = "firebrick",
              linetype = "dashed",linewidth = 0.7) +
    geom_line(data = nc2, aes(N1, N2), colour = "darkgreen",
              linetype = "dashed", linewidth = 0.7) +
    geom_path(data = orbits, aes(N1, N2, group = orbit),
              colour = "steelblue", linewidth = 0.3) +
    geom_point(data = stab, aes(N1, N2, fill = stable),
               shape = 21, size = 2, colour = "black", stroke = 1) +
    geom_text(data = stab, aes(N1, N2, label = equilibrium),
              nudge_x = 0.3, nudge_y = 1.23, size = 2) +
    scale_fill_manual(values = c(`TRUE` = "black",
                                  `FALSE` = "white"), guide = 'none') +
    coord_cartesian(xlim = c(0, N1_max), ylim = c(0, N2_max)) +
    labs(x = expression(N[1]), y = expression(N[2]), title = outcome,
         subtitle = paste( "r1 =", params$r1, ", r2 =",params$r2,
                           "| K1 =", params$K1,", K2 =", params$K2,
                           "| a12 =", params$a12, ", a21 =", params$a21
                           )) +
    annotate("text", x = N1_max - 0.3, y = N2_max - 0.3,
             label = paste(stab$label, collapse = "\n"),
             hjust = 0.8, vjust = 0.8, size = 2) +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 8))

  return(p)
}

#=======================================================================
#=======================================================================
# for testing in readme
# Generate appropriate set of initial states for given parameter set
# n is number of initial points we want to generate.

N_init_generate = function(params, n){

  K1  = params$K1
  K2  = params$K2
  a12 = params$a12
  a21 = params$a21

  N1_max = 1.2 * max(K1, K2 / a12)
  N2_max = 1.2 * max(K2, K1 / a21)

  # Trial: initial states normally distributed around E3 on phase plane.
  D = 1 - a12 * a21

  if (abs(D) > 1e-12) {
    E3 = c(N1 = (K1 - a21 * K2) / D, N2 = (K2 - a12 * K1) / D)
  } else {
    E3 = c(N1 = NA_real_, N2 = NA_real_)
  }
  E3_valid = !anyNA(E3) && E3["N1"] > 0 && E3["N2"] > 0


  if (E3_valid) {
    N1_init = rtnorm(n = n,
                     mean =  E3[1]/2 + 1, sd = E3[1],
                     lower = 0, upper = N1_max)
    N2_init = rtnorm(n = n,
                     mean =  E3[2]/2 + 1, sd = E3[2],
                     lower = 0, upper = N2_max)
  } else{
    N1_init = c(sample(0:N1_max, size = n, replace = TRUE), 0.02)
    N2_init = c(sample(0:N2_max, size = n, replace = TRUE), 0.03)
  }
  ####

  N_init = cbind(N1_init, N2_init)
  colnames(N_init) = c("N1", "N2")

  return(N_init)
}


N_init_generate_grid = function(params, n){

  K1  = params$K1
  K2  = params$K2
  a12 = params$a12
  a21 = params$a21

  N1_max = 1.2 * max(K1, K2 / a12)
  N2_max = 1.2 * max(K2, K1 / a21)

  N1_init = seq(0.01, N1_max, length.out = n)
  N2_init = seq(0.01, N2_max, length.out = n)

  grid_df = expand.grid(N1 = N1_init, N2 = N2_init) # data frame format
  grid_matrix = as.matrix(grid_df)

  return(grid_matrix)
}

#=======================================================================
# plot_phase_portrait_LV(states = c(2, 3),
#                     params = list(r1 = 1, r2 = 0.8, K1 = 10, K2 = 8,
#                               a12 = 0.4, a21 = 0.6))
