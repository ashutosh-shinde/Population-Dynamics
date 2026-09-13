library(deSolve)
library(ggplot2)
library(ggrepel)
library(msm)
#=======================================================================
# Preay-predator models

models_list = list(

  LV = function(t, states, params){
    with(as.list(c(states, params)), {
      dN = r*N - b*N*P
      dP = g*b*N*P - d*P
      return(list(c(N = dN, P = dP)))
    })
  },

  Vt = function(t, states, params){
    with(as.list(c(states, params)), {
      dN = r*N*(1 - N/K) - b*N*P
      dP = g*b*N*P - d*P
      return(list(c(N = dN, P = dP)))
    })
  },

  RM = function(t, states, params){
    with(as.list(c(states, params)), {
      dN = r*N*(1 - N/K) - b*N*P/(N0 + N)
      dP = g*b*N*P/(N0 + N) - d*P
      return(list(c(N = dN, P = dP)))
    })
  }
)


#=======================================================================
#=======================================================================
#Orbit


pp_orbit = function(model, states, params, t){

  ode_rhs = if (is.function(model)){
    model
  } else{
    stop("Unknown model", model)
  }


  if (is.null(dim(states))){
    states = data.frame(N = states[1], P = states[2])
  }else{
    states = as.data.frame(states)
  }


  orbit_structure = data.frame()
  n = nrow(states)

  for (i in 1:n){
    current_state = c(N = states[i, "N"], P = states[i, "P"])
    Traj = as.data.frame(ode(y = current_state, times = t,
                             func = ode_rhs, parms = params,
                             method = "radau"))

    Traj$orbit = i
    orbit_structure = rbind(orbit_structure, Traj)
  }
  return(orbit_structure)
}


#=======================================================================
#=======================================================================
# Equlibria

equilibria_LV = function(params){
  with(params, {
    eq = data.frame(
      name = c('e0', 'e1'),
      N = c(0, d/(g*b)),
      P = c(0, r/b))
    return(eq)
  })
}


equilibria_Vt = function(params){
  with(params, {
    N_star = d/(g*b)
    P_star = r/b * (1 - N_star/K)

    eq = data.frame(
      name = c("e0", "e1", "e2"),
      N = c(0, K, N_star),
      P = c(0, 0, P_star))
    eq$positive = eq$N >= 0 & eq$P >= 0
    eq = eq[eq$positive, ]
    return(eq)
  })
}


equilibria_RM = function(params){
  with(params, {
    N_star = d*N0/(g*b - d)
    P_star = (r/b)*(1 - N_star/K)*(N0 + N_star)

    eq = data.frame(
      name =c("e0", "e1", "e2"),
      N = c(0, K, N_star),
      P = c(0, 0, P_star))
    eq$positive = eq$N >= 0 & eq$P >= 0
    eq = eq[eq$positive, ]
    return(eq)
  })
}


#=======================================================================
# Jacobian

jacobian_LV = function(N, P, params){
  with(params, {
    J11 = r - b*P
    J12 = -b*N
    J21 = g*b*P
    J22 = g*b*N - d
    J = matrix(c(J11, J12,J21, J22), nrow = 2, byrow = TRUE)
    return(J)
  })
}


jacobian_Vt = function(N, P, params){
  with(params, {
    J11 = r*(1 - (2*N)/K) - b*P
    J12 = -b*N
    J21 = g*b*P
    J22 = g*b*N - d
    J = matrix(c(J11, J12, J21, J22), nrow = 2, byrow = TRUE)
    return(J)
  })
}


jacobian_RM = function(N, P, params){
  with(params, {
    J11 = r*(1 - (2*N)/K) - b*P*N0/(N0 + N)^2
    J12 = -b*N/(N0 + N)
    J21 = g*b*P*N0/(N0 + N)^2
    J22 = g*b*N/(N0 + N) - d
    J = matrix(c(J11, J12, J21, J22), nrow = 2, byrow = TRUE)
    return(J)
  })
}

#=======================================================================
#Nullclines

nullclines_LV = function(params, N_max, P_max, n = 200){
  with(params, {
    #nullcline of dN/dt = 0.
    # df stores (N, P, type) for constant line P = r/b.
    nl_N = data.frame(N = seq(0, N_max, length.out = n),
                        P = r/b,
                        type = "N nullcline", segment = '1')

    # Nullcline of dP/dt = 0
    nl_P = data.frame(N = d/(g*b),
                      P = seq(0, P_max, length.out = n),
                      type = "P nullcline", segment = '1')

    nl_N_0 = data.frame(N = 0, P = seq(0, P_max, length.out = n),
                        type = "N nullcline", segment = '0')
    nl_P_0 = data.frame(N = seq(0, N_max, length.out = n), P = 0,
                        type = "P nullcline", segment = '0')

    return(rbind(nl_N_0, nl_N, nl_P_0, nl_P))
  })
}

nullclines_Vt = function(params, N_max, P_max, n = 200){
  with(params, {
    Nseq = seq(0, N_max, length.out = n)
    nl_N = data.frame(N = Nseq,
                      P = r/b * (1 - Nseq/K),
                      type = "N nullcline", segment = '1')
    nl_N = nl_N[nl_N$P >= 0, ]   # drop NBS part of P =.... .

    nl_P = data.frame(N = d/(g*b),
                      P = seq(0, P_max, length.out = n),
                      type = "P nullcline", segment = '1')

    nl_N_0 = data.frame(N = 0, P = seq(0, P_max, length.out = n),
                        type = "N nullcline", segment = '0')
    nl_P_0 = data.frame(N = seq(0, N_max, length.out = n), P = 0,
                        type = "P nullcline", segment = '0')

    return(rbind(nl_N_0, nl_N, nl_P_0, nl_P))
  })
}

nullclines_RM = function(params, N_max, P_max, n = 200){
  with(params, {
    Nseq = seq(0, N_max, length.out = n)
    nl_N = data.frame(N = Nseq,
                      P = r/b * (1 - Nseq/K) * (N0 + Nseq),
                      type = "N nullcline", segment = '1')
    nl_N = nl_N[nl_N$P >= 0, ]

    nl_P = data.frame(N = d*N0/(g*b - d),
                      P = seq(0, P_max, length.out = n),
                      type = "P nullcline", segment = '1')

    nl_N_0 = data.frame(N = 0, P = seq(0, P_max, length.out = n),
                        type = "N nullcline", segment = '0')
    nl_P_0 = data.frame(N = seq(0, N_max, length.out = n), P = 0,
                        type = "P nullcline", segment = '0')

    return(rbind(nl_N_0, nl_N, nl_P_0, nl_P))
  })
}

#=======================================================================
# call-functions

get_equilibria = function(model, params){
  if (identical(model, models_list$LV)) {
    return(equilibria_LV(params))

  } else if (identical(model, models_list$Vt)) {
    return(equilibria_Vt(params))

  } else if (identical(model, models_list$RM)) {
    return(equilibria_RM(params))

  } else {
    stop("Unknown model")
  }
}


get_jacobian = function(model, N, P, params){
  if (identical(model, models_list$LV)) {
    return(jacobian_LV(N, P, params))

  } else if (identical(model, models_list$Vt)) {
    return(jacobian_Vt(N, P, params))

  } else if (identical(model, models_list$RM)) {
    return(jacobian_RM(N, P, params))
  } else {
    stop("Unknown model")
  }
}

get_nullclines = function(model, params, N_max, P_max){
  if (identical(model, models_list$LV)) {
    return(nullclines_LV(params, N_max, P_max))
  } else if (identical(model, models_list$Vt)) {
    return(nullclines_Vt(params, N_max, P_max))
  } else if (identical(model, models_list$RM)) {
    return(nullclines_RM(params, N_max, P_max))
  } else {
    stop("Unknown model")
  }
}


#=======================================================================
#=======================================================================
# Stability


pp_stability = function(model, params){

  classify = function(eig){

    tol = 1e-8

    if (any(abs(Im(eig)) > tol)) {
      re = Re(eig[1])
      if (re < -tol) {
        return('stable spiral')
      } else if (re > tol) {
        return('unstable spiral')
      } else {
        return("center (or other non-hyperbolic)")
      }
    } else {
      eig = Re(eig)
      if (all(eig < -tol)) {
        return('stable node')
      } else if (all(eig > tol)) {
        return("unstable node")
      } else if (any(eig < -tol) && any(eig > tol)) {
        return("saddle point")
      } else{
        return("(non-hyperbolic)")
      }
    }
  }


  eq = get_equilibria(model = model, params = params)

  n = nrow(eq)

  stability_table = data.frame(equilibrium = eq$name,
                               N = eq$N, P = eq$P,
                               eigenvalue_1 = NA, eigenvalue_2 = NA,
                               stable = FALSE, classification = "")


  for (i in 1:n){
    N = eq$N[i]
    P = eq$P[i]

    J = get_jacobian(model = model,N = N, P = P, params = params)
    eig = eigen(J, only.values = TRUE)$values
    label = classify(eig)

    stability_table$eigenvalue_1[i] = eig[1]
    stability_table$eigenvalue_2[i] = eig[2]
    stability_table$stable[i] = grepl("^stable", label)
    stability_table$classification[i] = label
  }

  return(stability_table)
}


#=======================================================================
#=======================================================================
#Plot

plot_phase_portrait_PP = function(model, model_name = NULL,
                                  states, params,
                                  t = seq(0, 50, length.out = 200)){

  orbits = pp_orbit(model = model , states = states,
                    params = params, t = t)
  stab = pp_stability(model, params)

  N_max = max(orbits$N, stab$N, na.rm = TRUE) * 1.25
  P_max = max(orbits$P, stab$P, na.rm = TRUE) * 1.25

  nullcline_df = get_nullclines(model, params, N_max, P_max)



  ##eq points, eigen values, calssification, params to print on graph
  fmt_eig = function(z) {
    re = Re(z)
    im = Im(z)
    if (abs(im) < 1e-8) {
      sprintf("%.2f", re)
    } else {
      sprintf("%.2f%s%.2fi", re, ifelse(im >= 0, "+", "-"), abs(im))
    }
  }

  stab_lines = sprintf(
    "%s: (N=%.2f, P=%.2f)  eig=(%s, %s)  [%s]",
    stab$equilibrium, stab$N, stab$P,
    sapply(stab$eigenvalue_1, fmt_eig),
    sapply(stab$eigenvalue_2, fmt_eig),
    stab$classification
  )
  stab_text = paste(stab_lines, collapse = "\n")

  param_text = paste(names(params), unlist(params),
                    sep = " = ", collapse = "\n ")
  ###

  p = ggplot() +
    geom_line(data = nullcline_df,
              aes(N, P, group = interaction(type, segment),
                  colour = type),
              linetype = "dashed", linewidth = 0.6) +
    scale_colour_manual(values = c("N nullcline" = "darkgreen",
                                   "P nullcline" = "firebrick"),
                        name = NULL) +
    geom_path(data = orbits, aes(N, P, group = orbit),
              colour = 'steelblue', linewidth = 0.3) +
    geom_point(data = stab, aes(N, P, fill = stable),
               shape = 21, size = 2, colour = "black", stroke = 1) +
    geom_text_repel(data = stab, aes(N, P, label = equilibrium),
                    size = 2, box.padding = 0.4, point.padding = 0.2,
                    min.segment.length = 0) +
    scale_fill_manual(values = c(`TRUE` = "black",
                                 `FALSE` = "white"), guide = 'none') +
    coord_cartesian(xlim = c(0, N_max) , ylim = c(0, P_max)) +
    labs(x = "Prey population (N)",y = "Predator population (P)",
         title = model_name,
         subtitle = stab_text) +
    annotate("text", x = N_max, y = P_max, label = param_text,
             hjust = 1, vjust = 1, size = 2.8, family = "mono") +
    theme_minimal() +
    theme(plot.subtitle = element_text(size = 6))

  return(p)
}


#=======================================================================
#=======================================================================
#test


# params_1 = list(r = 2, b = 2, g = 0.8, d = 0.5)
# plot_phase_portrait_PP(model = models_list$LV, model_name = 'LV',
#                        states = c(2, 3),
#                        params = params_1)

# p = plot_phase_portrait_PP(models_list$LV, model_name = 'LV',
#                            states = matrix(c(0.5, 0.5,
#                                              1.5, 1.5,
#                                              2.0, 0.5,
#                                              0.3, 2.0),
#                                            ncol = 2,
#                                            byrow = TRUE,
#                                            dimnames =
#                                              list(NULL, c("N", "P"))),
#                            params = list(r = 1, b = 1, g = 1, d = 1),
#                            t = seq(0, 50,  length.out = 500))
# print(p)
#=======================================================================
#=======================================================================
#=======================================================================
