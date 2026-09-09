            #Finding equilibria and their stability#

# For model: f(x) = dx/dt = r_m * x * (1 - x) - H(x)
# Where H(x) is one of the Holling harvester functions (harvest_fn)

# fi: Computes growth rate with i-th type holling function harvesting.
# f_derivative: Numerical derivative f'(x)

# find_equilibria: Finds all roots of f(x) = 0 in [0, l], and classifies
# each one as stable or unstable using f'(x*).

# Idea of the code:
# Instead of solving f(x)=0 algebraically, we hunt for roots numerically
# 1. Build a fine grid and evaluate f on it.
# 2. Scan for sign changes. so we know root lie between the grid points.
# 3. Use uniroot() to pin down the roots i.e. fixed points
# 4. Use f_prime to find numerical derivative of f at x*
# 5. Find the stability by LST.
# 6. Plot each equlibrium points in red for stable equilibria,
#    grey for unstable eqilibria; against the variable parameter value :
#    h_0 for type 1, and h for type 2, 3.

#=======================================================================

# Holling functionals

f0 = function(x, r_m, c){
  # x   : Population density, x = N/K_0 (can range beyond 1 if N > K_0)
  # r_m : Intrinsic growth rate
  # c  : Constant total harvest rate.

  if (x <= 0) {
    return(0)
  } # To avoid NBS.

  return(r_m * x * (1 - x) - c)
}


f1 = function(x, r_m, h_0){
  # x   : Population density, x = N/K_0 (can range beyond 1 if N > K_0)
  # r_m : Intrinsic growth rate
  # h_0  : Constant per-capita harvest rate. (same unit as r_m)
  return((r_m*x*(1 - x) - h_0*x))
}


f23 = function(y, k, h, theta) {
  # y     : rescaled density, y = N/N0 (N0 = half-saturation constant)
  # k     : rescaled carrying capacity, k = K0/N0
  # h     : rescaled harvest rate, h = h0/(r_m*N0)
  # theta : Holling exponent
  # (theta = 1 -> Type II, theta > 1 -> Type III)
  return((y * (1 - y / k)) - (h * (y^theta) / (1 + y^theta)))
}

#=======================================================================

# Numerical derivative

# params := list of parameters that one shall pass with
# corresponding params.
# eps = error
f_derivative = function(fn, x, params, eps = 0.000001){
  f_right = do.call(fn, c(list(x + eps), params))
  f_left  = do.call(fn, c(list(x - eps), params))
  return((f_right - f_left) / (2 * eps))
}

#f_derivative(f23, 0.5, list(k = 2, h = 0.5, theta = 1))

#=======================================================================

# eqb points :

find_equilibria = function(fn, params, I_min = 0, I_max = 2,
                           n_grid = 2000){

  grid = seq(I_min, I_max, length.out = n_grid)

  # Evaluate f on grid points
  f_vals = sapply(grid, function(x){
    do.call(fn, c(list(x), params))})


  roots = numeric(0)   # Collecting equilibrium points here.

  # Scan neighbouring grid points for a sign change.
  for (i in 1:(n_grid - 1)) {

    if (is.na(f_vals[i]) || is.na(f_vals[i + 1])) {
      next   # Skip if either value is undefined (e.g. numerical issue)
    }

    if (f_vals[i] == 0) {
      roots = c(roots, grid[i])
      # If exact zero landed on grid point

    } else if (f_vals[i] * f_vals[i + 1] < 0) {
      # Sign flips -> uniroot()
      root_i = uniroot(function(x) do.call(fn, c(list(x), params)),
                     c(grid[i], grid[i + 1]))$root
      roots = c(roots, root_i)
    }
  }

  if (f_vals[n_grid] == 0) {
    roots = c(roots, grid[n_grid])
  } # Because we missed last point in the grid.

  roots = unique(round(roots, 6))
  # Drop duplicate roots up to 6th decimal place.


  stable = sapply(roots, function(x_star){
    f_derivative(fn, x_star, params) < 0
  }) # Store root as stable if f'(x) < 0 (LST).


  return(data.frame(x_star = roots, stable = stable ))
}


# find_equilibria(
#   f23,
#   list(k = 0.1, h =0.23, theta = 1),
#   I_min = 0,
#   I_max = 1
# )

#=======================================================================

plot_bifurcation = function(fn, var_param, range_param = seq(0, 2,
                            by = 0.0001), params, I_min = 0, I_max = 2){

  # Lists to store different equilibrium
  stable_eqb = list()
  unstable_eqb = list()

  for (i in range_param) {

    # Set parameter value
    p = params
    p[[var_param]] = i

    # Find equilibria
    eqb = find_equilibria(fn, p, I_min = I_min, I_max = I_max)
    #print(eqb)
    #print(eqb$stable)

    # Classify eqb values in stable, unstable points.
    stable = eqb$x_star[eqb$stable == TRUE]
    unstable = eqb$x_star[eqb$stable == FALSE]

    stable_eqb[[length(stable_eqb) + 1]] = c(i, stable)
    unstable_eqb[[length(unstable_eqb) + 1]] = c(i, unstable)
  }

  ##############################################################

  # Critical Harvest rate h_c :
  # Population will definately go exctinct for harvest rate > h_c.
  # I have calculate this algebrically and put here directly to plot it
  # on the graph.
  if (identical(fn, f0)) {
    h_c = params$r_m / 4

  } else if (identical(fn, f1)) {
    h_c = params$r_m

  } else if (identical(fn, f23) && params$theta == 1) {
    h_c = ((params$k + 1)^2) / (4 * params$k)

  } else if (identical(fn, f23) && params$theta > 1) {
    h_c = NULL
  }

  ################################################################

  # Base plot
  plot(NA, xlab = var_param, ylab = "Equilibrium",
       xlim = range(range_param), ylim = c(I_min, I_max),
       main = "Bifurcation diagram", cex.main = 0.7)


  # Plot stable eqb
  for (p in stable_eqb) {
    points(rep(p[1], length(p)-1), p[-1], col = "red", pch = 16,
           cex = 0.3)
  }


  # Plot unstable eqb
  for (p in unstable_eqb) {
    points(rep(p[1], length(p)-1), p[-1], col = "gray", pch = 16,
           cex = 0.3)
  }


  if (!is.null(h_c) && h_c >= min(range_param) &&
      h_c <= max(range_param)) {

    lines(rep(h_c, 2), c(I_min, I_max), col = "green",lwd = 1,
          cex = 0.5)
  }

  legend("topright", legend = c("Stable eqb", "Unstable eqb",
        "Critical harvest"), col = c("red", "gray", "green"),lty = 1,
        lwd = 2, cex = 0.5)
}



# plot_bifurcation(f23, var_param = "h",
#                  range_param = seq(0, 2.5, by = 0.001),
#                  params = list(k = 1.8, theta = 1), I_min = 0,
#                  I_max = 2.5)

# plot_bifurcation(f0, var_param = "c",
#                  range_param = seq(0, 3, by = 0.01),
#                  params = list( r_m = 2), I_min = 0,
#                  I_max = 3)

#=======================================================================
