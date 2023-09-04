# Code for Part II: Simulation

# Research Module in Econometrics & Statistics 
# Prof. Dr. Liebl & Dr. Christopher Walsh
# Winter 2021/22, M.Sc. Economics, Bonn University
# Xingyu Tao, Xuan Li, Sven Jacobs

################################################################################

#####################
#    Preparation    #
#####################

set.seed(123) # Seed for reproducibility

###################
#    Figure 14    #
###################

# Data

n <- 100
X <- runif(n)
m_X <- sin(2*pi*X)
epsilon <- rnorm(n, sd = 0.25)
Y <- m_X + epsilon
x <- seq(0, 1, length.out = n)
m_x <- sin(2*pi*x)

h <- 0.104 # Asymptotically optimal bandwidth (computed below)
    
# Estimation
    
m_hat_NW <- NW(x = x, X = X, Y = Y, K = epanechnikov, h = h)

# Plot

plot(X, Y, xlim = c(0, 0.15), ylim = c(-0.1, 1.25), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(x, m_x, lwd = 1.5, col = "dodgerblue")
lines(x, m_hat_NW, lwd = 1.5, col = "darkorange")
rect(0, -2, h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
legend("top", legend = c(paste0("Observations (n =", " ", n, ")"), "True curve", "Nadaraya-Watson"),
       pch = c(1, NA, NA), lty = c(NA, 1, 1), lwd = 1.5, col = c("black", "dodgerblue", "darkorange"),
       bty = "n", cex = 1.25)
text(0.05, -0.1, paste0("[0, h =", " ", h, ")"), cex = 1.25)
clip(0, h, -2,  2)
polygon(c(x, rev(x)), c(m_x, rev(m_hat_NW)), lty = 0, col = rgb(1, 1, 153/255, alpha = 0.4))

##############################
#    Regression functions    #
##############################

m_1 <- function(x) {sin(2*pi*x)}
m_2 <- function(x) {2 - 2*x + 2*exp(- (x - 0.5)^2 / 0.01)}
m_3 <- function(x) {2 - 2*x + 2*exp(- (x - 0)^2 / 0.01)}

#############################
#    Bandwidth functions    #
#############################

# Note:
#
# The bandwidth functions assume that the regression error is homoscedastic
# and that the design density is supported on [0, 1].

# Kernel constants for Epanechnikov kernel
R_epa <- 3/5
kappa2_epa <- 1/5

# AMISE optimal bandwidth for the local linear estimator

# Bandwidth constant
h_opt_LL_const_fun <- function(R_K, kappa2_K, sd_error, m_fun, f_fun) {
    
    # input: - R_K: integrated squared kernel (scalar)
    #        - kappa2_K: kernel variance (scalar)
    #        - sd_error: standard deviation of error (scalar)
    #        - m_fun: regression function (function)
    #        - f_fun: design density (function)
    
    # output: - h_opt_LL_const: bandwidth constant for LL (scalar)
    
    der2_m_fun <- Deriv(m_fun, nderiv = 2)
    integrand <- function(x) {der2_m_fun(x)^2 * f_fun(x)}
    
    numerator <- R_K * sd_error^2
    denominator <- kappa2_K^2 * integrate(integrand, 0, 1)$value
    
    h_opt_LL_const <- (numerator/denominator)^(1/5)
    
}

# Bandwidth function
h_opt_LL_fun <- function(const, n) {const * n^(-1/5)}

# AMISE optimal bandwidth for the Nadaraya-Watson estimator

# Bandwidth constant
h_opt_NW_const_fun <- function(R_K, kappa2_K, sd_error, m_fun, f_fun) {
    
    # input: - R_K: integrated squared kernel (scalar)
    #        - kappa2_K: kernel variance (scalar)
    #        - sd_error: standard deviation of error (scalar)
    #        - m_fun: regression function (function)
    #        - f_fun: design density (function)
    
    # output: - h_opt_NW_const: bandwidth constant for NW (scalar)
    
    der1_m_fun <- Deriv(m_fun, nderiv = 1)
    der2_m_fun <- Deriv(m_fun, nderiv = 2)
    der1_f_fun <- Deriv(f_fun, nderiv = 1)
    integrand <- function(x) {(0.5*der2_m_fun(x) + der1_f_fun(x)*der1_m_fun(x)/f_fun(x))^2 * f_fun(x)}
    
    numerator <- R_K * sd_error^2
    denominator <- kappa2_K^2 * 4 * integrate(integrand, 0, 1)$value
    
    h_opt_NW_const <- (numerator/denominator)^(1/5)
    
}

# Bandwidth function
h_opt_NW_fun <- function(const, n) {const * n^(-1/5)}

##################################
#    Main simulation function    #
##################################

# Note 1:
#
# If the regressor X is uniformly distributed (in our case on [0, 1]),
# then the bandwidths for NW and LL coincide since for uniform designs
# the asymptotic bias is the same (design bias is absent).
#
# 
# Note 2:
#
# For small sample sizes (here n = 50 and n = 100) two problems occur when
# computing the Integrated Squared Error (ISE).
# 
# Problem I:
# In a very few cases it happens that an interval (x - h, x + h) contains less
# than two observations. Then, the LL estimator at x is not defined as the 
# sum of weights (i.e. the denominator) equals 0.
#
# Problem II:
# In a few cases the computation of the integral over the boundary region fails.
# Essentially, this affects only the boundary-adjusted NW estimator.
#
# Due to the large number of Monte-Carlo repetitions (10000) we remove these
# few iterations when computing the Mean Integrated Squared Error (MISE), i.e.
# NA values are stripped by setting na.rm = TRUE.

simulation_fun <- function(sd_error, m_fun, h_opt_LL_const, h_opt_NW_const, RNG) {
    
    # input: - sd_error: standard deviation of error (scalar)
    #        - m_fun: regression function (function)
    #        - h_opt_LL_const: bandwidth constant for LL (scalar)
    #        - h_opt_NW_const: bandwidth constant for NW (scalar)
    #        - RNG: random number generator for design density (function)
    
    # output: - changes_table: summary of changes in MISE (data frame)
    
    set.seed(123) # Seed for reproducibility
    
    n <- c(50, 100, 250, 500, 1000, 2500, 5000) # Sample sizes
    
    # Pre-allocated containers
    ISE_boundary_NW <- matrix(NA, 10000, length(n))
    ISE_boundary_LL <- matrix(NA, 10000, length(n))
    ISE_boundary_NW_adjusted <- matrix(NA, 10000, length(n))
    MISE_boundary_NW <- vector("numeric", length(n))
    MISE_boundary_LL <- vector("numeric", length(n))
    MISE_boundary_NW_adjusted <- vector("numeric", length(n))
    
    # AMISE optimal bandwidths
    h_opt_LL <- h_opt_LL_fun(h_opt_LL_const, n = n)
    h_opt_NW <- if (h_opt_LL_const == h_opt_NW_const) {h_opt_LL} else {h_opt_NW_fun(h_opt_NW_const, n = n)}
    
    for (i in 1:length(n)) {
        
        B <- ifelse(n[i] < 1000, 10000, 1000) # Monte-Carlo repetitions
        
        # Functions to compute the Squared Errors (SE)
        SE_NW <- function(x) {(m_fun(x) - NW(x, X = X, Y = Y, K = epanechnikov, h = h_opt_NW[i]))^2}
        SE_LL <- function(x) {(m_fun(x) - LL(x, X = X, Y = Y, K = epanechnikov, h = h_opt_LL[i]))^2}
        SE_NW_adjusted <- function(x) {(m_fun(x) - NW_boundary(x, X = X, Y = Y, h = h_opt_NW[i],
                                                               K_interior = epanechnikov,
                                                               K_left = epanechnikov_left, K_right = epanechnikov_right,
                                                               boundary_left = 0, boundary_right = 1))^2}
        
        for (j in 1:B) {
            
            X <- RNG(n[i])
            epsilon <- rnorm(n[i], sd = sd_error)
            Y <- m_fun(X) + epsilon
        
            # Integrated Squared Errors (ISE) for left boundary region
            ISE_boundary_NW[j, i] <- tryCatch(
                expr = integrate(SE_NW, lower = 0, upper = h_opt_NW[i])$value,
                error = function(err) {NA}
            )
            ISE_boundary_LL[j, i] <- tryCatch(
                expr = integrate(SE_LL, lower = 0, upper = h_opt_NW[i])$value,
                error = function(err) {NA}
            )
            ISE_boundary_NW_adjusted[j, i] <- tryCatch(
                expr = integrate(SE_NW_adjusted, lower = 0, upper = h_opt_NW[i])$value,
                error = function(err) {NA}
            )
            
        }
        
        # Mean Integrated Squared Errors (MISE) for left boundary region
        MISE_boundary_NW[i] <- mean(ISE_boundary_NW[, i], na.rm = TRUE)
        MISE_boundary_LL[i] <- mean(ISE_boundary_LL[, i], na.rm = TRUE)
        MISE_boundary_NW_adjusted[i] <- mean(ISE_boundary_NW_adjusted[, i], na.rm = TRUE)
        
    }
    
    # Percent changes in MISE relative to standard Nadaraya-Watson estimator 
    change_LL <- round((MISE_boundary_LL - MISE_boundary_NW)/MISE_boundary_NW * 100, 2)
    change_NW_adjusted <- round((MISE_boundary_NW_adjusted - MISE_boundary_NW)/MISE_boundary_NW * 100, 2)
    
    # Boundary regions
    boundary_region <- vector("character", length(n))
    for (i in 1:length(n)) {boundary_region[i] <- paste0("[0,", " ", round(h_opt_NW[i], 3), ")")}
    
    # Data frame with results
    changes_table <- setNames(data.frame(boundary_region, change_LL, change_NW_adjusted, row.names = n),
                              c("Boundary region", "Change: LL", "Change: NW adjusted"))

}

#################
#    Designs    #
#################

# Uniform design 
f_uniform <- function(x) {1} # Support only on [0, 1]

# Clustered (non-uniform) design

# Density
f_clustered <- function(x) {3 * (5/12 - (x - 0.5)^2)} # Support only on [0, 1]

# RNG for clustered design density (via inverse transform method)
rclustered <- function(n) {
    
    # input: - n: sample size (vector)
    
    # output: - r: random numbers distributed according to clustered design density (vector)
    
    U <- runif(n)
    f <- function(x, u) {5/4*x - (x - 0.5)^3 - 1/8 - u} # First part of function is the CDF of f_clustered
    
    r <- vector("numeric", length(n))
    
    for (i in 1:n) {r[i] <- uniroot(f, c(0, 1), u = U[i])$root}
    
    r
    
}

###########################
#    Tables from paper    #
###########################

# Computation of bandwidth constants

h_opt_LL_const_table_2a <- h_opt_LL_const_fun(R_epa, kappa2_epa, 0.25, m_1, f_uniform)
h_opt_LL_const_table_2b <- h_opt_LL_const_fun(R_epa, kappa2_epa, 0.1, m_1, f_uniform)
h_opt_LL_const_table_2c <- h_opt_LL_const_fun(R_epa, kappa2_epa, 0.5, m_1, f_uniform)
h_opt_LL_const_table_2d <- h_opt_LL_const_fun(R_epa, kappa2_epa, 0.25, m_1, f_clustered)
h_opt_NW_const_table_2d <- h_opt_NW_const_fun(R_epa, kappa2_epa, 0.25, m_1, f_clustered)
h_opt_LL_const_table_3 <- h_opt_LL_const_fun(R_epa, kappa2_epa, 0.25, m_2, f_uniform)
h_opt_LL_const_table_4 <- h_opt_LL_const_fun(R_epa, kappa2_epa, 0.25, m_3, f_uniform)

# Tables

# Note:
#
# The code below runs approximately 15 hours (exact: 52780.938 sec) on a MacBook Air 2020 with M1 chip.
# One should be careful when increasing the sample size as the run time grows exponentially in the sample size.

table_2a <- simulation_fun(0.25, m_1, h_opt_LL_const_table_2a, h_opt_LL_const_table_2a, runif)
table_2b <- simulation_fun(0.1, m_1, h_opt_LL_const_table_2b, h_opt_LL_const_table_2b, runif)
table_2c <- simulation_fun(0.5, m_1, h_opt_LL_const_table_2c, h_opt_LL_const_table_2c, runif)
table_2d <- simulation_fun(0.25, m_1, h_opt_LL_const_table_2d, h_opt_NW_const_table_2d, rclustered)
table_3 <- simulation_fun(0.25, m_2, h_opt_LL_const_table_3, h_opt_LL_const_table_3, runif)
table_4 <- simulation_fun(0.25, m_3, h_opt_LL_const_table_4, h_opt_LL_const_table_4, runif)