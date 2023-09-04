# Functions used during the project

# Research Module in Econometrics & Statistics 
# Prof. Dr. Liebl & Dr. Christopher Walsh
# Winter 2021/22, M.Sc. Economics, Bonn University
# Xingyu Tao, Xuan Li, Sven Jacobs

################################################################################

################################################################################
#                                 Preparation                                  #
################################################################################

rm(list = ls()) # Remove all objects from current workspace

################################################################################
#                                 Functions                                    #
################################################################################

#############################
#    Epanechnikov kernel    #
#############################

epanechnikov <- function(u) {ifelse(abs(u) <= 1, 3/4 * (1 - u^2), 0)}

# Popular alternative choice: Gaussian kernel
gaussian <- dnorm

###################################
#    Nadaraya-Watson estimator    #
###################################

NW <- function(x, X, Y, K, h) {
    
    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h: bandwidth (scalar)
    
    # output: - estimates: estimates for the regression function at the evaluation points (vector)

    absolute_weights <- rbind(sapply(X, function(X_i) {K((X_i - x)/h)}))
    
    effective_kernels <- absolute_weights/rowSums(absolute_weights)
    
    estimates <- as.vector(effective_kernels %*% Y)
    
}

##########################################################
#    Nadaraya-Watson estimator with additional output    #
##########################################################

NW_2 <- function(x, X, Y, K, h) {
    
    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h: bandwidth (scalar)
    
    # output: - list containing
    #               - estimates: estimates for the regression function at the evaluation points (vector)
    #               - effective_kernels: effective kernels at the evaluation points (matrix)
    
    absolute_weights <- rbind(sapply(X, function(X_i) {K((X_i - x)/h)}))
    
    effective_kernels <- absolute_weights/rowSums(absolute_weights)
    
    estimates <- as.vector(effective_kernels %*% Y)
    
    list(estimates = estimates, effective_kernels = effective_kernels)
    
}

################################
#    Local linear estimator    #
################################

LL <- function(x, X, Y, K, h) {
    
    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h: bandwidth (scalar) 
    
    # output: - estimates for the regression function at the evaluation points (vector)

    sapply(x, function(x_i) {
        
        X_matrix <- cbind(1, X - x_i)
        W_matrix <- diag(K((X - x_i)/h))
        
        estimate <- t(c(1, 0)) %*% solve((t(X_matrix)%*%W_matrix%*%X_matrix)) %*% t(X_matrix) %*% W_matrix %*% Y
        
    })

}
    
#######################################################
#    Local linear estimator with additional output    #
#######################################################

LL_2 <- function(x, X, Y, K, h) {
    
    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h: bandwidth (scalar)

    # output: - list containing
    #               - estimates: estimates for the regression function at the evaluation points (vector)
    #               - effective_kernels: effective kernels at the evaluation points (matrix)
    #               - slopes: estimates for the slope of the regression function at the evaluation points (vector)
    
    list_complete <- sapply(x, function(x_i) {
        
        X_matrix <- cbind(1, X - x_i)
        W_matrix <- diag(K((X - x_i)/h))
        
        auxiliary_matrix <- solve((t(X_matrix)%*%W_matrix%*%X_matrix)) %*% t(X_matrix) %*% W_matrix
        
        effective_kernel <- t(c(1, 0)) %*% auxiliary_matrix
        slope <- t(c(0, 1)) %*% auxiliary_matrix %*% Y
        estimate <- effective_kernel %*% Y 
        
        list(estimate, effective_kernel, slope)
        
    })
    
    estimates <- unlist(list_complete[1, ])
    effective_kernels <- matrix(unlist(list_complete[2, ]), nrow = length(x), ncol = length(Y), byrow = TRUE) 
    slopes <- unlist(list_complete[3, ])
    
    list(estimates = estimates, effective_kernels = effective_kernels, slopes = slopes)
        
}

#######################################
#    Epanechnikov boundary kernels    #
#######################################

# Note:
# 
# The following boundary kernels rely on Müller (1991).
# However, the left boundary kernels of Müller are right boundary kernels for us, and vice versa.
# This is because of the different definition of the kernel estimator concerning the argument of the kernel.
# Müller uses K(x - X), whereas we define K(X - x).

# Right boundary kernels (Müller, 1991, Table 1)

epanechnikov_right <- function(u, rho) {
    
    ifelse(u >= -1 & u <= rho, 6*(1+u)*(rho-u)*(1/(1+rho)^3)*(1 + 5*((1-rho)/(1+rho))^2 + 10*((1-rho)/(1+rho)^2)*u), 0)
    
}

# Left boundary kernels

epanechnikov_left <- function(u, rho) {epanechnikov_right(-u, rho)}

#####################################################
#    Boundary-adjusted Nadaraya-Watson estimator    #
#####################################################

NW_boundary <- function(x, X, Y, h, K_interior, K_left, K_right, boundary_left, boundary_right) {
    
    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - h: bandwidth (scalar)
    #        - K_interior: kernel used in the interior (function)
    #        - K_left: left boundary kernels (function)
    #        - K_right: right boundary kernels (function)
    #        - boundary_left: lower boundary of the support of X (scalar)
    #        - boundary_right: upper boundary of the support of X (scalar)
    
    # output: - estimates: estimates for the regression function at the evaluation points (vector)
    
    absolute_weights <- matrix(NA, length(x), length(Y))
    
    for (i in 1:length(x)) {
        
        if (x[i] >= boundary_left & x[i] < boundary_left + h) {
            
            absolute_weights[i, ] <- K_left(u = (X - x[i])/h, rho = (x[i] - boundary_left)/h)
            
        } else if (x[i] >= boundary_left + h & x[i] <= boundary_right - h) {
            
            absolute_weights[i, ] <- K_interior((X - x[i])/h)
            
        } else if (x[i] > boundary_right - h & x[i] <= boundary_right) {
            
            absolute_weights[i, ] <- K_right(u = (X - x[i])/h, rho = (boundary_right - x[i])/h)
            
        }
    
    }
    
    effective_kernels <- absolute_weights/rowSums(absolute_weights)
    
    estimates <- as.vector(effective_kernels %*% Y)
    
}

############################################################################
#    Boundary-adjusted Nadaraya-Watson estimator with additional output    #
############################################################################

NW_boundary_2 <- function(x, X, Y, h, K_interior, K_left, K_right, boundary_left, boundary_right) {
    
    # input: - x: evaluation points (vector)
    #        - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - h: bandwidth (scalar)
    #        - K_interior: kernel used in the interior (function)
    #        - K_left: left boundary kernels (function)
    #        - K_right: right boundary kernels (function)
    #        - boundary_left: lower boundary of the support of X (scalar)
    #        - boundary_right: upper boundary of the support of X (scalar)
    
    # output: - list containing
    #               - estimates: estimates for the regression function at the evaluation points (vector)
    #               - effective_kernels: effective kernels at the evaluation points (matrix)
    
    absolute_weights <- matrix(NA, length(x), length(Y))
    
    for (i in 1:length(x)) {
        
        if (x[i] >= boundary_left & x[i] < boundary_left + h) {
            
            absolute_weights[i, ] <- K_left(u = (X - x[i])/h, rho = (x[i] - boundary_left)/h)
            
        } else if (x[i] >= boundary_left + h & x[i] <= boundary_right - h) {
            
            absolute_weights[i, ] <- K_interior((X - x[i])/h)
            
        } else if (x[i] > boundary_right - h & x[i] <= boundary_right) {
            
            absolute_weights[i, ] <- K_right(u = (X - x[i])/h, rho = (boundary_right - x[i])/h)
            
        }
        
    }
    
    effective_kernels <- absolute_weights/rowSums(absolute_weights)
    
    estimates <- as.vector(effective_kernels %*% Y)
    
    list(estimates = estimates, effective_kernels = effective_kernels)
    
}

#####################
#    LOOCV error    #
#####################

CV_error_fun <- function(X, Y, K, h, method,
                         K_left = epanechnikov_left, K_right = epanechnikov_right,
                         boundary_left = 0, boundary_right = 1) {
    
    # input: - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h: bandwidth (scalar)
    #        - method: "NW" for Nadaraya-Watson, "LL" for local linear regression,
    #                  "NW boundary" for boundary-adjusted Nadaraya-Watson (character)
    #        - K_left: left boundary kernels (function), with default
    #        - K_right: right boundary kernels (function), with default
    #        - boundary_left: lower boundary of the support of X (scalar), with default
    #        - boundary_right: upper boundary of the support of X (scalar), with default
    
    # output: - CV_error: LOOCV error for bandwidth h (scalar)
    
    if (method == "NW") {
        
        NW_output <- NW_2(x = X, X = X, Y = Y, K = K, h = h)
        estimates <- NW_output$estimates
        effective_kernels <- NW_output$effective_kernels
    
    } else if (method == "LL") {
        
        LL_output <- LL_2(x = X, X = X, Y = Y, K = K, h = h)
        estimates <- LL_output$estimates
        effective_kernels <- LL_output$effective_kernels
        
    } else if (method == "NW boundary") {
        
        NW_boundary_output <- NW_boundary_2(x = X, X = X, Y = Y, h = h,
                                            K_interior = K, K_left = K_left, K_right = K_right,
                                            boundary_left = boundary_left, boundary_right = boundary_right)
        estimates <- NW_boundary_output$estimates
        effective_kernels <- NW_boundary_output$effective_kernels
        
    }
    
    CV_error <- 1/length(Y) * sum(((Y - estimates)/(1 - diag(effective_kernels)))^2)
    
}

##############################
#    CV optimal bandwidth    #
##############################

h_CV_fun <- function(X, Y, K, h_grid, method,
                     K_left = epanechnikov_left, K_right = epanechnikov_right,
                     boundary_left = 0, boundary_right = 1,
                     plot = FALSE) {
    
    # input: - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h_grid: bandwidth grid to find minimum of CV error (vector)
    #        - method: "NW" for Nadaraya-Watson, "LL" for local linear regression,
    #                  "NW boundary" for boundary-adjusted Nadaraya-Watson (character)
    #        - K_left: left boundary kernels (function), with default
    #        - K_right: right boundary kernels (function), with default
    #        - boundary_left: lower boundary of the support of X (scalar), with default
    #        - boundary_right: upper boundary of the support of X (scalar), with default
    #        - plot: plot CV error over bandwidth grid (boolean), with default
    
    # output: - h_CV: CV optimal bandwidth (scalar)
    #         - plot of CV errors if plot == TRUE (plot)
    
    CV_errors <- sapply(h_grid, function(h) {
        
        CV_error_fun(X = X, Y = Y, K = K, h = h, method = method,
                     K_left = K_left, K_right = K_right,
                     boundary_left = boundary_left, boundary_right = boundary_right)
        
    })
    
    h_CV <- h_grid[which.min(CV_errors)]
    
    if (h_CV == h_grid[1] | h_CV == h_grid[length(h_grid)]) {
        
        warning("Selected bandwidth equals the smallest or largest grid value.
                You may wish to expand the grid after consulting the plot of the CV error.")
        
    }
    
    if (plot == TRUE) {
        
        plot(h_grid, CV_errors,
             xlab = "h", ylab = "CV error", cex.lab = 1.25, cex.axis = 1.25, cex = 0.75)
        rug(h_grid, ticksize = 0.015)
        abline(v = h_CV, col = "red")
        legend("top", legend = bquote(paste("h"["CV"]*" = ", .(round(h_CV, 3)))), bty = "n", cex = 1.25)
        
    }
    
    return(h_CV)
    
}

####################################################################
#    95% asymptotic confidence bands for local linear estimator    #
####################################################################

# Note:
#
# Confidence bands for the Nadaraya-Watson estimator can be computed the same way.

confidence_interval_LL <- function(X, Y, K, h) {
    
    # input: - X: data for the regressor (vector)
    #        - Y: data for the regressand (vector)
    #        - K: kernel (function)
    #        - h: bandwidth (scalar)
    
    # output: - list containing
    #               - confidence_interval_LL_upper: upper points for the 95% asymptotic confidence interval (vector)
    #               - confidence_interval_LL_lower: lower points for the 95% asymptotic confidence interval (vector)
    
    LL_output <- LL_2(x = X, X = X, Y = Y, K = K, h = h)
    estimates <- LL_output$estimates
    effective_kernels <- LL_output$effective_kernels
    
    squared_prediction_errors <- ((Y - estimates)/(1 - diag(effective_kernels)))^2
    error_matrix <- diag(squared_prediction_errors)
    
    estimates_variance <- sapply(X, function(x_i) {
        
        X_matrix <- cbind(1, X - x_i)
        W_matrix <- diag(K((X - x_i)/h))
        
        auxiliary_matrix <- solve((t(X_matrix)%*%W_matrix%*%X_matrix))
        estimate_covariance_matrix <- auxiliary_matrix %*% t(X_matrix)%*%W_matrix%*%error_matrix%*%W_matrix%*%X_matrix %*% auxiliary_matrix
        estimate_variance <- estimate_covariance_matrix[1, 1]
        
    })
    
    confidence_interval_LL_upper <- estimates + qnorm(p = 1 - 0.05/2)*sqrt(estimates_variance)
    confidence_interval_LL_lower <- estimates - qnorm(p = 1 - 0.05/2)*sqrt(estimates_variance)
    
    list(confidence_interval_LL_upper = confidence_interval_LL_upper, confidence_interval_LL_lower = confidence_interval_LL_lower)
    
}