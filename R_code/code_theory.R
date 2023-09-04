# Code for Part I: Theory

# Research Module in Econometrics & Statistics 
# Prof. Dr. Liebl & Dr. Christopher Walsh
# Winter 2021/22, M.Sc. Economics, Bonn University
# Xingyu Tao, Xuan Li, Sven Jacobs

################################################################################

#####################
#    Preparation    #
#####################

set.seed(123) # Seed for reproducibility

##################
#    Figure 1    #
##################

# Data for example

m_fun <- function(x) {sin(2*pi*x)} # True regression function
n <- 100
X <- seq(0, 1, length.out = n) # For illustration fixed design. Similar to runif(n).
m_X <- m_fun(X)
epsilon <- rnorm(n, sd = 0.25)
Y <- m_X + epsilon

h <- 0.2

# Estimation

m_hat_NW <- NW(x = X, X = X, Y = Y, K = epanechnikov, h = h)
m_hat_NW_0.5 <- NW(x = 0.5, X = X, Y = Y, K = epanechnikov, h = h) # Interior
m_hat_NW_0 <- NW(x = 0, X = X, Y = Y, K = epanechnikov, h = h) # Boundary

# Figure 1a

plot(X, Y, ylab = "", xaxt = "n", cex.lab = 1.25, cex.axis = 1.25)
axis(1, at = c(0, 0.2, 0.5 - h, 0.4, 0.5, 0.6, 0.5 + h, 0.8, 1), cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW, lwd = 1.5, col = "darkorange")
points(0.5, m_hat_NW_0.5, pch = 16, col = "red")
segments(0.5 - h, m_hat_NW_0.5, 0.5 + h, m_hat_NW_0.5, col = "red")
rect(0.5 - h, -2, 0.5 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
curve(3/4 * (1 - ((0.5 - x)/h)^2), from = 0.5 - h, to = 0.5 + h, lty = "dashed", add = TRUE)
legend("topright", legend = c("True curve", "Nadaraya-Watson", "Kernel weights", "Estimate at 0.5"),
       pch = c(NA, NA, NA, 16), lty = c(1, 1, 2, NA), lwd = 1.5,
       col = c("dodgerblue", "darkorange", "black", "red"), bty = "n", cex = 1.25)
text(0.5, -1.5, paste0("h = ", h), cex = 1.25)

# Figure 1b

plot(X, Y, ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW, lwd = 1.5, col = "darkorange")
points(0, m_hat_NW_0, pch = 16, col = "red")
segments(0 - h, m_hat_NW_0, 0 + h, m_hat_NW_0, col = "red")
rect(0 - h, -2, 0 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
curve(3/4 * (1 - ((0 - x)/h)^2), from = 0, to = 0 + h, lty = "dashed", add = TRUE)
legend("topright", legend = c("True curve", "Nadaraya-Watson", "Kernel weights", "Estimate at 0"),
       pch = c(NA, NA, NA, 16), lty = c(1, 1, 2, NA), lwd = 1.5,
       col = c("dodgerblue", "darkorange", "black", "red"), bty = "n", cex = 1.25)
text(0.5, -1.5, paste0("h = ", h), cex = 1.25)

##################
#    Figure 8    #
##################

# Data

n_example_II <- 100
X_example_II <- sort(runif(n_example_II))
m_X_example_II <- -1 + 2*X_example_II
epsilon_example_II <- 0
Y_example_II <- m_X_example_II + epsilon_example_II

h_example_II <- 0.2

# Estimation

m_hat_NW_example_II <- NW(x = X_example_II, X = X_example_II, Y = Y_example_II,
                          K = epanechnikov, h = h_example_II)

# Plot

plot(X_example_II, Y_example_II, xlab = "X", ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X_example_II, m_X_example_II, lwd = 1.5, col = "dodgerblue")
lines(X_example_II, m_hat_NW_example_II, lwd = 1.5, col = "darkorange")
legend("top", legend = c("True curve", "Nadaraya-Watson"),
       lty = c(1, 1), lwd = 1.5, col = c("dodgerblue", "darkorange"), bty = "n", cex = 1.25)
text(0.5, -0.95, paste0("h = ", h_example_II), cex = 1.25)

##################
#    Figure 9    #
##################

# Estimation

m_hat_LL <- LL(x = X, X = X, Y = Y, K = epanechnikov, h = h)

# Plot

plot(X, Y, ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW, lwd = 1.5, col = "darkorange")
lines(X, m_hat_LL, lwd = 1.5, col = "forestgreen")
legend("topright", legend = c("True curve", "Nadaraya-Watson", "Local linear"),
       lty = c(1, 1, 1), lwd = 1.5, col = c("dodgerblue", "darkorange", "forestgreen"),
       bty = "n", cex = 1.25)
text(0.5, -1.5, paste0("h = ", h), cex = 1.25)

##################
#    Figure 2    #
##################

# Estimation

output_LL_0 <- LL_2(x = 0, X = X, Y = Y, K = epanechnikov, h = h)
m_hat_LL_0 <- output_LL_0$estimates
slope_m_hat_LL_0 <- output_LL_0$slopes

# Plot

plot(X, Y, xlim = c(-0.2, 0.25), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_LL, lwd = 1.5, col = "forestgreen")
points(0, m_hat_LL_0, pch = 16, col = "red")
rect(0 - h, -2, 0 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(v = 0, lty = "dashed")
legend("bottomright", legend = c("True curve", "Local linear", "Estimate at 0"),
       pch = c(NA, NA, 16), lty = c(1, 1, NA), lwd = 1.5, col = c("dodgerblue", "forestgreen", "red"),
       bty = "n", cex = 1.25)
text(-0.1, -1.5, paste0("h = ", h), cex = 1.25)
clip(0 - h, 0 + h, -2,  2)
abline(m_hat_LL_0, slope_m_hat_LL_0, col = "red")

##################
#    Figure 3    #
##################

# Estimation

output_NW_0.5 <- NW_2(x = 0.5, X = X, Y = Y, K = epanechnikov, h = h)
output_LL_0.5 <- LL_2(x = 0.5, X = X, Y = Y, K = epanechnikov, h = h)
m_hat_LL_0.5 <- output_LL_0.5$estimates
effective_kernel_NW_0.5 <- output_NW_0.5$effective_kernels
effective_kernel_LL_0.5 <- output_LL_0.5$effective_kernels

output_NW_0 <- NW_2(x = 0, X = X, Y = Y, K = epanechnikov, h = h)
effective_kernel_NW_0 <- output_NW_0$effective_kernels
effective_kernel_LL_0 <- output_LL_0$effective_kernels

# Figure 3a

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(mat = layout_matrix, heights = c(0.95, 0.05))

plot(X, Y, xlim = c(0.25, 0.75), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW, lwd = 1.5, col = "darkorange")
points(0.5, m_hat_NW_0.5, pch = 16, col = "red")
points(X, effective_kernel_NW_0.5, pch = 18)
rect(0.5 - h, -2, 0.5 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(h = 0, lty = "dashed")

plot(X, Y, xlim = c(0.25, 0.75), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_LL, lwd = 1.5, col = "forestgreen")
points(0.5, m_hat_LL_0.5, pch = 16, col = "red")
points(X, effective_kernel_LL_0.5, pch = 18)
rect(0.5 - h, -2, 0.5 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(h = 0, lty = "dashed")

par(mar = c(0, 0, 0, 0))

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

legend("top", legend = c("True curve", "NW", "LL", "Estimate at 0.5", "Effective kernel at 0.5"),
       pch = c(NA, NA, NA, 16, 18), lty = c(1, 1, 1, NA, NA), lwd = 1.5,
       col = c("dodgerblue", "darkorange", "forestgreen", "red", "black"),
       bty = "n", horiz = TRUE, cex = 1.25)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Restore default

# Figure 3b

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(mat = layout_matrix, heights = c(0.95, 0.05))

plot(X, Y, xlim = c(0, 0.25), ylim = c(-0.1, 1.25), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW, lwd = 1.5, col = "darkorange")
points(0, m_hat_NW_0, pch = 16, col = "red")
points(X, effective_kernel_NW_0, pch = 18)
rect(0 - h, -2, 0 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(h = 0, lty = "dashed")

plot(X, Y, xlim = c(0, 0.25), ylim = c(-0.1, 1.25), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_LL, lwd = 1.5, col = "forestgreen")
points(0, m_hat_LL_0, pch = 16, col = "red")
points(X, effective_kernel_LL_0, pch = 18)
rect(0 - h, -2, 0 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(h = 0, lty = "dashed")

par(mar = c(0, 0, 0, 0))

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

legend("top", legend = c("True curve", "NW", "LL", "Estimate at 0", "Effective kernel at 0"),
       pch = c(NA, NA, NA, 16, 18), lty = c(1, 1, 1, NA, NA), lwd = 1.5,
       col = c("dodgerblue", "darkorange", "forestgreen", "red", "black"),
       bty = "n", horiz = TRUE, cex = 1.25)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Restore default

###################
#    Figure 10    #
###################

# Estimation

m_hat_NW_boundary <- NW_boundary(x = X, X = X, Y = Y, h = h,
                                 K_interior = epanechnikov, K_left = epanechnikov_left, K_right = epanechnikov_right,
                                 boundary_left = 0, boundary_right = 1)

# Plot

plot(X, Y, ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW, lwd = 1.5, col = "darkorange")
lines(X, m_hat_LL, lwd = 1.5, col = "forestgreen")
lines(X, m_hat_NW_boundary, lwd = 1.5, col = "chocolate4")
legend("topright", legend = c("True curve", "Nadaraya-Watson", "Local linear", "NW boundary-adjusted"),
       lty = c(1, 1, 1, 1), lwd = 1.5, col = c("dodgerblue", "darkorange", "forestgreen", "chocolate4"),
       bty = "n", cex = 1.25)
text(0.5, -1.5, paste0("h = ", h), cex = 1.25)

###################
#    Figure 11    #
###################

# Estimation

output_NW_boundary_0 <- NW_boundary_2(x = 0, X = X, Y = Y, h = h,
                                      K_interior = epanechnikov, K_left = epanechnikov_left, K_right = epanechnikov_right,
                                      boundary_left = 0, boundary_right = 1)
m_hat_NW_boundary_0 <- output_NW_boundary_0$estimates
effective_kernel_NW_boundary_0 <- output_NW_boundary_0$effective_kernels

# Plot

layout_matrix <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)
layout(mat = layout_matrix, heights = c(0.95, 0.05))

plot(X, Y, xlim = c(0, 0.25), ylim = c(-0.1, 1.25), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_NW_boundary, lwd = 1.5, col = "chocolate4")
points(0, m_hat_NW_boundary_0, pch = 16, col = "red")
points(X, effective_kernel_NW_boundary_0, pch = 18)
rect(0 - h, -2, 0 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(h = 0, lty = "dashed")

plot(X, Y, xlim = c(0, 0.25), ylim = c(-0.1, 1.25), ylab = "", cex.lab = 1.25, cex.axis = 1.25)
lines(X, m_X, lwd = 1.5, col = "dodgerblue")
lines(X, m_hat_LL, lwd = 1.5, col = "forestgreen")
points(0, m_hat_LL_0, pch = 16, col = "red")
points(X, effective_kernel_LL_0, pch = 18)
rect(0 - h, -2, 0 + h, 2, border = NA, col = rgb(1, 1, 153/255, alpha = 0.2))
abline(h = 0, lty = "dashed")

par(mar = c(0, 0, 0, 0))

plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")

legend("top", legend = c("True curve", "NW boundary-adjusted", "LL", "Estimate at 0", "Effective kernel at 0"),
       pch = c(NA, NA, NA, 16, 18), lty = c(1, 1, 1, NA, NA), lwd = 1.5,
       col = c("dodgerblue", "chocolate4", "forestgreen", "red", "black"),
       bty = "n", horiz = TRUE, cex = 1.25)

par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1)) # Restore default