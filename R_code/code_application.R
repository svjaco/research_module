# Code for Part III: Application

# Research Module in Econometrics & Statistics 
# Prof. Dr. Liebl & Dr. Christopher Walsh
# Winter 2021/22, M.Sc. Economics, Bonn University
# Xingyu Tao, Xuan Li, Sven Jacobs

################################################################################

##############
#    Data    #
##############

# Source: MHE Data Archive (https://economics.mit.edu/faculty/angrist/data1/mhe)

data <- "lee_2008.dta"

# Relevant variables for our analysis:
#
# difdemshare: Democratic vote share margin of victory, election t
# mdemwinnext: Probability of victory, election t+1 

##########################
#    Data preparation    #
##########################

# Read in data and keep only relevant variables
data_lee_2008 <- read_stata(data, col_select = c(use, difdemshare, mdemwinnext)) %>%
    filter(difdemshare >= -0.25 & difdemshare <= 0.25 & use == 1)

# Create bins of width = 0.005
data_lee_2008$interval <- cut2(data_lee_2008$difdemshare, seq(-0.25, 0.25, 0.005))
data_lee_2008 <- data_lee_2008[!duplicated(data_lee_2008$interval), ]
data_lee_2008 <- data_lee_2008[order(data_lee_2008$interval), ]
data_lee_2008$difdemshare_binned <- seq(-0.25, 0.245, 0.005)

difdemshare_binned <- data_lee_2008$difdemshare_binned
mdemwinnext <- data_lee_2008$mdemwinnext

###################
#    Figure 15    #
###################

# Note:
#
# The code replicates Fig. 5a of Lee (2008), without the parametric (logit) fit. 

plot(difdemshare_binned, mdemwinnext,
     xaxt = "n", xlab = "Democratic Vote Share Margin of Victory, Election t",
     yaxt = "n", ylab = "Probability of Victory, Election t+1",
     pch = 16, cex.lab = 1.25)
axis(1, at = seq(-0.25, 0.25, 0.05), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1.25)
abline(v = 0, lty = "dashed")
legend("topleft", legend = c("Local average (binsize = 0.005)"), pch = 16, bty = "n", cex = 1.25)

##################
#    Figure 5    #
##################

h_grid <- seq(0.05, 0.25, length.out = 200)

# Cross-validated bandwidths for each method and side of the threshold 0
# (including plots of CV error over bandwidth grid)

h_CV_NW_lower <- h_CV_fun(X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
                          K = epanechnikov, h_grid = h_grid, method = "NW", plot = TRUE)
h_CV_NW_upper <- h_CV_fun(X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
                          K = epanechnikov, h_grid = h_grid, method = "NW", plot = TRUE)

h_CV_LL_lower <- h_CV_fun(X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
                          K = epanechnikov, h_grid = h_grid, method = "LL", plot = TRUE)
h_CV_LL_upper <- h_CV_fun(X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
                          K = epanechnikov, h_grid = h_grid, method = "LL", plot = TRUE)

# h_CV_NW_adjusted_lower <- h_CV_fun(X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
#                           K = epanechnikov, h_grid = h_grid, method = "NW boundary",
#                           boundary_left = -0.25, boundary_right = -0.005,
#                           plot = TRUE)
# h_CV_NW_adjusted_upper <- h_CV_fun(X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
#                           K = epanechnikov, h_grid = h_grid, method = "NW boundary",
#                           boundary_left = 0, boundary_right = 0.245,
#                           plot = TRUE)

# Estimates for each method and side of the threshold 0

NW_lower <- NW(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
               K = epanechnikov, h = h_CV_NW_lower)
NW_upper <- NW(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
               K = epanechnikov, h = h_CV_NW_upper)

LL_lower <- LL(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
               K = epanechnikov, h = h_CV_LL_lower)
LL_upper <- LL(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
               K = epanechnikov, h = h_CV_LL_upper)

NW_adjusted_lower <- NW_boundary(x = difdemshare_binned[1:50], X = difdemshare_binned[1:50], Y = mdemwinnext[1:50],
                                 h = h_CV_NW_lower, K_interior = epanechnikov,
                                 K_left = epanechnikov_left, K_right = epanechnikov_right,
                                 boundary_left = -0.25, boundary_right = -0.005)
NW_adjusted_upper <- NW_boundary(x = difdemshare_binned[51:100], X = difdemshare_binned[51:100], Y = mdemwinnext[51:100],
                                 h = h_CV_NW_upper, K_interior = epanechnikov,
                                 K_left = epanechnikov_left, K_right = epanechnikov_right,
                                 boundary_left = 0, boundary_right = 0.245)

# Plot

plot(difdemshare_binned, mdemwinnext,
     xaxt = "n", xlab = "Democratic Vote Share Margin of Victory, Election t",
     yaxt = "n", ylab = "Probability of Victory, Election t+1",
     pch = 16, cex.lab = 1.25)
axis(1, at = seq(-0.25, 0.25, 0.05), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1.25)
abline(v = 0, lty = "dashed")
lines(difdemshare_binned[1:50], NW_lower, lwd = 1.5, col = "darkorange")
lines(difdemshare_binned[51:100], NW_upper, lwd = 1.5, col = "darkorange")
lines(difdemshare_binned[1:50], LL_lower, lwd = 1.5, col = "forestgreen")
lines(difdemshare_binned[51:100], LL_upper, lwd = 1.5, col = "forestgreen")
lines(difdemshare_binned[1:50], NW_adjusted_lower, lwd = 1.5, col = "chocolate4")
lines(difdemshare_binned[51:100], NW_adjusted_upper, lwd = 1.5, col = "chocolate4")
legend("topleft",
       legend = c("Local average (binsize = 0.005)", "Nadaraya-Watson", "Local linear", "NW boundary-adjusted"),
       pch = c(16, NA, NA, NA), lty = c(NA, 1, 1, 1), lwd = 1.5,
       col = c("black", "darkorange", "forestgreen", "chocolate4"),
       bty = "n", cex = 1.25)

##################
#    Figure 6    #
##################

# Confidence intervals (95% asymptotic)

confidence_interval_left <- confidence_interval_LL(X = difdemshare_binned[1:50], Y =  mdemwinnext[1:50], K = epanechnikov, h = h_CV_LL_lower)
confidence_interval_left_lower <- confidence_interval_left$confidence_interval_LL_lower
confidence_interval_left_upper <- confidence_interval_left$confidence_interval_LL_upper

confidence_interval_right <- confidence_interval_LL(X = difdemshare_binned[51:100], Y =  mdemwinnext[51:100], K = epanechnikov, h = h_CV_LL_upper)
confidence_interval_right_lower <- confidence_interval_right$confidence_interval_LL_lower
confidence_interval_right_upper <- confidence_interval_right$confidence_interval_LL_upper

# Plot

plot(difdemshare_binned, mdemwinnext,
     xaxt = "n", xlab = "Democratic Vote Share Margin of Victory, Election t",
     yaxt = "n", ylab = "Probability of Victory, Election t+1",
     pch = 16, cex.lab = 1.25)
axis(1, at = seq(-0.25, 0.25, 0.05), cex.axis = 1.25)
axis(2, at = seq(0, 1, 0.1), cex.axis = 1.25)
abline(v = 0, lty = "dashed")
lines(difdemshare_binned[1:50], LL_lower, lwd = 1.5, col = "forestgreen")
lines(difdemshare_binned[51:100], LL_upper, lwd = 1.5, col = "forestgreen")
polygon(x = c(difdemshare_binned[1:50], rev(difdemshare_binned[1:50])),
        y = c(confidence_interval_left_upper, rev(confidence_interval_left_lower)),
        col = adjustcolor("forestgreen", alpha.f = 0.1),
        border = NA)
polygon(x = c(difdemshare_binned[51:100], rev(difdemshare_binned[51:100])),
        y = c(confidence_interval_right_upper, rev(confidence_interval_right_lower)),
        col = adjustcolor("forestgreen", alpha.f = 0.1),
        border = NA)
legend("topleft",
       legend = c("Local average (binsize = 0.005)", "Local linear", "95% asymptotic confidence bands"),
       pch = c(16, NA, 15), lty = c(NA, 1, NA), lwd = 1.5,
       col = c("black", "forestgreen", adjustcolor("forestgreen", alpha.f = 0.2)),
       bty = "n", cex = 1.25)