# Packages used during the project

# Research Module in Econometrics & Statistics 
# Prof. Dr. Liebl & Dr. Christopher Walsh
# Winter 2021/22, M.Sc. Economics, Bonn University
# Xingyu Tao, Xuan Li, Sven Jacobs

################################################################################

# Packages

packages <- c(
    "Deriv",     # Computation of derivatives (symbolic)
    "haven",     # Reading Stata files
    "tidyverse", # Collection of packages for data science
    "Hmisc"      # Cutting a numeric variable into intervals
)

# Install only missing packages

installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
}

# Load packages (without showing output)

invisible(lapply(packages, library, character.only = TRUE))