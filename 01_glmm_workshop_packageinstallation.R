## ---------------------------
##
## Script name: 01_glmm_workshop_packageinstallation.R 
##
## Purpose of script: Code for installing packages required for GLMM workshop
##
## Author: James E Paterson
##
## Date Updated: 2024-03-21
##
## 
##
## ---------------------------

# A working internet connection is required to install new packages.
# You may be prompted to choose a 'mirror' for downloading.
# If so, try to choose the closest geographical option.

# check.packages function: install and load multiple R packages.
# Function from: https://gist.github.com/smithdanielle/9913897
# Check to see if packages are installed. Install them if they are not,
# then load them into the R session.
check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

# Make vector of needed packages
workshop.packages <- c("dplyr", "lme4", "ggplot2", 
                       "boot", "glmmTMB", "performance", 
                       "DHARMa", "lattice", "car")

# Apply function to packages
check.packages(workshop.packages)
# If prompted, easier to choose option that does NOT require compilation.