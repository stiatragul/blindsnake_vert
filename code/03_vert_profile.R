# 03_vert_profile.R
# 2024-03-25
# Adapted from Emma's script for Measuring and analysing vertebral length of sea snakes along the axial column
# October 2021


### READ ME! ###
# Comments are instructions
# make changes where indicated

##### DO FIRST! ######
# Set working directory
# setwd(" ") # paste between the quotation the path to working directory, OR


# Data --------------------------------------------------------------------
# Read in the three sets of replicate measurements for a single specimen
lig1 <- read.csv("data/x_ray/lig_R00050_ST_1.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
lig2 <- read.csv("data/x_ray/lig_R00050_ST_2.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
lig3 <- read.csv("data/x_ray/lig_R00050_ST_3.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 

# Plot  -------------------------------------------------------------------
source('code/utility/func_vert_profile_plotter.R')
vert_profile_plotter(lig1, lig2, lig3, "A. ligatus R00050")


