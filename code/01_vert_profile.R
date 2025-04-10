## METADATA ===============================================================
## Filename: 01_vert_profile.R
## Description: Adapted from Emma Sherratt's script for Measuring and analysing vertebral length of sea snakes along the axial column
## March 2024
## Author: Sarin Tiatragul (sarin.tiatragul@anu.edu.au)
##=======================================================================##

### READ ME! ###
# Comments are instructions
# make changes where indicated

##### DO FIRST! ######
# Set working directory
# setwd(" ") # paste between the quotation the path to working directory, OR

# Data --------------------------------------------------------------------
# Read in the three sets of replicate measurements for a single specimen

# Ligatus
lig1 <- read.csv("data/x_ray/lig_R00050_ST_1.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
lig2 <- read.csv("data/x_ray/lig_R00050_ST_2.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
lig3 <- read.csv("data/x_ray/lig_R00050_ST_3.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 

# Grypus
gry1 <- read.csv("data/x_ray/gry_R146673_ST_1.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
gry2 <- read.csv("data/x_ray/gry_R146673_ST_2.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
gry3 <- read.csv("data/x_ray/gry_R146673_ST_3.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 

# Diversus
div1 <- read.csv("data/x_ray/div_R51126_ST_1.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
div2 <- read.csv("data/x_ray/div_R51126_ST_2.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 
div3 <- read.csv("data/x_ray/div_R51126_ST_3.csv", header=TRUE, row.names = 1)[,(c("X", "Y"))] 

# Plot  -------------------------------------------------------------------
source('code/utility/func_vert_profile_plotter.R')

pdf(file = "output/vertebral_profile.pdf", width = 8.3, height = 11.7)
par(mfrow=c(3,1))
# vert_raw_plotter(lig1, lig2, lig3, "A. ligatus R00050")
vert_profile_plotter(lig1, lig2, lig3, "A. ligatus R00050", .error = F)
vert_profile_plotter(div1, div2, div3, "A. diversus R51126", .error = F)
vert_profile_plotter(gry1, gry2, gry3, "A. grypus R146673", .error = F)
dev.off()
