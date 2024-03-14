# 00_stat_for_paper.R
# Putter Tiatragul
# March 2024
# Statistics for methods section in the paper

# libraries ---------------------------------------------------------------
# Load packages with two lines

library(dplyr)
library(phytools); 
library(geiger);
library(ape); 
library(geomorph); 
library(phylolm)
'%notin%' <- Negate('%in%')

# source ------------------------------------------------------------------
# for prepared data from Anilios landmark and semi landmark
load('data/script_generated_data/dorsal_head_shape.rda')
load('data/script_generated_data/soil_bulk_density_data.rda')
load('data/script_generated_data/subset_mcmctree_shape.rda')

body_shape_df <- read.csv('data/script_generated_data/blindsnake_logbodyshape_ratio.csv')
PC_log_body_ratio_df <- read.csv('data/script_generated_data/blindsnake_sp_pc_logbodyrat.csv')
PC_log <- read.csv('data/script_generated_data/log_body_shape_ratio_pc_all.csv')


# Need in this script

# Full individual data
linear_df <- read.csv('data/script_generated_data/blindsnake_body_traits.csv')

# Summary data
anilios_data <- read.csv(file = 'data/script_generated_data/anilios_summary_data.csv', row.names = 1)

# Load full tree, subsetted tree, and subset 
load('data/script_generated_data/subset_mcmctree_shape.rda')
anilios_tree <- ape::drop.tip(sub_phy, tip = c("Ramphotyphlops_multilineatus","Acutotyphlops_subocularis"))
anilios_tree$tip.label

# Summarise linear measurement data ---------------------------------------

genus_species <- paste(linear_df$genus, linear_df$species, sep = "_")
linear_df$species <- genus_species

# Check vertebrae by length -----------------------------------------------

# Filter out rows with NA or 0 values in precloacal_vert
vert_df <- linear_df[!is.na(linear_df$precloacal_vert) & 
                       linear_df$precloacal_vert != 0 & 
                       linear_df$genus == "Anilios" &       # only Anilios species
                       linear_df$species != "Anilios_sp",]  # excluding unknown species

# Calculate total vertical length and vertical ratio
vert_df$total_vert <- vert_df$precloacal_vert + vert_df$postcloacal_vert
vert_df$vert_ratio <- vert_df$total_vert / vert_df$svl


# Material and methods ----------------------------------------------------

## Taxonoimic sampling

# We examined this many specimens
nrow(vert_df)

# Our data is presented by this many lineages
unique(vert_df$species)  |> length()

# Our sample size per species ranged from 
paste0("sample size ranged from ", min(table(vert_df$species)), " to ", max(table(vert_df$species)))

# Median | mean of species are:
paste("median =", median(table(vert_df$species)))
mean(table(vert_df$species))

## Institution we measured from
unique(vert_df$institution)






