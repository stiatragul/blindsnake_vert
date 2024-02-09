# 03_visualising_data.R
# Putter Tiatragul
# Feb 2024
# Visualise relationship between number of vertebrae and other traits
# 

# libraries ---------------------------------------------------------------
# Load packages with two lines

library(ggplot2)
# library(dplyr)
library(phytools); 
library(ape); 
library(geomorph); 
# library(geiger)
# library(R.utils)
library(phylolm)
# library(factoextra)

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
vert_df <- linear_df[!is.na(linear_df$precloacal_vert) & linear_df$precloacal_vert != 0, ]

# Calculate total vertical length and vertical ratio
vert_df$total_vert <- vert_df$precloacal_vert + vert_df$postcloacal_vert
vert_df$vert_ratio <- vert_df$total_vert / vert_df$svl

plot(x = vert_df$total_length, y = vert_df$total_vert, xlab = "Total vertebrae #", ylab = "Total length (mm)", bty = "n")

# No clear pattern of vertebrae and total length, indicates there may be some species that elongate their vertebrae and others that shorten their vertebrae relative to the mean. 

# Q1 Test for sexual dimorphism ----------------------------------------------

## Question: Which species show sexual dimorphism in vert/length ratio?

# Make a subset of data where sex has been identified
dimorph_df <- vert_df[vert_df$sex %in% c('f', 'm') & 
                        !vert_df$species %in% c('Anilios aspina', 'Anilios longissimus', 'Anilios margaretae',
                                                'Anilios sp', 'Anilios zonula', 'Ramphotyphlops multilineatus', 
                                                'Sundatyphlops polygrammicus', 'Anilios fossor', 'Anilios vagurima',
                                                'Anilios chamodracaena'), ]

# Vector of species where sample size is at least 3 for males and females
include_species <- dimorph_df %>% 
  dplyr::group_by(species, sex) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= 3) %>% dplyr::distinct(species) %>% 
  dplyr::filter(species %notin% c("Anilios_systenos")) %>%
  c()

# Subset data to test for sexual dimorphism for each trait using only species that have at least 3 males and 3 females
subset_dimorph <- dimorph_df[dimorph_df$species %in% include_species$species,]

# Fit linear models and check which species have sexual dimorphism
source('code/utility/func_sexual_dimorphism_tester.R')
d_vert_ratio <- sex_dimorphic_tester("vert_ratio ~ sex + species + sex:species", subset_dimorph)
d_vert_ratio

## ANSWER: 4/29 species are sexually dimorphic. 
d_vert_ratio$dimorph_sp


# QUESTION 2 --------------------------------------------------------------

## Question: How does number of vertebrae vary with total body length by species and sex

filtered_vert_df <- vert_df %>% dplyr::group_by(species) %>% dplyr::filter(n() > 4, sex %in% c("f","m")) 

ggplot(filtered_vert_df, aes(x = total_vert, y = total_length, colour = sex)) +
  geom_point() + 
  facet_wrap(~species) + theme_bw() +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  theme(strip.text = element_text(face = "italic")) +
  ylab("Total length (mm)") +
  xlab("# of vertebrae")  

# Fit linear models where total_length ~ total vert for each species
fit_lm_2 <- by(filtered_vert_df, filtered_vert_df$species, function(df) {
  lm(total_vert ~ total_length, data = df)
})

# View the model output
lapply(fit_lm_2, summary)

# Create a more human readable output 
coefficients_df <- data.frame()

# Loop through each model output
for (species in names(fit_lm_2)) {
  # Extract coefficients and p-values for total_vert
  coef_p <- summary(fit_lm_2[[species]])$coefficients["total_length", ]
  # Extract p-value
  p_value <- coef_p["Pr(>|t|)"]
  # Extract coefficient
  coefficient <- coef_p["Estimate"]
  # Determine significance
  sig <- ifelse(p_value < 0.05, TRUE, FALSE)

  coefficients <- data.frame(species = species, coefficient, p_value, sig)
  # Bind the coefficients to the dataframe
  coefficients_df <- rbind(coefficients_df, coefficients)
}

coefficients_df

table(coefficients_df$sig)["TRUE"]

# Fit linear regression to quantify
fit_procDlm_2 <- procD.lm(total_vert ~ species + total_length + sex + species:sex, data = filtered_vert_df, turbo = TRUE)
summary(fit_procDlm_2)

## Answer:
## In most species, we see variation in length but not in # vertebrae. 

## Interpretation:
## Model shows significant effects of species and total length on number of vertebrae, 
## but no significant effects from sex and interaction between species and sex. In general, 
## longer snakes do have more number of vertebrae. 
## fit_lm_2 shows that in 7/26 species (where we have enough samples to test dimorphism) 
## total length is correlated with total vert but not all. 

# Question 3 --------------------------------------------------------------

## Question: How does vertebrae/length ratio differ across the phylogeny

filtered_vert_df$length_by_vert <- filtered_vert_df$total_length/filtered_vert_df$total_vert

vert_data <- vert_df %>% dplyr::group_by(species) %>% 
  dplyr::summarise(mean_svl = mean(svl),
                   mean_total_length = mean(total_length),
                   mean_preclo = mean(precloacal_vert),
                   mean_post = mean(postcloacal_vert),
                   mean_tot_vert = mean(precloacal_vert + postcloacal_vert),
                   body_ratio = mean(midbody_diameter/svl))
vert_df
vert_data <- as.data.frame(vert_data)
rownames(vert_data) <- vert_data$species
vert_data$ratio <- vert_data$mean_tot_vert/ vert_data$mean_total_length

check_data <- name.check(phy = sub_phy, data = vert_data)
vert_data <- vert_data[which(rownames(vert_data) %notin% check_data$data_not_tree),]

plotTree.barplot(sub_phy, setNames(vert_data$ratio, rownames(vert_data)), 
                 args.barplot=list(xlab="Vertebrae/Length ratio"))

## Answer 
## There's variation across the phylogeny. As seen in plotTree.barplot.


# Question 4 --------------------------------------------------------------

## Question: Is there  a relationship between mean total vertebrae and body shape ratio (fineness)

# Calculate width ratio
width_ratio_cex <- vert_data$body_ratio*90
width_ratio_cex <- setNames(obj = width_ratio, rownames(vert_data))

dev.off()
plot(mean_total_length ~ mean_tot_vert, cex = width_ratio_cex, data = vert_data, bty="n", pch = 19)

## Figure caption
## Size of point is body fineness. We see no clear pattern between mean total length and mean total vert


# Question 5 --------------------------------------------------------------

## Question: What is the best-fitting evolutionary model for vertebrae ratio

# Subset data to ones we have vertebrae count data
anilios_vert <- vert_data[anilios_tree$tip.label, ]
anilios_data <- anilios_data[rownames(anilios_vert), ]
anilios_data$pre_cloa <- round(anilios_vert$mean_preclo, 0)
anilios_data$tot_vert <- anilios_vert$mean_tot_vert
anilios_data$ver_rati <- anilios_vert$ratio
anilios_data$width_ratio <- anilios_data$mean_width/anilios_data$svl

# Fit and compare what evolutionary model fit best
fit.phylolm.ev <- function(trait, .data, .phy){
  formu <- paste(trait, "~", "1")
  
  fit_evol_mod_BM <- phylolm(formu, .data, .phy, measurement_error=T,model="BM")
  fit_evol_mod_OUf<- phylolm(formu, .data, .phy, measurement_error=T,model="OUfixedRoot")
  fit_evol_mod_OUr<- phylolm(formu, .data, .phy, measurement_error=T,model="OUrandomRoot")
  fit_evol_mod_LD <- phylolm(formu, .data, .phy, measurement_error=F,model="lambda")
  fit_evol_mod_EB <- phylolm(formu, .data, .phy, measurement_error=T,model="EB",lower.bound=-10,upper.bound=10)
  
  fit_ev <- list(fit_evol_mod_BM, fit_evol_mod_OUf, fit_evol_mod_OUr, fit_evol_mod_LD, fit_evol_mod_EB)
  fit.obj <- lapply(fit_ev, AIC)
  names(fit.obj) <- c("BM", "OUFixedRoot", "OURandomRoot", "Lambda", "EB")
  AIC.tab <- data.frame(AIC = unlist(fit.obj))
  return(AIC.tab)
  
}

# Find best fitting model for total vert and vert ratio
tot_vert_AIC <- fit.phylolm.ev("tot_vert", anilios_data, anilios_tree) 
ver_rati_AIC <- fit.phylolm.ev("ver_rati", anilios_data, anilios_tree) 

## Answer
## Early burst model best fit these two traits best, followed by BM and lambda

# Question 6 --------------------------------------------------------------

## Question: Does vertebrae ratio show phylogenetic signal? 

## Phylogenetic signal from fitting phylolm
sig_tot_vert <- phylolm(tot_vert~1, anilios_data, anilios_tree,model="lambda")
sig_ver_rati <- phylolm(ver_rati~1, anilios_data, anilios_tree,model="lambda")

summary(sig_ver_rati)
# anilios_data$random <- runif(nrow(anilios_data), min = 0, max = 100)
# sig_ver_rati_random <- phylolm(random~width_ratio, anilios_data, anilios_tree,model="lambda")

## Answer
## Trait is evolving under a BM model. 

## Interpretation
# Trait variation among species is entirely explained by their shared evolutionary history



# Question 7 --------------------------------------------------------------

## Question: Does vertebrae ratio correlate with body width ratio? 

# Correlation between width ratio and vert ratio --------------------------------
vert_ratio <- anilios_data$ver_rati; names(vert_ratio) <- anilios_data$species
width_ratio <- anilios_data$width_ratio; names(width_ratio) <- anilios_data$species

### Using geomorph::procD.pgls 
vert.pgls <- geomorph::procD.pgls(vert_ratio ~ width_ratio, phy = anilios_tree)
summary(vert.pgls)
predict(vert.pgls)
coefficients(vert.pgls)
dev.off()
plot(vert_ratio ~ width_ratio, bty="n", cex = width_ratio_cex[names(width_ratio)], pch = 19)
abline(a = coefficients(vert.pgls)[1], b = coefficients(vert.pgls)[2])

### using phylolm
fit_width_EB <- phylolm(ver_rati ~ width_ratio, data = anilios_data, phy = anilios_tree, measurement_error=T, model="EB", lower.bound=-10, upper.bound=10)
fit_width_BM <- phylolm(ver_rati ~ width_ratio, data = anilios_data, phy = anilios_tree, measurement_error=T, model="BM")
fit_width_LB <- phylolm(ver_rati ~ width_ratio, data = anilios_data, phy = anilios_tree, model="lambda")

summary(fit_width_EB)
summary(fit_width_BM)
summary(fit_width_LB)

## Answer: 
## There is a statistically significant correlation between vertebrae ratio and fineness of the body. 
## Species that have more vertebrae per length (vertebrae ratio > 1) are usually more fine, whereas more robust species
## have lower number of vertebrae. From phylolm, when fitted with EB or BM/Lambda, still show significant effect of width_ratio.

## Possible explanation: 
## This could be that these species need more vertebrae to fit through crevices 
## in hard soil? More robust species are found in softer soil (from blindsnakemorpho project)


# Question 8 --------------------------------------------------------------

## Question: Does vertebrae ratio correlate with any environmental factors? 
## Environmental factors: 
# mean annual temperature, 
# aridity index, 
# soil compactness.

# Total vertebrae  -- no effect
fit_phylm_total_vert <- phylolm(tot_vert ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda")
summary(fit_phylm_total_vert)

# Vertebrae ratio -- full model
fit_phylm_total_ratio <- phylolm(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda")
summary(fit_phylm_total_ratio)

## There's at least an effect of temperature, but this might be because the variables might be correlated, we can try fitting it individually.


## PICK UP HERE ### 2024-02-09

# Stepwise model selection
phylostep(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree,
             direction = "both", k = 2)



# Visualise data
dev.off()
par(mfrow=c(1,3))
plot(ver_rati ~ temp_mean, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_2, lty = 3)

plot(ver_rati ~ mean_bulk, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_3, lty = 3)

plot(ver_rati ~ ARID, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_4, lty = 3)






fit_vert_2 <- phylolm(log(tot_vert) ~ temp_mean, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_vert_3 <- phylolm(log(tot_vert) ~ mean_bulk, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_vert_4 <- phylolm(log(tot_vert) ~ ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_vert_5 <- phylolm(log(tot_vert) ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)

sum_vert_2 <- summary(fit_vert_2); sum_vert_3 <- summary(fit_vert_3); sum_vert_4 <- summary(fit_vert_4); sum_vert_5 <- summary(fit_vert_5)

sum_vert_2 
sum_vert_3
sum_vert_4
sum_vert_5

par(mfrow=c(2,2))
plot(tot_vert ~ temp_mean, data = anilios_data)
plot(tot_vert ~ mean_bulk, data = anilios_data)
plot(tot_vert ~ ARID, data = anilios_data)

# fit_vrat_2 <- phylolm(ver_rati ~ temp_mean, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# fit_vrat_3 <- phylolm(ver_rati ~ mean_bulk, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# fit_vrat_4 <- phylolm(ver_rati ~ ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# fit_vrat_5 <- phylolm(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)

summary(fit_vrat_2)
summary(fit_vrat_3)
summary(fit_vrat_4)
summary(fit_vrat_5)




