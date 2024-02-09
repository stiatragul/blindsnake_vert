# 04_model_fitting.R
# Sarin Tiatragul 
# July 2022 modified Feb 2023
# Combine phenotypic data with environment data and plot

# libraries ---------------------------------------------------------------
# Load packages with two lines
pkgs = c("phytools", "ape", "geomorph", "ggplot2", "dplyr", "corHMM", "geiger") # package names
inst = lapply(pkgs, library, character.only = TRUE) # load them

library(R.utils)
library(phylolm)
library(factoextra)

'%notin%' <- Negate('%in%')

# source ------------------------------------------------------------------
# for prepared data from Anilios landmark and semi landmark
load('data/script_generated_data/dorsal_head_shape.rda')
load('data/script_generated_data/soil_bulk_density_data.rda')
load('data/script_generated_data/subset_mcmctree_shape.rda')

linear_df <- read.csv('data/script_generated_data/blindsnake_body_traits.csv')
body_shape_df <- read.csv('data/script_generated_data/blindsnake_logbodyshape_ratio.csv')
PC_log_body_ratio_df <- read.csv('data/script_generated_data/blindsnake_sp_pc_logbodyrat.csv')
PC_log <- read.csv('data/script_generated_data/log_body_shape_ratio_pc_all.csv')

# Load full tree, subsetted tree, and subset 
load('data/script_generated_data/subset_mcmctree_shape.rda')


anilios_tree <- ape::drop.tip(sub_phy, tip = c("Ramphotyphlops_multilineatus","Acutotyphlops_subocularis"))

anilios_data <- read.csv(file = 'data/script_generated_data/anilios_summary_data.csv', row.names = 1)
anilios_tree$tip.label

# Summarise linear measurement data ---------------------------------------

genus_species <- paste(linear_df$genus, linear_df$species, sep = "_")
linear_df$species <- genus_species

# Check vertebrae by length -----------------------------------------------

vert_df <- linear_df %>% 
  dplyr::filter(!is.na(precloacal_vert)) %>% 
  dplyr::filter(precloacal_vert != 0) %>% 
  dplyr::mutate(total_vert = precloacal_vert + postcloacal_vert,
                vert_ratio = total_vert/svl)

plot(x = vert_df$total_length, y = vert_df$total_vert, xlab = "Total vertebrae #", ylab = "SVL", bty = "n")

# No clear pattern of vertebrae and total length, indicates there may be some species that elongate their vertebrae and others that shorten their vertebrae relative to the mean. 

# Test for sexual dimorphism ----------------------------------------------

# Make a subset of data where sex has been identified
dimorph_df <- vert_df %>%
  dplyr::filter(sex %in% c('f', 'm')) %>% 
  dplyr::filter(!species %in% c('Anilios aspina', 'Anilios longissimus', 'Anilios margaretae',
                                'Anilios sp', 'Anilios zonula', 'Ramphotyphlops multilineatus', 
                                'Sundatyphlops polygrammicus', 'Anilios fossor', 'Anilios vagurima',
                                'Anilios chamodracaena'))

# Vector of species where sample size is at least 3 for males and females
include_species <- dimorph_df %>% dplyr::group_by(species, sex) %>% dplyr::summarise(n = n()) %>% 
  dplyr::filter(n >= 3) %>% dplyr::distinct(species) %>% 
  dplyr::filter(species %notin% c("Anilios_systenos")) %>%
  c()

# Subset data to test for sexual dimorphism for each trait using only species that have at least 3 males and 3 females
subset_dimorph <- dimorph_df[dimorph_df$species %in% include_species$species,]

# Fit linear models and check which species have sexual dimorphism
source('code/utility/func_sexual_dimorphism_tester.R')
d_vert_ratio <- sex_dimorphic_tester("vert_ratio ~ sex + species + sex:species", subset_dimorph)

d_vert_ratio

## Four species are sexually dimorphic. 

# How does this vary across the phylogeny?
filtered_vert_df <- vert_df %>% dplyr::group_by(species) %>% dplyr::filter(n() > 4, sex %in% c("f","m")) 

ggplot(filtered_vert_df, aes(x = total_vert, y = total_length, colour = sex)) +
  geom_point() + facet_wrap(~species) + theme_bw() +
  theme(strip.text = element_text(face = "italic")) +
  xlab("Total length (mm)") +
  ylab("# of vertebrae")  



library(ggtree)
filtered_vert_df$length_by_vert <- filtered_vert_df$total_length/ filtered_vert_df$total_vert

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
vert_data$ratio <- vert_data$mean_tot_vert/ vert_data$mean_total_length

check_data <- name.check(phy = sub_phy, data = vert_data)
vert_data <- vert_data[which(rownames(vert_data) %notin% check_data$data_not_tree),]

plotTree.barplot(sub_phy, setNames(vert_data$ratio, rownames(vert_data)), 
                 args.barplot=list(xlab="Vertebrae/Length ratio"))

### Visualise

length_vert_plot <- ggplot(data = vert_data, aes(x = mean_tot_vert, y = mean_total_length, size = body_ratio)) + 
  geom_point() +
  theme_classic() 

length_vert_plot


library(mvMORPH)
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


tot_vert_AIC <- fit.phylolm.ev("tot_vert", anilios_data, anilios_tree) 
ver_rati_AIC <- fit.phylolm.ev("ver_rati", anilios_data, anilios_tree) 

# Phylogenetic signal -----------------------------------------------------
sig_tot_vert <- phylolm(tot_vert~1, anilios_data, anilios_tree,model="lambda")
# sig_ver_rati <- phylolm(ver_rati~1, anilios_data, anilios_tree,model="lambda")
sig_ver_rati <- phylolm(ver_rati~width_ratio, anilios_data, anilios_tree,model="lambda")

# Correlation between width ratio and vert ratio --------------------------------
vert_ratio <- anilios_data$ver_rati; names(vert_ratio) <- anilios_data$species
width_ratio <- anilios_data$width_ratio; names(width_ratio) <- anilios_data$species


### geomorph::procD.pgls 
vert.pgls <- procD.pgls(vert_ratio ~ (width_ratio), phy = anilios_tree)
summary(vert.pgls)
predict(vert.pgls)
coefficients(vert.pgls)
dev.off()
plot(vert_ratio ~ width_ratio, bty="n")


### using phylolm
fit_width_EB <- phylolm(ver_rati ~ width_ratio, data = anilios_data, phy = anilios_tree, measurement_error=T, model="EB",lower.bound=-10,upper.bound=10)
fit_width_BM <- phylolm(ver_rati ~ width_ratio, data = anilios_data, phy = anilios_tree, measurement_error=T, model="BM")

summary(fit_width_BM)
summary(fit_width_EB)



# Does number of vertebrae correlate with anything? 

# Env predictors including soil compactness, percentage sand, mean temp, and mean prec

# Using total number of vertebrae
fit_total_vert <- phylolm(tot_vert ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda")



summary(fit_total_vert)




fit_total_ratio <- phylolm(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda")
summary(fit_total_ratio)

# Stepwise model


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

dev.off()
par(mfrow=c(1,3))
plot(ver_rati ~ temp_mean, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_2, lty = 3)

plot(ver_rati ~ mean_bulk, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_3, lty = 3)

plot(ver_rati ~ ARID, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_4, lty = 3)


