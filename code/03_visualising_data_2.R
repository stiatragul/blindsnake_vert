## METADATA ===============================================================
## Filename: 03_visualising_data_2.R
## Description: Visualise relationship between number of vertebrae and other traits. And test for correlations. THIS ONE WITH VERTEBRA RATIO calculated diffferently
## 
## R version: 4.4.0 for Windows
## Author: Putter Tiatragul (sarin.tiatragul@anu.edu.au)
##=======================================================================##

# libraries ---------------------------------------------------------------
library(ggplot2); library(dplyr); library(phytools); 
library(geiger); library(ape); library(geomorph); 
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

# Full individual data
linear_df <- read.csv('data/script_generated_data/blindsnake_body_traits.csv')
linear_df <- linear_df[linear_df$reg_no %notin% c("R.8515","R.127956"),] ### Dropped these samples because mis-id as troglodytes but in Queensland. Probably torresianus or broomi

# Summary data
anilios_data <- read.csv(file = 'data/script_generated_data/anilios_summary_data.csv', row.names = 1)

# Load full tree, subsetted tree, and subset 
load('data/script_generated_data/subset_mcmctree_shape.rda')
anilios_tree <- ape::drop.tip(sub_phy, tip = c("Ramphotyphlops_multilineatus","Acutotyphlops_subocularis"))
anilios_tree$tip.label

b_tree <- sub_phy
b_tree$tip.label <- gsub(pattern = "Anilios_", replacement = "A. ", b_tree$tip.label)
b_tree$tip.label <- gsub(pattern = "Acutotyphlops_", replacement = "Acu. ", b_tree$tip.label)
b_tree$tip.label <- gsub(pattern = "Ramphotyphlops_", replacement = "R. ", b_tree$tip.label)

plotTree(b_tree, ftype="i", lwd = 2, fsize=1, type = "fan")

# Summarise linear measurement data ---------------------------------------

genus_species <- paste(linear_df$genus, linear_df$species, sep = "_")
linear_df$species <- genus_species


# Check vertebrae by length -----------------------------------------------

# Filter out rows with NA or 0 values in precloacal_vert
vert_df <- linear_df[!is.na(linear_df$precloacal_vert) & linear_df$precloacal_vert != 0 & linear_df$reg_no != "R152714", ]

# Calculate total vertebral length and vertebral ratio
vert_df$total_vert <- vert_df$precloacal_vert + vert_df$postcloacal_vert
vert_df$vert_ratio <- vert_df$total_vert / (vert_df$total_length - vert_df$headlength_xray)

plot(x = vert_df$total_length, y = vert_df$total_vert, xlab = "Total vertebrae #", ylab = "Total length (mm)", bty = "n")

# No clear pattern of vertebrae and total length, indicates there may be some species that elongate their vertebrae and others that shorten their vertebrae relative to the mean. 

# Total number of vertebrae -----------------------------------------------

source('code/utility/func_plotTreeboxplot.R')

# Total number of vertebrae
total_vertebrae <- setNames(vert_df$total_vert, nm = vert_df$species)
total_vertebrae <- total_vertebrae[which(names(total_vertebrae) %in% sub_phy$tip.label)]
geiger::name.check(phy = sub_phy, data = total_vertebrae)
spp<-factor(names(total_vertebrae),untangle(ladderize(sub_phy),"read.tree")$tip.label)

summary(vert_df$total_vert)
vert_df[order(vert_df$total_vert), c("species", "total_vert")]

## Summary SE for species
Vertebra_summary_table <- vert_df %>% 
  dplyr::group_by(species) %>% 
  summarise(mean_vert = mean(total_vert),
            max_vert = max(total_vert),
            min_vert = min(total_vert),
            median_vert = median(total_vert),
            mean_preclo = mean(precloacal_vert),
            max_preclo = max(precloacal_vert),
            min_preclo = min(precloacal_vert),
            median_preclo = median(precloacal_vert),
  )

write.csv(Vertebra_summary_table, file = "assets/tables/vertebra_summary.csv", row.names = FALSE)

Vertebra_summary_table_by_sex <- vert_df %>% 
  dplyr::group_by(species, sex) %>% 
  dplyr::filter(sex %in% c("f","m")) %>% 
  summarise(sample_size = n(),
            mean_vert = mean(total_vert),
            max_vert = max(total_vert),
            min_vert = min(total_vert),
            median_vert = median(total_vert),
            mean_preclo = mean(precloacal_vert),
            max_preclo = max(precloacal_vert),
            min_preclo = min(precloacal_vert),
            median_preclo = median(precloacal_vert),
            mean_postclo = mean(postcloacal_vert),
            max_postclo = max(postcloacal_vert),
            min_postclo = min(postcloacal_vert),
            median_postclo = median(postcloacal_vert),
  )

Vertebra_summary_table_by_sex

## Head length x-ray
skull_summary <- vert_df %>% 
  dplyr::group_by(species, sex) %>% 
  dplyr::filter(sex %in% c("f","m")) %>% 
  summarise(sample_size = n(),
            skull_length = mean(headlength_xray),
            rel_skull = mean(headlength_xray/svl))

write.csv(skull_summary, file = "assets/tables/skull_summary.csv", row.names = FALSE)


# Q1 Test for sexual dimorphism ----------------------------------------------
## Question: Which species show sexual dimorphism in vert/length ratio?

# Make a subset of data where sex has been identified
dimorph_df <- vert_df[vert_df$sex %in% c('f', 'm') & 
                        !vert_df$species %in% c('Anilios aspina', 'Anilios longissimus', 'Anilios margaretae',
                                                'Anilios sp', 'Anilios zonula', 'Ramphotyphlops multilineatus', 
                                                'Sundatyphlops polygrammicus', 'Anilios fossor', 'Anilios vagurima',
                                                'Anilios chamodracaena'), ]

# Vector of species where sample size is at least 3 for males and females
exclude_species <- dimorph_df %>% 
  dplyr::group_by(species, sex) %>% 
  dplyr::summarise(n = dplyr::n()) %>% 
  dplyr::filter(n <= 3) %>% dplyr::distinct(species) %>% 
  dplyr::filter(species %notin% c("Anilios_systenos")) %>%
  c()

# Subset data to test for sexual dimorphism for each trait using only species that have at least 3 males and 3 females
subset_dimorph <- dimorph_df[dimorph_df$species %notin% exclude_species$species,]

# Fit linear models and check which species have sexual dimorphism
source('code/utility/func_sexual_dimorphism_tester.R')
d_vert_ratio <- sex_dimorphic_tester("vert_ratio ~ sex + species + sex:species", subset_dimorph)
d_vert_ratio

## ANSWER: 3/14 species are sexually dimorphic in vertebra ratio.
d_vert_ratio$dimorph_sp

postcloacal_dimorph <- sex_dimorphic_tester("postcloacal_vert ~ sex + species + sex:species", subset_dimorph)
postcloacal_dimorph$dimorph_sp
postcloacal_dimorph$table

precloacal_dimorph <- sex_dimorphic_tester("precloacal_vert ~ sex + species + sex:species", subset_dimorph)
precloacal_dimorph$dimorph_sp
precloacal_dimorph$table


# QUESTION 2 --------------------------------------------------------------

## Question: How does number of vertebrae vary with total body length by species and sex

filtered_vert_df <- vert_df %>% dplyr::group_by(species) %>% dplyr::filter(n() > 4, sex %in% c("f","m"), reg_no != "R152714") 

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
                   mean_vert_ratio = mean(vert_ratio, na.rm = T), ### Modified to use number of (tbl - head length) / number of vertebrae
                   body_ratio = mean(midbody_diameter/total_length),
                   # aspect = max(total_length/midbody_diameter),
                   aspect = mean(total_length/midbody_diameter),
                   max_total_vert = max(total_vert),
                   max_preclo = max(precloacal_vert),
                   max_svl = max(svl),
                   max_tbl = max(total_length),
                   max_vert_ratio = max(vert_ratio, na.rm = T),)

vert_data <- as.data.frame(vert_data)
rownames(vert_data) <- vert_data$species



### Another way is to use the data from the biggest total length individual in that species
slice_df <- vert_df %>% dplyr::group_by(species) %>% 
  dplyr::slice(which.max(total_length)) %>% 
  dplyr::ungroup()

slice_df$ratio <- slice_df$total_vert /  (slice_df$total_length - slice_df$headlength_xray)  ## MODIFIED TO USE TOTAL LENGTH - HEAD LENGTH as denominator
slice_df$aspect <- slice_df$total_length/slice_df$midbody_diameter

# Calculating vertebrae number to skeleton_length ratio
vert_data$ratio <- vert_data$mean_vert_ratio

check_data <- geiger::name.check(phy = sub_phy, data = vert_data)
vert_data_full <- vert_data
vert_data <- vert_data[which(rownames(vert_data) %notin% check_data$data_not_tree),]

plotTree.barplot(sub_phy, setNames(vert_data$ratio, rownames(vert_data)), 
                 args.barplot=list(xlab="Number of vertebrae/neck-to-tail length (1/mm)"))


# Aspect ratio as boxplot
aspect_v <- setNames(vert_df$total_length / vert_df$midbody_diameter, nm = vert_df$species)
aspect_v <- aspect_v[which(names(aspect_v) %in% sub_phy$tip.label)]

geiger::name.check(phy = sub_phy, data = aspect_v)
spp<-factor(names(aspect_v),untangle(ladderize(sub_phy),"read.tree")$tip.label)

## Statistics about aspect ratio
summary(aspect_v)
aspect_v[(which.min(aspect_v))]
aspect_v[(which.max(aspect_v))]


table_s3 <- vert_df %>% 
  dplyr::group_by(species) %>% 
  dplyr::summarise(mean_aspect = mean(total_length/midbody_diameter),
                   max_aspect = max(total_length/midbody_diameter),
                   min_aspect = min(total_length/midbody_diameter),
                   sample_n = n()) 

write.csv(table_s3, file = "assets/tables/aspect_ratios_summary.csv", row.names = F)

## This is using old code where we can't adjust the ylim

# Figure 1: phylo + boxplots 

# pdf(file = "output/tree_boxplot_aspectRatio.pdf", width = 20, height = 12)
plotTree.boxplot(sub_phy,x=aspect_v~spp,
                 args.boxplot = list(xlab="Aspect ratio (total length/body width)",
                                     # ylim=c(0,200),
                                     col="grey99"))
# dev.off()

# Fig vertebrae number boxplot
vert_num <- setNames(vert_df$total_vert, nm = vert_df$species)
vert_num <- vert_num[which(names(vert_num) %in% sub_phy$tip.label)]

pdf(file = "output/tree_boxplot_totalVert.pdf", width = 10, height = 12)
plotTree.boxplot(sub_phy,x=vert_num~spp,
                 args.boxplot = list(xlab="Total number of vertebrae",
                                     # ylim=c(20,170),
                                     col="grey99"))
dev.off()


## Answer 
## There's variation across the phylogeny. As seen in plotTree.barplot.


# Question 4 --------------------------------------------------------------

## Question: Is there  a relationship between mean total vertebrae and body shape ratio (fineness)

# Calculate width ratio
width_ratio_cex <- vert_data$body_ratio*90
width_ratio_cex <- setNames(obj = width_ratio_cex, rownames(vert_data))

dev.off()
plot(mean_total_length ~ mean_tot_vert, cex = width_ratio_cex, data = vert_data, bty="n", pch = 19,
     xlab = "Mean total vertebrae", ylab = "Mean total length (mm)")
# text(y = vert_data$mean_total_length, x = vert_data$mean_tot_vert, labels = names(width_ratio_cex), pos = 4, adj = 0.5)

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
## Early burst model best fit for length/vertebrae ratio but closely followed by BM (not 2 AIC score different) 


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

## Correlation between total vertebrae number positively correlates with body size

anilios_data$max_vert <- anilios_vert$max_total_vert
anilios_data$mean_vert <- anilios_vert$mean_tot_vert
anilios_data$max_svl <- anilios_vert$max_svl

mean_vert <- anilios_data$mean_vert; names(mean_vert) <- anilios_data$species 
max_vert <- anilios_data$max_vert; names(max_vert) <- anilios_data$species 
max_svl <- anilios_data$max_svl; names(max_svl) <- anilios_data$species
max_tbl <- setNames(anilios_vert$max_tbl, anilios_vert$species)

colnames(anilios_data)
colnames(vert_data_full)

### Data for species not in the phylogeny
not_tree_data <- vert_data_full[which(rownames(vert_data_full) %notin% c(anilios_tree$tip.label, "Anilios_sp")),]
not_tree_labs <- gsub(pattern = "Anilios_", replacement = "", x = not_tree_data$species)

## Test pleomerism for species in tree
pleomerism.pgls <- geomorph::procD.pgls(max_tbl ~ max_vert, phy = anilios_tree)
summary(pleomerism.pgls)
coefficients(pleomerism.pgls)
sp_labs <- gsub(pattern = "Anilios_", replacement = "", x = anilios_data$species)

## Physignal
# max_tbl.physig <- physignal(A = max_tbl, phy = anilios_tree, iter = 999)
# max_tbl.physig
max_tbl.physigz <- physignal.z(A = max_tbl, phy = anilios_tree, iter = 999, lambda = "mean")
summary(max_tbl.physigz)

# max_vert.physig <- physignal(A = max_vert, phy = anilios_tree, iter = 999)
max_vert.physigz <- physignal.z(A = max_vert, phy = anilios_tree, iter = 999, lambda = "mean")
summary(max_vert.physigz)

fit_pleomerism <- phylolm(max_tbl ~ max_vert, phy = anilios_tree, model="lambda")
summary(fit_pleomerism)

# Using phylosig
phylosig_maxtbl<-phytools::phylosig(tree = anilios_tree, x = max_tbl, method = "lambda", test=TRUE)
# phylosig_maxtbl<-phytools::phylosig(tree = anilios_tree, x = log(max_tbl), method = "lambda", test=TRUE)
phylosig_maxtbl

phylosig_maxvert<-phytools::phylosig(tree = anilios_tree, x = max_vert, method = "lambda", test=TRUE)
phylosig_maxvert

### PLOTTING
# Log max body length against log max number of vertebrae



### FIGURE 2: PGLS tbl v vert number
### Emma suggest plotting as phylomorphospace and add species that are not on the tree overlaid

pdf(file = "output/phylomorpho-tbl-vert.pdf", width = 10, height = 7.5)

phylomorphospace(anilios_tree, anilios_vert[,c("max_total_vert", "max_tbl")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Max. total body length (mm)", xlab = "Maximum number of total vertebrae")
# Add slope from pgls (but this is just for species in tree)
abline(a = coefficients(pleomerism.pgls)[1], b = coefficients(pleomerism.pgls)[2])
summary(pleomerism.pgls)[1]

# add species not in tree to overlay text next to it
points(y = not_tree_data$max_tbl, x = not_tree_data$max_total_vert, cex = 2, col="blue", pch = 20)
text(y = not_tree_data$max_tbl, x = not_tree_data$max_total_vert, labels = not_tree_labs, pos = 1, cex = 0.7, col="blue")

dev.off()

## Including Head & Polly 2007 typhlopid data
# typh_data <- read.csv("data/HeadPollyData/typhlopid_data.csv")
# typh_data
# 
# plot(max_tbl~max_vert, bty="n", pch = 19,
#      ylab = "Maximum SVL", xlab = "Maximum total vertebrae", 
#      xlim=c(100,500), ylim=c(100,800))
# points(x = typh_data$max_vert, y = typh_data$TBL, col = "red", pch = 19)
# 
# text(y = max_tbl, x = max_vert, labels = sp_labs, pos = 1, cex = 0.7)
# text(y = typh_data$TBL, x = typh_data$max_vert, labels = typh_data$Taxon, pos = 1, cex = 0.7)

# Question 8 --------------------------------------------------------------

## Question: Does vertebrae ratio correlate with body width ratio? 

# Correlation between width ratio and vert ratio --------------------------------
vert_ratio <- anilios_vert$ratio; names(vert_ratio) <- anilios_data$species
width_ratio <- anilios_vert$body_ratio; names(width_ratio) <- anilios_data$species
aspect_ratio <- anilios_vert$aspect; names(aspect_ratio) <- anilios_data$species

### Using geomorph::procD.pgls 
vert.pgls <- geomorph::procD.pgls(vert_ratio ~ aspect_ratio, phy = anilios_tree)
summary(vert.pgls)
coefficients(vert.pgls)

dev.off()
# Plotting subset just Anilios species
hist(aspect_ratio)
hist(vert_ratio)

# Vertebrae ratio against aspect ratio

plot(vert_ratio ~ aspect_ratio, bty="n", 
     # cex = width_ratio_cex[names(width_ratio)], 
     pch = 19, xlab = "Mean aspect ratio", ylab = "Mean vertebra ratio (1/mm)")
text(x = aspect_ratio, y = vert_ratio, labels = sp_labs, pos = 1, cex = 0.7)
abline(a = coefficients(vert.pgls)[1], b = coefficients(vert.pgls)[2])


# Mean aspect ratio (y) against max. number of vertebrae
ratio_maxvert.pgls <- geomorph::procD.pgls(max_vert ~ aspect_ratio, phy = anilios_tree)
summary(ratio_maxvert.pgls)
coefficients(ratio_maxvert.pgls)

plot(max_vert ~ aspect_ratio, bty="n", pch = 19, 
     ylab = "Max number of vertebrae", xlab = "Mean aspect ratio")
text(y = max_vert, x = aspect_ratio, labels = sp_labs, pos = 1, cex = 0.7)
abline(a = coefficients(ratio_maxvert.pgls)[1], b = coefficients(ratio_maxvert.pgls)[2])

# Mean aspect ratio (y) against mean number of vertebrae
ratio_meanvert.pgls <- geomorph::procD.pgls(mean_vert ~ aspect_ratio, phy = anilios_tree)
summary(ratio_meanvert.pgls)
coefficients(ratio_meanvert.pgls)

### Figure 3: COMBINED Phylomorphosapce 
pdf(file = "output/phylomorpho-vert-aspect-total2.pdf", width = 9, height = 12)
par(mfrow=c(2,1))

## Vertebrae ratio v mean aspect ratio
phylomorphospace(anilios_tree, anilios_vert[,c("aspect", "ratio")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Mean vertebrae ratio (1/mm)", xlab = "Mean aspect ratio",
                 xlim=c(0,200))
# Add slope from pgls (but this is just for species in tree)
abline(a = coefficients(vert.pgls)[1], b = coefficients(vert.pgls)[2])
# add species not in tree to overlay text next to it
points(y = not_tree_data$ratio, x = not_tree_data$aspect, cex = 2, col="#9252FF", pch = 20)
text(y = not_tree_data$ratio, x = not_tree_data$aspect, labels = not_tree_labs, pos = 1, cex = 0.7, col="#9252FF")

## Mean number of vertebrae v mean aspect ratio

phylomorphospace(anilios_tree, anilios_vert[,c("aspect", "mean_tot_vert")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Mean number of vertebrae", xlab = "Mean aspect ratio",
                 xlim=c(0,200))
# Add slope from pgls (but this is just for species in tree)
abline(a = coefficients(ratio_meanvert.pgls)[1], b = coefficients(ratio_meanvert.pgls)[2])
# add species not in tree 
points(y = not_tree_data$mean_tot_vert, x = not_tree_data$aspect, cex = 2, col="#9252FF", pch = 20)
text(y = not_tree_data$mean_tot_vert, x = not_tree_data$aspect, labels = not_tree_labs, pos = 1, cex = 0.7, col="#9252FF")

dev.off()

### Using maximum total length individual per species

anilios_slice <- slice_df[which(slice_df$species %in% anilios_tree$tip.label),]

vert_asp.pgls <- geomorph::procD.pgls(setNames(anilios_slice$vert_ratio, anilios_slice$species) ~ setNames(anilios_slice$aspect, anilios_slice$species), phy = anilios_tree)
summary(vert_asp.pgls)
coefficients(vert_asp.pgls)

plot(y=anilios_slice$vert_ratio, x=anilios_slice$aspect, bty="n", 
     # cex = width_ratio_cex[names(width_ratio)], 
     pch = 19,
     xlab = "Max aspect ratio", ylab = "Vertebrae ratio")
text(y=anilios_slice$vert_ratio, x=anilios_slice$aspect, labels = anilios_slice$species, pos = 1, cex = 0.7)
abline(a = coefficients(vert_asp.pgls)[1], b = coefficients(vert_asp.pgls)[2])

# plot(anilios_vert$max_total_vert ~ anilios_vert$aspect, bty = "n")


### Will need to change the column names and make sure this works

### using phylolm to fit this regression but using different evolutionary models
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
## Our BM model had the lowest AIC score.

## Possible explanation: 
## This could be that these species need more vertebrae to fit through crevices 
## in hard soil? More robust species are found in softer soil (from blindsnakemorpho project). Or building larger vertebrae?


# Question 9 --------------------------------------------------------------

## Question: Does vertebrae ratio correlate with any environmental factors? 
# Environmental factors:
# mean annual temperature, 
# soil compactness.

# We are using mean total vertebrae/mean total length for each species

## PGLS using geomorph to test effect of mean annual temperature

annual_mean_temp <- setNames(anilios_data$temp_mean, anilios_data$species)
bulk_density <- setNames(anilios_data$max_bulk, anilios_data$species)

### Using geomorph::procD.pgls 
temp.pgls <- geomorph::procD.pgls(vert_ratio ~ annual_mean_temp, phy = anilios_tree)
summary(temp.pgls)
coefficients(temp.pgls)

### Using geomorph::procD.pgls 
soil.pgls <- geomorph::procD.pgls(vert_ratio ~ bulk_density, phy = anilios_tree)
summary(soil.pgls)
coefficients(soil.pgls)

### Mean total vertebrae
### Using geomorph::procD.pgls 
number_temp.pgls <- geomorph::procD.pgls(mean_vert ~ annual_mean_temp, phy = anilios_tree)
summary(number_temp.pgls)
coefficients(number_temp.pgls)

number_soil.pgls <- geomorph::procD.pgls(mean_vert ~ bulk_density, phy = anilios_tree)
summary(number_soil.pgls)
coefficients(number_soil.pgls)

# Figure 4: Combined results vertebrae ratio and max total vertebrae against ecological variables
pdf(file = "output/phylomorpho-vert-ecology2.pdf", width = 12, height = 9)
par(mfrow=c(2,2))


### VERT RATIO VS ECOLOGY

## Vert ratio v Annual mean temp
phylomorphospace(anilios_tree, anilios_data[,c("temp_mean", "ver_rati")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Vertebra ratio", xlab = "")
abline(a = coefficients(temp.pgls)[1], b = coefficients(temp.pgls)[2])

## VERT RATIO VS Soil bulk density
phylomorphospace(anilios_tree, anilios_data[,c("max_bulk", "ver_rati")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Vertebra ratio", xlab = "")
abline(a = coefficients(soil.pgls)[1], b = coefficients(soil.pgls)[2])

## VERTEBRAE NUMBER VS ECOLOGY

# Vertebrae number to anilios_data
anilios_data$aspect <- aspect_ratio[anilios_data$species]

# vert vs temperature
phylomorphospace(anilios_tree, anilios_data[,c("temp_mean", "tot_vert")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Mean vertebra number", xlab = "Mean annual temperature (°C)")
abline(a = coefficients(number_temp.pgls)[1], b = coefficients(number_temp.pgls)[2])

phylomorphospace(anilios_tree, anilios_data[,c("max_bulk", "tot_vert")], 
                 bty="n", label="horizontal", node.size=c(0.4,1.3),
                 ylab = "Mean vertebra number", xlab = expression(Max~soil~bulk~density~(g/cm^3)))
abline(a = coefficients(number_soil.pgls)[1], b = coefficients(number_soil.pgls)[2])

dev.off()




### Against aspect ratio

# Temperature
aspect_temp.pgls <- geomorph::procD.pgls(aspect_ratio ~ annual_mean_temp, phy = anilios_tree)
summary(aspect_temp.pgls)
coefficients(aspect_temp.pgls)

# Soil bulk density
aspect_soil.pgls <- geomorph::procD.pgls(aspect_ratio ~ bulk_density, phy = anilios_tree)
summary(aspect_soil.pgls)
coefficients(aspect_soil.pgls)

plot(aspect_ratio ~ annual_mean_temp, bty="n",
     pch = 19,
     xlab = "Mean annual temperature (°C)", ylab = "Mean aspect ratio")
text(x = annual_mean_temp, y = aspect_ratio, labels = sp_labs, pos = 1, cex = 0.7)
abline(a = coefficients(aspect_temp.pgls)[1], b = coefficients(aspect_temp.pgls)[2])

## Soil
plot(aspect_ratio ~ bulk_density, bty="n", pch = 19,
     xlab = expression(Max~soil~bulk~density~(g/cm^3)), ylab = "")
text(x = bulk_density, y = aspect_ratio, labels = sp_labs, pos = 1, cex = 0.7)
abline(a = coefficients(aspect_soil.pgls)[1], b = coefficients(aspect_soil.pgls)[2])


dev.off()




# Total vertebrae  -- full model - no effect
fit_phylm_total_vert <- phylolm(tot_vert ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda")
summary(fit_phylm_total_vert)

# Vertebrae ratio -- full model
fit_phylm_total_ratio <- phylolm(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda")
summary(fit_phylm_total_ratio)

## There's at least an effect of temperature, but this might be because the variables might be correlated, we can try fitting it individually.

# Stepwise model selection
phylostep(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, direction = "both", k = 2)

## Step wise fitting also indicate temperature provides the best fit. 

fit_vrat_2 <- phylolm(ver_rati ~ temp_mean, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_vrat_3 <- phylolm(ver_rati ~ mean_bulk, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_vrat_4 <- phylolm(ver_rati ~ ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
fit_vrat_5 <- phylolm(ver_rati ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)

summary(fit_vrat_2)
summary(fit_vrat_3)
summary(fit_vrat_4)
summary(fit_vrat_5)

# Visualise data
dev.off()
# par(mfrow=c(1,3))
plot(ver_rati ~ temp_mean, data = anilios_data, bty = "n", pch = 19,
     ylab = "Mean # of Vertebrae / total length (mm)", xlab = "Mean annual temperature (°C)")
abline(fit_vrat_2, lty = 3)

plot(ver_rati ~ mean_bulk, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_3, lty = 3)

plot(ver_rati ~ ARID, data = anilios_data, bty = "n", ylab = "# Vertebrae / total length (mm)")
abline(fit_vrat_4, lty = 3)

# fit_vert_2 <- phylolm(log(tot_vert) ~ temp_mean, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# fit_vert_3 <- phylolm(log(tot_vert) ~ mean_bulk, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# fit_vert_4 <- phylolm(log(tot_vert) ~ ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# fit_vert_5 <- phylolm(log(tot_vert) ~ mean_bulk + temp_mean + ARID, data = anilios_data, phy = anilios_tree, model = "lambda", boot = 500)
# 
# sum_vert_2 <- summary(fit_vert_2); sum_vert_3 <- summary(fit_vert_3); sum_vert_4 <- summary(fit_vert_4); sum_vert_5 <- summary(fit_vert_5)
# 
# sum_vert_2 
# sum_vert_3
# sum_vert_4
# sum_vert_5

# par(mfrow=c(2,2))
# plot(tot_vert ~ temp_mean, data = anilios_data)
# plot(tot_vert ~ mean_bulk, data = anilios_data)
# plot(tot_vert ~ ARID, data = anilios_data)

## Answer: vertebrae ratio correlate with temperature only, not aridity or soil compactness

