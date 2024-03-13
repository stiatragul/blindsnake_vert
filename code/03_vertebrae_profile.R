# 03_vertebrae_profile.R
# Putter Tiatragul
# March 2024
# Visualise vertebrae length by numnber of vertebrae

# load data ---------------------------------------------------------------

ligatus_1 <- read.csv('data/x_ray/vert_lengths.csv')
grypus_1 <- read.csv('data/x_ray/vert_length_grypus.csv')
 

# Processing --------------------------------------------------------------

ligatus_1 <- ligatus_1[c("Label", "Length")]
ligatus_1$Label <- sub("\\.dcm$", "", ligatus_1$Label)
ligatus_1$vert_no <- ave(seq_along(ligatus_1$Label), ligatus_1$Label, FUN = seq_along)
ligatus_1$rel_length <- ligatus_1$Length / mean(ligatus_1$Length)


grypus_1 <- grypus_1[c("Label", "Length")]
grypus_1$Label <- sub("\\.dcm$", "", grypus_1$Label)
grypus_1$rel_length <- grypus_1$Length / mean(grypus_1$Length)
grypus_1$vert_no <- ave(seq_along(grypus_1$Label), grypus_1$Label, FUN = seq_along)

lig_1 <- ligatus_1[ligatus_1$Label == unique(ligatus_1$Label)[1], ]
gry_1 <- grypus_1[grypus_1$Label == unique(grypus_1$Label)[1], ]



# Plotting ----------------------------------------------------------------
par(mfrow=c(1,2))

plot(rel_length ~ vert_no, data = lig_1, pch = 20, xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n",
     ylim = c(0, 1.6))

plot(rel_length ~ vert_no, data = gry_1, pch = 20, xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n",
     ylim = c(0, 1.6))

