# archive_03_vertebrae_profile.R
# Putter Tiatragul
# March 2024
# Visualise vertebrae length by numnber of vertebrae

# load data ---------------------------------------------------------------

ligatus_1 <- read.csv('data/x_ray/vert_lengths.csv')
grypus_1 <- read.csv('data/x_ray/vert_length_grypus.csv')

# ligatus_2 <- read.csv('data/x_ray/vert_lengths_2.csv')
# vert_df2 <- rbind(ligatus_2)


# Processing --------------------------------------------------------------

vert_df <- vert_df[c("Label", "Length")]
vert_df$Label <- sub("\\.dcm$", "", vert_df$Label)
vert_df$vert_no <- ave(seq_along(vert_df$Label), vert_df$Label, FUN = seq_along)
vert_df$rel_length <- vert_df$Length / mean(vert_df$Length)
vert_df$trial <- 1

vert_df2 <- vert_df2[c("Label", "Length")]
vert_df2$Label <- sub("\\.dcm$", "", vert_df2$Label)
vert_df2$rel_length <- vert_df2$Length / mean(vert_df2$Length)
vert_df2$vert_no <- ave(seq_along(vert_df2$Label), vert_df2$Label, FUN = seq_along)
vert_df2$trial <- 2

vert_df <- rbind(vert_df, vert_df2)

# Subset by sample
lig_1 <- vert_df[vert_df$Label == unique(vert_df$Label)[1], ]
lig_2 <- vert_df[vert_df$Label == unique(vert_df$Label)[2], ]
lig_3 <- vert_df[vert_df$Label == unique(vert_df$Label)[3], ]
gry_1 <- vert_df[vert_df$Label == unique(vert_df$Label)[4], ]
gry_2 <- vert_df[vert_df$Label == unique(vert_df$Label)[5], ]

dev.off()
par(mfrow=c(1,3))

plot(rel_length ~ vert_no, data = lig_1, 
     col = setNames(viridis::viridis(2), unique(lig_1$trial)), pch = c(16, 17),
     xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n")

plot(rel_length ~ vert_no, data = lig_2, 
     col = setNames(viridis::viridis(2), unique(lig_2$trial)), pch = c(16, 17),
     xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n")

plot(rel_length ~ vert_no, data = lig_3, 
     # col = setNames(viridis::viridis(2), unique(lig_3$trial)), pch = c(16, 17),
     xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n")


plot(rel_length ~ vert_no, data = gry_1, 
     xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n")

plot(rel_length ~ vert_no, data = gry_2, 
     xlab = "Vertebrae no", ylab = "Rel. vert. length", bty = "n")

# Colours -----------------------------------------------------------------

# Generate a vector of colours using some palette 
num_colours <- length(unique(vert_df$Label))
uniq_colours <- viridis::viridis(6)

# Assign colours based on unique values in column

colour_map <- setNames(uniq_colours, unique(vert_df$Label)) 
vert_df$colours <- colour_map[vert_df$Label]

# Plot --------------------------------------------------------------------

plot(Length ~ vert_no, data = vert_df, 
     col = vert_df$colours,
     pch = setNames(c(15, 16, 17), unique(vert_df$Label)),
     xlab = "Vertebrae no", ylab = "Vertebrae length", bty = "n")


lig_1 <- vert_df[vert_df$Label == unique(vert_df$Label)[1], ]


plot(Length ~ vert_no, data = vert_df[vert_df$Label == unique(vert_df$Label)[1], ], 
     col = vert_df$colours,
     pch = setNames(c(15, 16, 17), unique(vert_df$Label)),
     xlab = "Vertebrae no", ylab = "Vertebrae length", bty = "n")


### Choosing one to plot 
par(mfrow = c(1,3))

for (i in 1:length((unique(vert_df$Label)))) {
  
  subset_df <- vert_df[which(vert_df$Label == unique(vert_df$Label)[i]),]
  
  # Plot the data with same scale for x and y axes
  plot(subset_df$Length ~ subset_df$vert_no, 
       col =  vert_df$colours,
       pch = 15,
       xlab = "Vertebrae no", ylab = "Vertebrae length", bty = "n",
       ylim = range(vert_df$Length), xlim = range(vert_df$vert_no))
  
}

