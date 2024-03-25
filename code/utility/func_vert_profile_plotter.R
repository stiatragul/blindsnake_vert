# func_vert_profile_plotter.R

# Function to plot vertebral profile

vert_profile_plotter <- function(.data1, .data2, .data3, .title = NULL){
  
  rep1 <- .data1
  rep2 <- .data2
  rep3 <- .data3
  
  # Calculate the lengths of each vertebra from the landmark coordinates
  lindist1 <- diag(as.matrix(dist(rep1))[-1,]) 
  lindist2 <- diag(as.matrix(dist(rep2))[-1,])
  lindist3 <- diag(as.matrix(dist(rep3))[-1,])
  
  # Combine into a single matrix
  lindist <- cbind(lindist1, lindist2, lindist3)
  
  # Average the vertebrae lengths across the three replicates
  avg <- apply(lindist, 1, mean)
  
  # Calculate the standard deviation for each vertebra length
  sdev <- apply(lindist, 1, sd)
  
  # Calculate the standard error for each vertebra length
  sdevm <- sdev / sqrt(nrow(lindist))
  mean(sdevm) # calculate your repeatability score (closest to 0 wins!)
  
  # Calculate the body length, which is the sum of the measurements you've taken
  sum(avg)
  
  # Get the number of vertebrae you've measured
  length(avg)
  
  # Plot a graph to show the pattern of intracolumnar variation in vertebrae size
  plot(avg, pch=19, xlab = "vertebra number", 
       ylab= "vertebra length (mm)", bty = "n",
       ylim=c(0, max(avg+sdev)))
  # Add error bars
  arrows(c(1:length(avg)), avg-sdevm, 
         c(1:length(avg)), avg+sdevm, 
         length=0.05, angle=90, code=3)
  # Add a title
  title(paste(.title)) 
  # change "My Name & Species name & body length & Number Vertebrae" to your info
  
}