# func_vert_profile_plotter.R

# Function to plot vertebral profile

vert_raw_plotter <- function(.data1, .data2, .data3, .title = NULL){
  
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


vert_profile_plotter <- function(.data1, .data2, .data3, .title = NULL, .error = FALSE){
  
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
  
  # Calculate the standard error for each relative vertebra length
  sdevm <- (sdev/sum(sdev)) / sqrt(nrow(lindist))
  
  # Calculate relative size
  relative_size <- avg / sum(avg)
  
  
  # Fit polynomial 4 degree
  polynomial_fit <- lm(relative_size ~ poly(1:length(relative_size), degree = 4, raw = TRUE))
  
  # Generate equally spaced points along the x-axis
  x_values <- seq(1, length(relative_size), length.out = 100)
  
  # Predict corresponding y-values using the fitted polynomial
  predicted_points <- predict(polynomial_fit, newdata = data.frame(x_values))
  
  # Plot the profile
  plot(relative_size, xlab = "Vertebra number", ylab = "Relative vertebra size", ylim=c(0,0.01), xlim=c(0, 400))
  title(paste(.title)) 
  
  # Add error bars if .error is TRUE
  if (.error) {
    arrows(c(1:length(relative_size)), relative_size-sdevm, 
           c(1:length(relative_size)), relative_size+sdevm, 
           length=0.05, angle=90, code=3)
  }
  
  lines(predicted_points, type = "l", col = "red")
  # legend("topleft", legend = c(" Profile", "Fitted Polynomial"), col = c("black", "red"), lty = c(NA, 1), pch = c(1, NA))
  
}




