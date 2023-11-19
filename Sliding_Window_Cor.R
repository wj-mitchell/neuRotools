# Sliding_Window_Cor.R || v2023.11.15

Sliding_Window_Cor <- function(x,
                               y,
                               window_size,
                               sigma = 3,
                               step_size = 1,
                               cor_use = "complete.obs",
                               cor_method = "spearman"){
  
# ----- QA CHECKS -----
  if (length(x) != length(y)){
    stop("The length of x and y do not match. Please submit two arrays of equal length and try again")
  }
  
  if (any(!is.numeric(x) | !is.numeric(y))){
    stop("X and Y must be numeric arrays of equal length. Please consider coercing or reformatting your values and try again")
  }
  
  if (!is.numeric(window_size) | window_size < 2){
    stop("Window size must be a numeric value greater than 1. Please correct this and try again")
  }
  
  if (!is.numeric(sigma) | sigma <= 0){
    stop("Sigma must be a non-zero numeric value. Please correct this and try again")
  }
  
  if (!is.numeric(step_size) | step_size < 1){
    stop("Step size must be a non-zero numeric value. Please correct this and try again")
  }

# ----- SETUP -----
  # Checking for pacman
  if (require("pacman") == FALSE){
    install.packages("pacman")
  }
  
  # Loading in my packages with my pacman manager
  pacman::p_load(stats)
  
  # Adding custom functions
  source("circle_shift.R")

# ----- DATA PREP -----

  # Calculating how many volumes there are 
  nVols <- length(x)

# ----- GENERATING WINDOW -----
  
# Generating the median and series indices 
  if ((nVols %% 2) != 0){
    median <- ceiling(nVols/2)
    series_index <- 0:nVols
  }
  if ((nVols %% 2) == 0){
    median <- nVols/2
    series_index <- 0:(nVols-1)
  }
  
# Defining the radius of the sliding window
  window_radius <- round(window_size/2)

# Creating a Gaussian distribution around the median of the series
  gauss_dist <- exp(-((c(series_index-median)^2) / (2 * sigma^2)))

# Creating an empty series on which to apply the window
  series_window <- rep(0, nVols)
  
# Applying the window (i.e., changing 0's to 1's)}
  series_window[(median - window_radius + 1):(median + window_radius)] <- 1

# Convolving the window we want to target with a Gaussian distribution
  convol <- convolve(y = gauss_dist, 
                     x = series_window, 
                     type = "open")
  
# Standardizing the convolution so that it peaks at 1
  convol <- convol/max(convol)                      
  
# Centering the series so that the peak is at the median
  convol <- convol[(median + 1):(length(convol) - median + 1)] 
  
# Shaving off any excess timepoints, in the event some are at the edges
  convol <- convol[1:nVols]             
  
# ----- APPLYING THE WINDOW -----                     
                     
# Defining how many unique window that we will have, defined as the number of volumes, minus how large the window radius is and divided by how far apart different iterations of windows will be
# Note that I'm using the window radius instead of the window_size; without doing that, the beginning and end timepoints are downweighted much more severely than any other timepoint on average. By using the radius, the beginning begins at a weight around 1 and end ends at a weight around 1. It might be easier to visualize with plot() if you're having trouble imagining what I'm saying 
  nWindow <- (nVols - window_radius) / step_size
  
# Creating an empty dataframe to house the sliding window 
  cor_sw <- rep(NA, nWindow)

# Identifying the indices around which each iteration of the window should center
  indices <- (seq(1, nWindow, step_size) + (window_radius/2))
  
# Iterating through the different windows
  for (WINDOW in indices){
    
    # Shifting the window iteratively
    convol_shift <- circle_shift(convol, 
                                 median - WINDOW)
    
    # Identify where the window cutoffs should occur
    cutoffs <- c(WINDOW - window_radius, WINDOW + window_radius)
    if (any(cutoffs < 1)){
      cutoffs[which(cutoffs < 1)] <- 1
    }
    if (any(cutoffs > nVols)){
      cutoffs[which(cutoffs > nVols)] <- nVols
    }
    
    # If the convol is spilling over from the front to the back
    # If the iteration we're on is in the first half of indices, but there are non-zero values at the tail
    if (WINDOW < median & convol_shift[length(convol_shift)] != 0){
      convol_shift[median:length(convol_shift)] <- 0
      convol_shift <- convol_shift * (sum(convol)/sum(convol_shift))
    }
    # If the convol is spilling over from the back to the front
    # If the iteration we're on is in the latter half of indices, but there are non-zero values at the start
    if (WINDOW > median & convol_shift[1] != 0){
      convol_shift[1:median] <- 0
      convol_shift <- convol_shift * (sum(convol)/sum(convol_shift))
    }
    
    # Convolving data and shaving off 0's
    data_x <- x * convol_shift 
    data_x <- data_x[cutoffs[1]:cutoffs[2]]
    data_y <- y * convol_shift 
    data_y <- data_y[cutoffs[1]:cutoffs[2]]
    
    # Generating Correlations
    cor_sw[which(indices == WINDOW)] <- cor(x = data_x,
                                            y = data_y,
                                            method = cor_method,
                                            use = cor_use) 
  }
  
# ----- GENERATING OUTPUT -----
  
# Return the correlation values
  return(cor_sw)
}
  
  
