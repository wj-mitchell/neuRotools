# HRF_convolve.R || v 2023.11.15

# I tried using the convolve function from that stats package to complete a linear convolution of an HRF function with a parametric modulator, but I found that the results were a little confusing. When I did the hand-calculations to check, the results were not consistent with a straight-forward linear convolution so I figured it would be simple enough to just make my own function to do this. If you feed it a series and some weights, it will add an HRF event to every non-zero entry in the series. Overlapping HRF functions are additive. The standard weights sum to 1, but will be scaled to the series value, so if the series value that launches the HRF event is 2, then the value of each weight will be doubled before being applied to the series. Lastly, if we have a TR that is different than 1, we can resample these weights by modifying either p or q. For example, if I change q to 2, reflecting a TR of 2 seconds, then these weights will be downsampled to 8 values rather than the standard 15. If I change q to 0.5, I will upsample to a weight series of 30 values. 

HRF_convolve <- function(series, 
                         weights = c(0, 0.000354107958396228, 0.0220818694830938, 0.116001537001027, 0.221299059999514, 0.242353095826523, 0.186831750619196, 0.113041009515928, 0.0572809597709863, 0.0253492394574127, 0.0100814114758446, 0.00367740475539297, 0.00124901102357508, 0.000399543113110135, 0),
                         resample_p = 1,
                         resample_q = 1){
  
  # ----- QA CHECKS -----
  
  # Making sure that series is a numeric array
  if (any(!is.numeric(series))){
    stop("Series must be an array of numeric values. Please try again.")
  }
  
  # Making sure that weights is a numeric array
  if (any(!is.numeric(weights))){
    stop("The convolution weights must be an array of numeric values. Please try again.")
  }
  
  # Making sure that p is is a numeric value
  if (!is.numeric(resample_p)){
    stop("Resample_p must be a single numeric value. Please try again.")
  }
  
  # Making sure that q is is a numeric value
  if (!is.numeric(resample_q)){
    stop("Resample_q must be a single numeric value. Please try again.")
  }
  
  # ----- DEPENDENCIES -----
  if (require("pacman") == FALSE){
    install.packages("pacman")
  }
  
  # Loading in my packages with my pacman manager
  pacman::p_load(signal)
  
  # ----- RESAMPLING -----
  
  # If we have an incongurent p and q ...
  if (resample_p != resample_q){
    
    # Resample our weights
    weights <- resample(weights, 
                        p = resample_p, 
                        q = resample_q)
  }

  # Check that we don't have more weights than we do things to apply them to
  if (length(weights) > length(series)){
    stop("You have more convolutional weights than you have values in your series. Please add a larger series or reduce the number of weights that you have")
  }
  
  # Identifying which indices in the series are calling for an HRF convolution
  indices <- which(series != 0)
  
  # Pulling the values of those indices
  values <- series[indices]
  
  # Creating an empty array of zeroes on which to add the convolution
  convolution <- rep(0, length(series))
  
  # Iterating through each of the non-zero series events 
  for (INDEX in 1:length(indices)){
    
    # Calculating how much to modulate the weights by, according to the value of the event
    weights_mod <- values[INDEX] * sum(weights)
    
    # Noting which timepoints in the series will have their value modified by this convolution
    window <- indices[INDEX]:(indices[INDEX] + length(weights) - 1)
    
    # However, if the max value in the window is longer than the length of the series 
    if (max(window) > length(convolution)){
      
      # Reduce the length of the window to match the length of the series
      window <- window[-(which(window > length(convolution)))]
    }
    
    # Then iterate through each timepoint within the window
    for (TIMEPOINT in window){
      
      # and make it equal to it's current value plus the weight for that position times its modifier
      convolution[TIMEPOINT] <- convolution[TIMEPOINT] + (weights[which(window == TIMEPOINT)] * weights_mod)
    }
  }
  
  # Return the convolutional series
  return(convolution)
}
  