#'
#' @title sliding_window_cor.R || v2025.06.12
#'
#' @description This function calculates a time-resolved correlation between two signals using a sliding window with optional Gaussian weighting and parallelization.
#'
#' @param x (numeric vector) First time series to be correlated; must be the same length as `y`.
#' @param y (numeric vector) Second time series to be correlated; must be the same length as `x`.
#' @param window (integer) The width of the sliding window (in timepoints) over which the correlation is computed.
#' @param sigma (numeric) Standard deviation of the Gaussian kernel used to weight values within each window; controls temporal smoothing.
#' @param step (integer) Step size between the centers of successive windows; determines the resolution of the time-resolved correlation.
#' @param use (character) Argument passed to the `use` parameter in the `cor()` function; typically "complete.obs" or "pairwise.complete.obs" to handle missing data.
#' @param method (character) Correlation method to use; one of "pearson", "spearman", or "kendall".
#' @param zeroes_to_na (logical) If TRUE, correlation will return NA in windows where either `x` or `y` has zero variance; if FALSE, such cases will return 0 or the result from `cor()`.
#' @param n_cores (integer) Number of processor cores to use for parallel computation; defaults to half the available cores.
#' @param parallelize (logical) If TRUE (default), uses parallel computation with `foreach`; if FALSE, runs in serial.
#'
#' @return (numeric vector) A vector of correlation values, one for each sliding window, aligned to the center of each window segment.
#'
#' @export
sliding_window_cor <- function(x,
                               y,
                               window,
                               sigma = 3,
                               step = 1,
                               use = "complete.obs",
                               method = "spearman",
                               zeroes_to_na = FALSE,
                               n_cores = parallel::detectCores() / 2,
                               parallelize = FALSE) {
  
  require(stats)
  source("https://raw.githubusercontent.com/wj-mitchell/neuRotools/main/circle_shift.R", local = FALSE)
  
  # ----- QUALITY CONTROL CHECKS -----
  
  # Confirm x and y are equal in length
  if (length(x) != length(y)){
    stop("The length of x and y do not match. Please submit two arrays of equal length and try again")
  }
  
  # Confirm that x and y are numeric
  if (any(!is.numeric(x) | !is.numeric(y))){
    stop("x and y must be numeric arrays of equal length. Please consider coercing or reformatting your values and try again")
  }
  
  # Confirm that window size is numeric and > 1
  if (!is.numeric(window) | window < 2){
    stop("Window size must be a numeric value greater than 1. Please correct this and try again")
  }
  
  # Confirm that sigma is numeric and > 0
  if (!is.numeric(sigma) | sigma <= 0){
    stop("Sigma must be a non-zero numeric value. Please correct this and try again")
  }
  
  # Confirm that step size is numeric and >= 1
  if (!is.numeric(step) | step < 1){
    stop("Step size must be a non-zero numeric value. Please correct this and try again")
  }
  
  # Confirm that zeroes_to_na is logical
  if (!is.logical(zeroes_to_na)){
    stop("zeroes_to_na must be either TRUE or FALSE. Please correct this and try again")
  }
  
  # ----- PREP: DEFINE CORE VARIABLES -----
  
  # Determine how many timepoints we have
  nVols <- length(x)
  
  # Define the radius of the window (half its size)
  window_radius <- round(window / 2)
  
  # Define the median of the series for alignment
  median <- ceiling(nVols / 2)
  
  # ----- CREATE A GAUSSIAN WINDOW FOR WEIGHTING -----
  convol <- generate_convolution_window(nVols, window_radius, sigma)
  
  # ----- DEFINE THE WINDOW CENTERS AND ITERATION PLAN -----
  
  # Defining how many unique window that we will have, defined as the number of volumes, minus how large the window radius is and divided by how far apart different iterations of windows will be
  # Note that I'm using the window radius instead of the window; without doing that, the beginning and end timepoints are downweighted much more severely than any other timepoint on average. By using the radius, the beginning begins at a weight around 1 and end ends at a weight around 1. It might be easier to visualize with plot() if you're having trouble imagining what I'm saying 
  nWindow <- floor((nVols - window_radius) / step)
  
  # Define the center indices of each window
  indices <- seq(1, by = step, length.out = nWindow) + (window_radius / 2)
  
  # ----- SLIDING WINDOW CORRELATION -----
  if (parallelize) {
    # Register parallel cluster
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    
    # Compute correlations in parallel
    cor_sw <- foreach::foreach(i = seq_along(indices), .combine = c, .packages = "stats") %dopar% {
      WINDOW <- indices[i]
      shift_amt <- median - WINDOW
      compute_weighted_correlation(x, y, convol, shift_amt, WINDOW, window_radius, median,
                                   zeroes_to_na, method, use)
    }
    
    # Shut down the parallel cluster
    parallel::stopCluster(cl)
    
  } else {
    # Compute correlations in serial
    cor_sw <- numeric(length(indices))
    
    for (i in seq_along(indices)) {
      WINDOW <- indices[i]
      shift_amt <- median - WINDOW
      cor_sw[i] <- compute_weighted_correlation(x, y, convol, shift_amt, WINDOW, window_radius, median,
                                                zeroes_to_na, method, use)
    }
  }
  
  # ----- RETURN THE SLIDING WINDOW CORRELATION SERIES -----
  return(cor_sw)
}

#' @keywords internal
#' @description Helper function to generate Gaussian-weighted convolution window
#' @return (numeric vector) normalized Gaussian convolution window of length nVols
generate_convolution_window <- function(nVols, window_radius, sigma) {
  # Define the center of the signal
  median <- ceiling(nVols / 2)
  
  # Generate a Gaussian distribution centered on the median index
  gauss_dist <- exp(-((0:(nVols - 1) - median)^2) / (2 * sigma^2))
  
  # Create a binary window of 1s centered around the median
  series_window <- rep(0, nVols)
  series_window[(median - window_radius + 1):(median + window_radius)] <- 1
  
  # Convolve the binary window with the Gaussian to apply weighting
  convol <- convolve(gauss_dist, series_window, type = "open")
  
  # Normalize so the peak value is 1
  convol <- convol / max(convol)
  
  # Trim the convolution to match the original series length
  convol <- convol[(median + 1):(length(convol) - median + 1)][1:nVols]
  
  return(convol)
}

#' @keywords internal
#' @description Helper function to apply weighted convolution and compute correlation
#' @return (numeric) correlation value for a given window
compute_weighted_correlation <- function(x, y, convol, shift_amt, window_center,
                                         window_radius, median, zeroes_to_na,
                                         method, use) {
  # Determine length of the signal
  nVols <- length(x)
  
  # Apply a circular shift to center the convolution
  convol_shift <- circle_shift(convol, shift_amt)
  
  # If the window spills off the back end (early in the series)
  if (window_center < median && convol_shift[nVols] != 0) {
    convol_shift[median:nVols] <- 0
    convol_shift <- convol_shift * (sum(convol) / sum(convol_shift))
  }
  
  # If the window spills off the front end (late in the series)
  if (window_center > median && convol_shift[1] != 0) {
    convol_shift[1:median] <- 0
    convol_shift <- convol_shift * (sum(convol) / sum(convol_shift))
  }
  
  # Apply the shifted convolution window to both signals
  data_x <- x * convol_shift
  data_y <- y * convol_shift
  
  # Trim both signals to only the cutoff region
  cutoffs <- c(max(1, window_center - window_radius), min(nVols, window_center + window_radius))
  data_x <- data_x[cutoffs[1]:cutoffs[2]]
  data_y <- data_y[cutoffs[1]:cutoffs[2]]
  
  # Check for complete cases and enough data to compute correlation
  complete_idx <- complete.cases(data_x, data_y)
  
  if (sum(complete_idx) >= 2 &&
      sd(data_x[complete_idx]) != 0 &&
      sd(data_y[complete_idx]) != 0) {
    
    # Compute the correlation on complete, non-constant data
    return(cor(data_x[complete_idx], data_y[complete_idx], method = method, use = use))
    
  } else if (!zeroes_to_na) {
    
    # If user doesn't want to convert zeroes to NA, try anyway (at their own risk)
    return(cor(data_x, data_y, method = method, use = use))
    
  } else {
    
    # Not enough data or all-zero, return NA
    return(NA)
  }
}
