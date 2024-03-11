lag_finder <- function(time, # An array of numeric values which you wish to apply the lags to
                       x, # An array of numeric values to add a time-shift and correlate with your y variable
                       y, # An array of numeric values to correlate with your x variable
                       lags = seq(-40,40,1), # all of the possible lags to test
                       window_size = NA, # the size of the sliding window that we should apply  
                       method = "spearman") # Which method that you want cor() to use
{ 
  
  
  # ----- QA CHECKS -----
  
  if (length(x) != length(y)){
    stop("The length of x and y do not match. Please submit two arrays of equal length and try again")
  }
  
  if (any(!is.numeric(x) | !is.numeric(y))){
    stop("X and Y must be numeric arrays of equal length. Please consider coercing or reformatting your values and try again")
  }
  
  if ((!is.numeric(window_size) | window_size < 2) & !is.na(window_size)){
    stop("Window size must be NA OR a numeric value greater than 1. Please correct this and try again")
  }
  
  if (!is.numeric(time)){
    stop("Time must be a numerically structured array. Please try coerecing them into the correct format and try again")
  }
  
  if (!is.numeric(lags) | length(lags) < 2){
    stop("The lags you submit must be a numerically structured array of 2 or more values. Please correct this and try again")
  }
  
  if (method != "spearman" & method != "pearson" & method != "kendall"){
    stop("This function depends upon the standard cor() function in base r. Acceptable methods include any of the methods accepted by that function. The method you submitted does not match any recognized options. Please review cor() documentation, make a correction, and try again.")
  }
  
  # ----- SETUP -----
  
  # Checking for pacman
  # if (require("pacman") == FALSE){
  #   install.packages("pacman")
  # }
  
  # Loading in my packages with my pacman manager
  # pacman::p_load()
  
  # If we want to use the sliding window function
  if (!is.na(window_size)){ 
    
    ## Adding custom function
    source("https://raw.githubusercontent.com/wj-mitchell/neuRotools/main/Sliding_Window_Cor.R", local= T)
  
    }
  
  # ----- FUNCTION -----
  
  # If window has no value
  if (is.na(window_size)){
  
    # Creating an empty array to house correlative values
    values <- NULL
  
  }
  
  # If we want to use the sliding window function
  if (!is.na(window_size)){
    
    # Creating an empty dataframe to house correlative values
    rows <- (length(x) - (window_size / 2) - abs(lags[which.max(abs(lags))]))
    cols <- length(lags)
    values <- as.data.frame(matrix(data = NA, 
                                   nrow = rows, 
                                   ncol = cols,
                                   dimnames = list(paste0("Win_",1:rows),
                                                   paste0("Lag_", lags))))
  }
    
  # Noting time mins and maxes
  max_time <- max(time)
  min_time <- min(time)

  # Creating a dataframe to house the Y variables
  data_y <- data.frame(time = time,
                       y = y)
  
  # Iterating through all possible lags
  for (LAG in 1:length(lags)){
    
    # Creating scratch dataframes to work within
    data_temp_x <- data.frame(time = time,
                              x = x)
    
    # Adding the lag to the behavioral data
    data_temp_x[["time"]] <- data_temp_x[["time"]] + lags[LAG]
    
    # If there actually was a lag
    if (lags[LAG] != 0){
      
      # Shaving off any observations that now fall outside of the Window range
      data_temp_x <- data_temp_x[-which(data_temp_x[["time"]] < min_time | data_temp_x[["time"]] > max_time),] 
    }
    
    # Merging out behavioral and neural data together
    df <- merge(data_temp_x, 
                data_y)
    
    # If window has no value
    if (is.na(window_size)){
      
      # Generating a correlative value for this region
      values[LAG] <- cor(df[["x"]], 
                         df[["y"]], 
                         method = method)
    }
    
    # If we want to use the sliding window approach
    if (!is.na(window_size)){
      
      # Generate a series of correlative value for this region
      values[,LAG] <- Sliding_Window_Cor(x = df[["x"]], 
                                         y = df[["y"]],
                                         window_size = window_size,
                                         cor_method = method)
      
    }
  }
  return(values)
}

# ----- TEST -----
# If this works like it should, we should see a correlation of 1 at point 0 and mirrored, randomly determined correlations on either side of lag 0, regardless of what values x takes.

# time <- 1:20
# lags <- seq(-4,4,1)
# x <- sample(1:1000,20,replace = TRUE)
#  
# values <- lag_finder(time = time,
#                      x = x,
#                      y = x,
#                      lags = lags,
#                      method = "pearson")
