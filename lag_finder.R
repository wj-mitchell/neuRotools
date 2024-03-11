lag_finder <- function(lags = seq(-40,40,1), # all of the possible lags to test
                       data_x, # The dataframe that this data exists within
                       data_y, # The dataframe that this data exists within
                       time, # The time variable which you wish to apply the lags to
                       x, # The variable to add a time-shift and correlate with your y variable
                       y, # The variable to correlate with your x variable
                       window_size = NA, # the size of the sliding window that we should apply  
                       method = "spearman") # Which method that you want cor() to use
{ 
  
  
  # ----- QA CHECKS -----
  
  if (length(data_x[[x]]) != length(data_y[[y]])){
    stop("The length of x and y do not match. Please submit two arrays of equal length and try again")
  }
  
  if (any(!is.numeric(data_x[[x]]) | !is.numeric(data_y[[y]]))){
    stop("X and Y must be numeric arrays of equal length. Please consider coercing or reformatting your values and try again")
  }
  
  if ((!is.numeric(window_size) | window_size < 2) & !is.na(window_size)){
    stop("Window size must be NA OR a numeric value greater than 1. Please correct this and try again")
  }
  
  if (!is.numeric(data_x[[time]]) | !is.numeric(data_y[[time]])){
    stop("Time must be a numerically structured variable in dataframes x and y. Please try coerecing them into the correct format and try again")
  }
  
  if (!is.numeric(lags) | length(lags) < 2){
    stop("The lags you submit must be a numerically structured array of 2 or more values. Please correct this and try again")
  }
  
  if (method != "spearman" | method != "pearson"| method != "kendall"){
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
    rows <- (nrow(data_x) - (window_size / 2))
    cols <- length(lags)
    values <- as.data.frame(matrix(data = NA, 
                                   nrow = rows, 
                                   ncol = cols,
                                   dimnames = list(paste0("Win_",1:rows),
                                                   paste0("Lag_", lags))))
  }
  
  # Noting which column is the time column
  col_time <- which(names(data_x) == time)
  
  # Noting time mins and maxes
  max_time <- max(data_x[,col_time])
  min_time <- min(data_x[,col_time])
  
  # Iterating through all possible lags
  for (LAG in 1:length(lags)){
    
    # Creating scratch dataframes to work within
    data_temp_x <- data_x
    
    # Adding the lag to the behavioral data
    data_temp_x[,col_time] <- data_temp_x[,col_time] + lags[LAG]
    
    # If there actually was a lag
    if (lags[LAG] != 0){
      
      # Shaving off any observations that now fall outside of the Window range
      data_temp_x <- data_temp_x[-which(data_temp_x[,col_time] < min_time | data_temp_x[,col_time] > max_time),] 
    }
    
    # Merging out behavioral and neural data together
    df <- merge(data_temp_x, 
                data_y)
    
    # Noting which column are the correlation columns
    col_x <- which(names(df) == x)
    col_y <- which(names(df) == y)
    
    # If window has no value
    if (is.na(window_size)){
      
      # Generating a correlative value for this region
      values[LAG] <- cor(df[,col_x], 
                         df[,col_y], 
                         method = method)
    }
    
    # If we want to use the sliding window approach
    if (!is.na(window_size)){
      
      # Generate a series of correlative value for this region
      values[,LAG] <- Sliding_Window_Cor(x = df[,col_x], 
                                         y = df[,col_y],
                                         window_size = window_size,
                                         cor_method = method)
      
    }
  }
  return(values)
}
