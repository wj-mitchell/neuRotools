# A function which will generate a gamma distribution as an array of a specific length
# Constructed specifically to be applied to lag a parametric modulator when the construct the modulator represents might occur across undetermined timepoints prior to the trial of interest

gammify <- function(bins = 10, # The length of the outputted array
                    n = 100000, # From the rgamma function, the number of datapoint samples generated for the distribution temple 
                    shape = 1.5, # From the rgamma function, specifies the shape of the distribution template 
                    rate = 0.5, # From the rgamma function, specifies the rate of the distribution template 
                    seed = 2, # Specifying a seed to reliably generate the distribution template
                    invert = T, # If true, distribution will be negatively skewed 
                    shave_tail = T, # Whether a proportion of tail values should be shaved off to focus specifically on the majority of the distribution
                    cutoff_tail = 0.65) # A value between 0 and 1 indicating what proportion of the tail shove be shaved off
{
  # Setting a seed for replication
  if (!is.na(seed)){
    set.seed(seed)
  }

  # Generating a robust Gamma Distribution
  distribution <- rgamma(n = n, 
                         shape = shape, 
                         rate = rate)
  
  # Inverting the distribution
  if (invert == TRUE){
    distribution <- ceiling(max(distribution)) - distribution
  }

  # Shave insignificant values
  if (shave_tail == TRUE & invert == TRUE){
    distribution <- distribution[distribution > ceiling(max(distribution)) * cutoff_tail]
  }
  if (shave_tail == TRUE & invert == FALSE){
    distribution <- distribution[distribution < ceiling(max(distribution)) * (1 - cutoff_tail)]
  }
  
  # Generating a new dataframe to house transformation
  df <- as.data.frame(matrix(NA, 
                             nrow = bins + 2,
                             ncol = 2,
                             dimnames = list(1:(bins + 2), 
                                             c("value", "prop"))))
  
  # Creating an array of equally sized bins which compose the entire range of our distribution
  df$value[2:nrow(df)] <- seq(floor(min(distribution)), 
                              ceiling(max(distribution)),
                              (ceiling(max(distribution)) - floor(min(distribution)))/bins)
  
  # Iterating through each of our bins and identifying the proportion of values from the distribution which fit in that bin
  for (BIN in 2:(nrow(df) - 1)){
    df$prop[BIN] <- length(which(distribution >= df$value[BIN] & distribution < df$value[BIN + 1])) / length(distribution)
  }

  # Specifying the first and last values should be 0's
  df$prop[1] <- 0
  df$prop[bins + 2] <- 0
  
  # Changes the values of value to be meaningful
  df$value <- -(nrow(df) - 1):0
  
  # Normalizing Proportion values such that the average is 1
  df$prop <- nrow(df) * df$prop
  
  return(df$prop)
}

# A demonstration of the plot of the function and also that the average value is 1 
# plot(gammify())
# mean(gammify())
