# 
#' @title circle_shift.R | v 2025.04.30
#'
#' @description Efficiently performs a circular (cyclic) shift on a numeric vector
#'
#' @param series (numeric vector) a numeric vector to be circularly shifted
#'
#' @param shift (integer) an integer indicating how many positions to shift (positive = right, negative = left)
#'
#' @return (numeric vector) a shifted numeric vector
#'
#' @export
#'

circle_shift <- function(series, shift) {
  
  # Get the length of the input vector
  n <- length(series)
  
  # Ensure the shift is within bounds (e.g., shift of 10 on a 5-element vector becomes 0)
  shift <- shift %% n
  
  # If no shift is needed, return the original series immediately
  if (shift == 0) return(series)
  
  # Preallocate a numeric vector of the same length for performance
  shifted <- numeric(n)
  
  # Assign the "end" of the series to the beginning of the shifted version
  # For example, if shift = 2, move the last (n - 2) elements to the front
  shifted[1:(n - shift)] <- series[(shift + 1):n]
  
  # Assign the "beginning" of the series to the end of the shifted version
  # For shift = 2, move the first 2 elements to the end
  shifted[(n - shift + 1):n] <- series[1:shift]
  
  # Return the shifted vector
  return(shifted)
}