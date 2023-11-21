# circle_shift.R | v 2023.11.15

circle_shift <- function(series, 
                         shift) {
    
    # Ensuring that the shift is not greater than the series length
    shift <- shift %% length(series)
    
    # Concatenating the last section of the series to the first section of the series, indexed by its beginning, end, and where the cut is occurring
    c(tail(series, length(series) - shift), 
      head(series, shift))
  }