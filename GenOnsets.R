# GenOnsets | v2022.09.03

GenOnsets <- function(PIDs,
                      Tasks = c("3_task-1", "5_task-2"),
                      TR = 2,
                      RawDir = "/data/Uncertainty/data/raw/",
                      DerivDir = "/data/Uncertainty/data/deriv/pipeline_1/fmriprep"){
 
  # QA Checks
  # Checking TR
  if (TR <= 0 | !is.numeric(TR) | is.na((TR))){
    stop(paste("TR refers to the time of repetition for your scans. It must be a numeric value greater than zero. You've entered", 
               TR, ". Please enter a new value and try again."))
  }
  
  # Checking Directories
  if (!dir.exists(RawDir)){
    stop(paste("RawDir refers to the parent directory containing the raw DICOM files of each participant (one level above the individual participant folders). It must be valid path, but your entry (", 
               RawDir, ") could not be located. Please enter a new path and try again."))
  }
  
  if (!dir.exists(DerivDir)){
    stop(paste("DerivDir refers to the parent directory containing the preprocessed NiFTi files of each participant (one level above the individual participant folders). It must be valid path, but your entry (", 
               DerivDir, ") could not be located. Please enter a new path and try again."))
  }
  
  # Loading packages
  library(tidyverse)

  # Creating a For Loop that will Generate Our Three Column Files
  # For each participant listed ...
  for (PID in PIDs){
  
    # Checking whether filepaths exist for participant
    if (!dir.exists(paste0(RawDir, "/sub-", PID))){
      stop(paste("No DICOM directory could be found for participant", PID, ". Perhaps their data has not been downloaded yet. Please check your filepaths or PIDs and try again." ))
    }
    if (!dir.exists(paste0(DerivDir, "/sub-", PID))){
      stop(paste("No NiFTi directory could be found for participant", PID, ". Perhaps their data has not been pre-processed yet. Please check your filepaths or PIDs and try again." ))
    }
      
    # ... and for each task they completed
    for (Task in Tasks){
      
      # Calculate how many files they have in their raw directory
      # If we find any files at this path, proceed ...
      if (length(list.files(paste0(RawDir, "/sub-", PID, "/",  Task, "/DICOM/"))) != 0){
        nFiles <- length(list.files(paste0(RawDir, "/sub-", PID, "/",  Task, "/DICOM/")))
      }
      
      # ... but if that directory doesn't work, try this other one.
      if (length(list.files(paste0(RawDir, "/sub-", PID, "/",  Task, "/DICOM/"))) == 0){
        nFiles <- length(list.files(paste0(RawDir, "/sub-", PID, "/scans/",  Task, "/DICOM/")))
      }
      
      # If the participant's uncertainty task has 759 files 
      if (nFiles == 759){
        
        # Create an onset sequence that removes the first 90 and last 90 seconds 
        onset <- seq(107, (nFiles * TR) - 90, 60)
      }
      
      # If the participant's uncertainty task has 729 files
      if (nFiles == 729){
        
        # Create an onset sequence that removes the first 60 and last 60 seconds
        onset <- seq(77, (nFiles * TR) - 60, 60)
      }
      
      # Create a duration sequence of 60 across the board
      duration <- rep(60, length(onset))
      
      # Set parametric modulation to 1 across the board
      paramod <- rep(1, length(onset))
      
      # Concatenate onset, duration and parametric modulation into a dataframe
      df_temp <- data.frame(onset, duration, paramod)
      
      # Create a new directory in the participant's raw files called "Onset"
      dir.create(paste0(DerivDir, "/sub-", PID, "/","onset"))
      
      # Set our working directory to that onset directory
      setwd(paste0(DerivDir, "/sub-", PID, "/","onset"))
      
      # If we're working with the first half video ...
      if (Task == "3_task-1"){
        
        # Save our dataframe as a text file with this name
        write.table(df_temp,
                    paste0("sub-", PID, "_task-uncertainty_run-1_timing.txt"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
      }
      
      # If we're working with the second half video ...
      if (Task == "5_task-2"){
        
        # Save our dataframe as a text file with a slightly different name
        write.table(df_temp,
                    paste0("sub-", PID, "_task-uncertainty_run-2_timing.txt"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
      }
      
      # Cleaning Space
      rm(df_temp, nFiles, onset, paramod, duration)
    }
  }
}
