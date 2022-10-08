# GenOnsets | v2022.09.03

GenOnsets <- function(PIDs,
                      Tasks = c("3_task-1", "5_task-2"),
                      TR = 2,
                      RawDir = "/data/Uncertainty/data/raw/",
                      BehavDir = "/data/Uncertainty/data/behav/",
                      DerivDir = "/data/Uncertainty/data/deriv/pipeline_1/fmriprep",
                      TrialLevel = T,
                      ParaMod = T){
 
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
  
  if (!dir.exists(BehavDir)){
    stop(paste("BehavDir refers to the parent directory containing the behavioral data to be parametrically regressed onto each participant's BOLD data (one level above the individual participant folders). It must be valid path, but your entry (", 
               BehavDir, ") could not be located. Please enter a new path and try again."))
  }
  
  # Loading packages
  library(tidyverse)
  source("https://github.com/wj-mitchell/neuRotools/blob/main/rucleaner.R?raw=TRUE", local = T)

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
        onset <- seq(107, (nFiles * TR) - 90 - TR, TR)
                
      }
   
      # If the participant's uncertainty task has 729 files
      if (nFiles == 729){
        
        # Create an onset sequence that removes the first 60 and last 60 seconds
        onset <- seq(77, (nFiles * TR) - 60 - TR, TR)
      }
      
      # Create a duration sequence equal to the TR across the board
      duration <- rep(TR, length(onset))
      
      # Import the dataframe containing this participants behavioral correlate
      behav_file <- list.files(path = BehavDir,
                               full.names = F,
                               pattern = paste0("^certainty_neuro_SR-", sprintf("%04d", PID), ".*\\.csv$"))
      
      if (ParaMod == T & ((str_detect(behav_file, "condB") & str_detect(Task, "task-2")) | (str_detect(behav_file, "condA") & str_detect(Task, "task-1")))){
        paramod <- rucleaner(file = behav_file,
                              dir = ReadDir,
                              unit_secs = 2,
                              shave_secs = 17) %>%
                    subset(!str_detect(.$Video, "Control"), select = (CertRate)) %>%
                    abs()
      }

      # Set parametric modulation to the behavioral correlate
      if (ParaMod == F | ((str_detect(behav_file, "condB") & str_detect(Task, "task-1")) | (str_detect(behav_file, "condA") & str_detect(Task, "task-2")))){
        paramod <- rep(1, length(onset))
      }

      # Concatenate onset, duration and parametric modulation into a dataframe
      df_temp <- data.frame(onset, duration, paramod)
      
      # Create a new directory in the participant's raw files called "Onset"
      dir.create(paste0(DerivDir, "/sub-", PID, "/","onset"))
      
      # Set our working directory to that onset directory
      setwd(paste0(DerivDir, "/sub-", PID, "/","onset"))
     
      rows <- seq(1,nrow(df_temp), 60/TR)

      # Iterate through each row in the new dataframe
      for (TRIAL in rows){   

        # If we're working with the first half video ...
        if (Task == "3_task-1"){

          # Save only the target row (which is a single observation) of our dataframe as a text file with this name
          write.table(df_temp[rows[TRIAL]:rows[(TRIAL + 29)],],
                      paste0("sub-", PID, "_task-uncertainty_run-1_min-", TRIAL,"_timing.txt"),
                      sep = "\t",
                      row.names = FALSE,
                      col.names = FALSE)
        }

        # If we're working with the second half video ...
        if (Task == "5_task-2"){

          # Save the target row of our dataframe as a text file with a slightly different name
          write.table(df_temp[rows[TRIAL]:rows[(TRIAL + 29)],],
                      paste0("sub-", PID, "_task-uncertainty_run-2_min-", TRIAL ,"_timing.txt"),
                      sep = "\t",
                      row.names = FALSE,
                      col.names = FALSE)
        }
       }
      
       # Writing an onset file for the spinning checkerboard 
        write.table(data.frame(x=c(seq(30, 
                                       60,
                                       TR),
                                   seq((nFiles * TR) - 60,
                                       (nFiles * TR) - 30, 
                                       TR)), 
                               y=TR,
                               z=1),
                    paste0("sub-", PID, "_task-uncertainty_CB_timing.txt"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
      
      # Cleaning Space
      rm(df_temp, nFiles, onset, paramod, duration)
    }
  }
}
