GenOnsets <- function(PIDs,  # An array of participant IDs to Process
                      Tasks = c("3_task-1", "5_task-2"), # An array of the task names that appear on the DICOM files
                      RawDir = "/data/Uncertainty/data/raw", # The directory in which your DICOM files are stored
                      BehavDir = "/data/Uncertainty/data/behav/", # The directory in which your MRI behavioral data is stored
                      DerivDir = "/data/Uncertainty/data/deriv/pipeline_1/fmriprep", # The directory in which your preprocessed data is stored
                      TR = 2, # The length of your repetition time in seconds
                      Shave_Length = 0, # How much time should be shaved from the beginning of your data
                      Checkerboard = F, # Whether to output onset files for the spinning checkerboard 30s before and after video
                      Suffix = "", # A suffix to add to your onset files to better differentiate them from one another
                      SeparateFiles = F, # An argument as to whether each trial should be saved as a separate onset file
                      ParaMod = T, # Whether you'd like to use behavioral data as a parametric modulator
                      Method = c("CPA", "Inflections", "Bins"),
                      Demean = T,
                      Detrend = F, # Whether that parametric modulation should specifically identify changes in behavior by subtracting the previous trials value from the current trial
                      ZScore = T,
                      Threshold = 1000,
                      OffsetLength = NA, # How many trials the parametric modulator should be offset by. Negative values will lag a given parametric modulation value behind it's associated trial, and positive value will make a parametric value precede its trial.
                      # GammifyBins = 6, # [IN DEVELOPMENT] Identify how many trials prior to the target trail should be affected by the gamma distribution
                      BufferBefore = 0, # [IN DEVELOPMENT] How long before the inflection the event should include in seconds
                      BufferAfter = 0, # [IN DEVELOPMENT] How long after the inflection the event should include in seconds
                      # BinLength = 60, # The length of each of your trials in seconds
                      # BinNum = 22, # The number of trials you have
                      # Stim_Length # [IN DEVELOPMENT] An argument to use in lieu of specifying how many trials we have. It will be divided by the Trial_Length Argument and automatically calcluate the number of trials you have
                      
                      
){
  
  # QA Checks [IN DEVELOPMENT] 
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
  if (require("pacman") == FALSE){
    install.packages("pacman")
  }
  
  # Loading in my packages with my pacman manager
  pacman::p_load(changepoint,
                 tidyverse)
  
  # Loading custom functions
  source("https://github.com/wj-mitchell/neuRotools/blob/main/rucleaner.R?raw=TRUE", local = T)
  # source("https://github.com/wj-mitchell/neuRotools/blob/main/gammify.R?raw=TRUE", local = T)  
  
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
      
      # if we want to create a parametric modulator
      if (ParaMod == T){
      
        # Import the dataframe containing this participants behavioral correlate
        behav_file <- list.files(path = BehavDir,
                                 full.names = F,
                                 pattern = paste0("^certainty_neuro_SR-", PID, ".*\\.csv$"))
        
        # if we don't have the behavioral data
        if (is_empty(behav_file)){
          
          # Move on to the next iteration
          next
        }
      
        # If the behavioral correlate file exists
        if (!is_empty(behav_file)){
          
          # and if this is the run that participants actually gave ratings for
          if ((str_detect(behav_file, "condB") & str_detect(Task, "task-2")) | (str_detect(behav_file, "condA") & str_detect(Task, "task-1"))){

            # Turn the time course in that behavioral file into parametrid modulator
            paramod <- rucleaner(file = behav_file,
                                 dir = BehavDir,
                                 unit_secs = Trial_Length,
                                 shave_secs = Shave_Length) %>%
              
              # Remove observations for the control video
              subset(!str_detect(.$Video, "Control"), select = (CertRate)) %>%
              
              # Get the absolute value of certainty 
              abs() %>% 
              
              # And reformat the values to an array of numeric values
              unlist %>%
              as.numeric()
            
              # We may later transform this array and non-zero values will be hard to distinguish from zero values. 
              # Therefore, we are identifying a zero value now for reference
              zero_point <- which(paramod == 0)[1]
          }
          
          # If we want to detrend the data ...
          if (Detrend == TRUE){
            
            # Use the difference function to basically subtract from any given datapoint the value that preceded it
            paramod <- c(0, diff(paramod))
          }
          
          # If we want to use change point analysis ...
          if (Method == "CPA"){
            
            # We'll identify inflections using the cpts command
            # NOTE: We MUST use CPA before scaling
            inflections <- cpts(cpt.mean(paramod, method = "PELT"))
          }
          
          # But if we want to use all inflections 
          if (Method == "Inflections"){
            
            # We'll identify any array value that's non-zero
            inflections <- which(paramod != paramod[zero_point])
          }
          
          # Scaling our array through demeaning, standardizing, or both (or neither!)
          paramod <- as.numeric(scale(paramod, scale = ZScore, center = Demean))
          
          # Any datapoints less than threshold will be converted to a value equal to zero
          # First I'm removing those timepoints from the inflections variable, if they exist there
          for (CENSORED in which(paramod < Threshold)){
            if (any(inflections == CENSORED)){
              inflections <- inflections[-which(inflections == CENSORED)]
            }
          }
          
          # Now we'll actually censor those timepoints
          paramod[paramod < Threshold] <- paramod[zero_point]
          
          
            
            # Use inflections as events. Take the value of any given inflection and add it to all timepoints around that time point
            if (Inflection == TRUE){
              for (ITEM in 1:length(paramod)){
                if (paramod[ITEM] != zerovalue){
                  for (TRIAL in (ITEM - (Inflection_Before/TR)):(ITEM + (Inflection_After/TR))){
                    if (TRIAL != ITEM & TRIAL > 0 & TRIAL <= length(paramod)){
                      paramod[TRIAL] <- paramod[ITEM] + paramod[TRIAL]
                    }
                  }
                }
              }
            }
          }
          
          rm(zerovalue, zeropoint)
          
          # Offset
          
          
          # Offset ratings by a certain number of trials 
          if (ParaMod_Offset == TRUE){
            paramod <- paramod[(1 + ParaMod_OffsetLength):(length(paramod)  + ParaMod_OffsetLength)]
            paramod[is.na(paramod)] <- 0 
          }
        }
          
          # But if this is another run ....
          if ((str_detect(behav_file, "condB") & str_detect(Task, "task-1")) | (str_detect(behav_file, "condA") & str_detect(Task, "task-2"))){
            # just set it to '1'
            paramod <- rep(1, length(onset))
          }
        }
      
        # Also, if we don't want to use the behavioral correalte as the parametric modulator
        if (ParaMod == F){
          # just set it to '1'
          paramod <- rep(1, length(onset))
        }
      }
      
      
      
      
      
      
      
      
      
      
      
      
      
      if (Method == "CPA" | Method == "Inflections"){
        
      }
      
      # MUST BE ALTERED FOR CPA APPROACH
      {
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
          onset <- seq(107, (nFiles * TR) - 90 - Trial_Length, Trial_Length)
        }
        
        # If the participant's uncertainty task has 729 files
        if (nFiles == 729){
          # Create an onset sequence that removes the first 60 and last 60 seconds
          onset <- seq(77, (nFiles * TR) - 60 - Trial_Length, Trial_Length)
        }
        
        # Create a duration sequence equal to the TR across the board
        duration <- rep(Trial_Length, length(onset))
        
      }
      
      # Or just set it to '1'
      if (is_empty(behav_file) | ParaMod == F | ((str_detect(behav_file, "condB") & str_detect(Task, "task-1")) | (str_detect(behav_file, "condA") & str_detect(Task, "task-2")))){
        paramod <- rep(1, length(onset))
      }
      
      # Concatenate onset, duration and parametric modulation into a dataframe
      df_temp <- data.frame(onset, duration, paramod)
      
      # If an onset directory doesn't already exist
      if (!dir.exists(paste0(DerivDir, "/sub-", PID, "/","onset"))){
        # Create a new directory in the participant's raw files called "Onset"
        dir.create(paste0(DerivDir, "/sub-", PID, "/","onset"))
      }
      
      # Set our working directory to that onset directory
      setwd(paste0(DerivDir, "/sub-", PID, "/","onset"))
      
      if (SeparateFiles == T){
        #Tracking which rows denote the start of a new trial 
        rows <- seq(1, nrow(df_temp), nrow(df_temp) / Trial_Num)
      }
      
      if (SeparateFiles == F){
        #Tracking which rows denote the start of a new trial 
        rows <- 1
      }
      
      # Iterate through each row that starts a new trial in the new dataframe
      for (TRIAL in 1:length(rows)){   
        # If we're working with the first half video ...
        if (Task == "3_task-1"){
          if (length(rows) != 1){
            # Save only the target row (which is a single observation) of our dataframe as a text file with this name
            write.table(df_temp[rows[TRIAL]:(rows[TRIAL] + ((nrow(df_temp) / Trial_Num) - 1)),],
                        paste0("sub-", PID, "_task-uncertainty_run-1_min-", TRIAL, Suffix ,"_timing.txt"),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE)
          }
          
          if (length(rows) == 1){
            # Save only the target row (which is a single observation) of our dataframe as a text file with this name
            write.table(df_temp,
                        paste0("sub-", PID, "_task-uncertainty_run-1", Suffix ,"_timing.txt"),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE)
          }
        }
        
        # If we're working with the second half video ...
        if (Task == "5_task-2"){
          if (length(rows) != 1){
            # Save the target row of our dataframe as a text file with a slightly different name
            write.table(df_temp[rows[TRIAL]:(rows[TRIAL] + ((nrow(df_temp) / Trial_Num
            ) - 1)),],
            paste0("sub-", PID, "_task-uncertainty_run-2_min-", TRIAL , Suffix ,"_timing.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
          }
          if (length(rows) == 1){
            # Save the target row of our dataframe as a text file with a slightly different name
            write.table(df_temp,
                        paste0("sub-", PID, "_task-uncertainty_run-2", Suffix ,"_timing.txt"),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE)
          }
        }
      }
      
      if (Checkerboard == T){
        # Writing an onset file for the spinning checkerboard 
        write.table(data.frame(x=c(seq(30, 
                                       60 - TR,
                                       TR),
                                   seq((((nFiles - 1) * TR) + 1) - 60,
                                       (((nFiles - 1) * TR) + 1) - 30 - TR, 
                                       TR)), 
                               y=TR,
                               z=1),
                    paste0("sub-", PID, "_task-uncertainty_CB", Suffix ,"_timing.txt"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
      }
      
      if (Shave_Length > 0){
        # Writing an onset file for the shaved data
        write.table(data.frame(x=60, 
                               y=Shave_Length,
                               z=1),
                    paste0("sub-", PID, "_task-uncertainty_Shaved", Suffix ,"_timing.txt"),
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
      }
      # Cleaning Space
      rm(df_temp, nFiles, onset, paramod, duration)
    }
  }
}
