ConditionSorter <- function(PIDs, # An array of participant IDs to Process
                            Tasks = c("3_task-1", "5_task-2"), # An array of the different run name(s) that appear on the DICOM files
                            RawDir = "/data/Uncertainty/data/raw", # The directory in which your DICOM files are stored
                            BehavDir = "/data/Uncertainty/data/behav/", # The directory in which your MRI behavioral data is stored
                            DerivDir = "/data/Uncertainty/data/deriv/pipeline_1/fmriprep", # The directory in which your pre-processed data is stored
                            ZeroValue = 0, # The value that is equivalent to zero or non-movement; could be different from zero if z-scoring is used
                            ParaMod = T, # Whether you'd like to use behavioral data as a parametric modulator
                            Components = c("Control", "Test"), # The study component we'd like to export as the parametric modulator
                            Suffix # A suffix to add to your onset files to better differentiate them from one another
){
  
  # Creating a For Loop that will Generate Our Three Column Files
  # For each participant listed ...
  for (PID in PIDs){
    
    # ... and for each task they completed
    for (TASK in Tasks){
      
      # and for each component of the study
      for (COMPONENT in Components){
      
        # Import the dataframe containing this participants behavioral correlate
        behav_file <- list.files(path = BehavDir,
                                 full.names = F,
                                 pattern = paste0("^certainty_neuro_SR-", PID, ".*\\.csv$"))
        
        # If the behavioral correlate file exists
        if (!is_empty(behav_file)){
          
          # and if this is the run that participants actually gave ratings for
          if ((str_detect(behav_file, "condB") & str_detect(TASK, "task-2")) | (str_detect(behav_file, "condA") & str_detect(TASK, "task-1"))){
            
            # If we're looking at a test component
            if (COMPONENT == "Test"){
            
              # Read in their data
              if (TASK == "5_task-2"){
                df <- read.table(paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-run-2_", Suffix,"_timing.txt"))
              }
              
              if (TASK == "3_task-1"){
                df <- read.table(paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-run-1_", Suffix,"_timing.txt"))
              }
            }
            
            # If we're looking at a control component
            if (COMPONENT == "Control"){
              
              # Read in their data
              df <- read.table(paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-control_", Suffix,"_timing.txt"))
            }
            
            # Create a new empty column named Condition
            df$Condition <- NA
            
            # If the first value of column 3 is 0
            if (df$V3[1] == ZeroValue){
              
              # Set the value of Condition to "NoChange"      
              df$Condition[1] <- "NoChange"
            }
            
            # If the first value of column 3 is not 0
            if (df$V3[1] != ZeroValue){
              
              # Set the value of Condition to "Increase"
              df$Condition[1] <- "Increase"
            }
            
            # Iterate through each subsequent row in the dataframe        
            for (row in 2:nrow(df)){
              
              # If the value of our given row is the same as the previous
              if (df$V3[row] == df$V3[row - 1]){
                
                # Set the value of that row as "NoChange"    
                df$Condition[row] <- "NoChange"
              }
              
              # If the value of our given row is less than the previous
              if (df$V3[row] < df$V3[row - 1]){
                
                # Set the value of that row as "Decrease"    
                df$Condition[row] <- "Decrease"
              }
              
              # If the value of our given row is more than the previous
              if (df$V3[row] > df$V3[row - 1]){
                
                # Set the value of that row as "Increase"    
                df$Condition[row] <- "Increase"
              }
            }
            
            # If we want to include parametric modulators
            if (ParaMod == TRUE){
              df$paramod <- df$V3
            }
            
            # If we don't want to include parametric modulators
            if (ParaMod == FALSE){
              # Define our parametric modualtor as 1
              df$paramod <- 1
            }
            
            # Iterate through the various conditions each row could have taken
            for (Cond in c("NoChange", "Increase", "Decrease")){
                
              # Create a temporary data frame
              df_temp <- df %>%
                
                # Fill it with observations that meet our currently-iterated condition
                subset(.$Condition == Cond) %>%
                
                # Only use the onset, duration, and parametric modulator columns
                subset(select = c("V1", "V2", "paramod"))
              
              # If we're looking at a test component
              if (COMPONENT == "Test"){
                
                # If we're looking at a second-half video rater
                if (TASK == "5_task-2"){
                  
                  # Save the file under this filename
                  write.table(df_temp,
                              paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-run-2_", Suffix,"_condition-",Cond,"_timing.txt"),
                              sep = "\t",
                              row.names = FALSE,
                              col.names = FALSE)
                }
                
                # If we're looking at a first-half video rater
                if (TASK == "3_task-1"){  
                  
                  # Save the file under this filename
                  write.table(df_temp,
                              paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-run-1_", Suffix,"_condition-",Cond,"_timing.txt"),
                              sep = "\t",
                              row.names = FALSE,
                              col.names = FALSE)
                }
              }
              
              # If we're looking at a control component
              if (COMPONENT == "Control"){
                
                # Save the file under this filename
                write.table(df_temp,
                            paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-control_", Suffix,"_condition-",Cond,"_timing.txt"),
                            sep = "\t",
                            row.names = FALSE,
                            col.names = FALSE)
              }
            }
          }
        }
      }
    }
  }
}
