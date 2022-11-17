ConditionSorter <- function(PIDs,
                            RawDir = "/data/Uncertainty/data/raw",
                            BehavDir = "/data/Uncertainty/data/behav/",
                            DerivDir = "/data/Uncertainty/data/deriv/pipeline_1/fmriprep",
                            Tasks = c("3_task-1", "5_task-2"),
                            Suffix)
{
  
  # Creating a For Loop that will Generate Our Three Column Files
  # For each participant listed ...
  for (PID in PIDs){
    
    # ... and for each task they completed
    for (Task in Tasks){
      
      # Import the dataframe containing this participants behavioral correlate
      behav_file <- list.files(path = BehavDir,
                               full.names = F,
                               pattern = paste0("^certainty_neuro_SR-", PID, ".*\\.csv$"))
      
      # If the behavioral correlate file exists
      if (!is_empty(behav_file)){
        
        # and if this is the run that participants actually gave ratings for
        if ((str_detect(behav_file, "condB") & str_detect(Task, "task-2")) | (str_detect(behav_file, "condA") & str_detect(Task, "task-1"))){
          
          # Read in their data
          if (Task == "5_task-2"){
            df <- read.table(paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-uncertainty_run-2", Suffix,"_timing.txt"))
            
          }
          
          if (Task == "3_task-1"){
            df <- read.table(paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-uncertainty_run-1", Suffix,"_timing.txt"))
          }
          
          # Create a new empty column named Condition
          df$Condition <- NA
          
          # If the first value of column 3 is 0
          if (df$V3[1] == 0){
            
            # Set the value of Condition to "NoChange"      
            df$Condition[1] <- "NoChange"
            
          }
          
          # If the first value of column 3 is not 0
          if (df$V3[1] != 0){
            
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
          
          # Define our parametric modualor as 1
          df$paramod <- 1
          
          # Iterate through the various conditions each row could have taken
          for (Cond in c("NoChange", "Increase", "Decrease")){
            df_temp <- df %>%
              subset(.$Condition == Cond) %>%
              subset(select = c("V1", "V2", "paramod"))
            if (Task == "5_task-2"){
              write.table(df_temp,
                          paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-uncertainty_run-2", Suffix,"_Condition-",Cond,"_timing.txt"),
                          sep = "\t",
                          row.names = FALSE,
                          col.names = FALSE)
            }
            if (Task == "3_task-1"){  
              write.table(df_temp,
                          paste0("/data/Uncertainty/data/deriv/pipeline_1/fmriprep/sub-", PID, "/onset/sub-", PID, "_task-uncertainty_run-1", Suffix,"_Condition-",Cond,"_timing.txt"),
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
