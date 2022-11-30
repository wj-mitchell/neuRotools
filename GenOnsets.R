GenOnsets <- function(PIDs,  # An array of participant IDs to Process
                      Tasks = c("3_task-1", "5_task-2"), # An array of the different run name(s) that appear on the DICOM files
                      RawDir = "/data/Uncertainty/data/raw/", # The directory in which your DICOM files are stored
                      BehavDir = "/data/Uncertainty/data/behav/", # The directory in which your MRI behavioral data is stored
                      DerivDir = "/data/Uncertainty/data/deriv/pipeline_1/fmriprep/", # The directory in which your preprocessed data is stored
                      TR = 2, # The length of your repetition time in seconds
                      ShaveLength = 17, # How much time should be shaved from the beginning of your data in seconds
                      ShavedFile = F, # Whether we want to generate a separate onset file for the shaved segment
                      Checkerboard = F, # Whether to output onset files for the spinning checkerboard 30s before and after video
                      Suffix = "", # A suffix to add to your onset files to better differentiate them from one another
#                       SeparateFiles = F, # An argument as to whether each trial should be saved as a separate onset file
                      ParaMod = T, # Whether you'd like to use behavioral data as a parametric modulator
                      Components = c("Control", "Test"), # The study component we'd like to export as the parametric modulator
                      Method = c("CPA", "Inflections", "Bins"), # Whether you'd like events in your parametric modulator to be defined by evenly-spaced bins (e.g., a 1200s video could be 20 trials each of 60s length), by any changes in ratings, or by using a PELT method change point analysis
                      BinLength = 30, # If using the Bin Method, the size of each bin in seconds
                      Demean = T, # Whether your parametric modulator should be demeaned (i.e., calculate the average and subtract it from each data point in the time course such that data on either side of the mean are balanced)
                      Detrend = T, # Whether your parametric modulator should be detrended (i.e., subtract the value of the subsequent time point from each time point so that only changes deviate from 0)
                      ZScore = T, # Whether your parametric modulator should be standardized (i.e., converted to standard deviation units to make it more comparable across individuals and studies)
                      BufferBefore = 10, # The time in seconds prior to each inflection point that should take the same value as the inflection point (to be used when the onset/duration of the event is probabilistic or unknown)
                      BufferAfter = 10, # The time in seconds following each inflection point that should take the same value as the inflection point (to be used when the onset/duration of the event is probabilistic or unknown)
                      Smoothing = T, # Whether inflections within the same overlapping bufferzones should be smoothed together (i.e., averaged across all inflection points)
                      Threshold = 1.5, # The value in standard deviation units below which parametric modulator inflections should be ignored
                      OffsetLength = 0, # How many TRs the parametric modulator should be offset by (Negative values will lag a given parametric modulation value behind it's associated trial, and positive value will make a parametric value precede its trial).
                      Override = NA, # Generate a cluster-based parametric modulator, but then override the values (useful to generate mean EVs) 
                      UseConditionSorter = T # Whether to use the custom condition sorter function which will break parametric modulations into three onset types: increases, decreases, or no changes
                      # GammifyBins = 6, # [IN DEVELOPMENT] Identify how many trials prior to the target trail should be affected by the gamma distribution
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
    if (!dir.exists(paste0(RawDir, "sub-", PID))){
      stop(paste("No DICOM directory could be found for participant", PID, ". Perhaps their data has not been downloaded yet. Please check your filepaths or PIDs and try again." ))
    }
    
    if (!dir.exists(paste0(DerivDir, "sub-", PID))){
      stop(paste("No NiFTi directory could be found for participant", PID, ". Perhaps their data has not been pre-processed yet. Please check your filepaths or PIDs and try again." ))
    }
    
    # ... and for each task they completed
    for (TASK in Tasks){
      
      # ... and for each component of the study 
      for (COMPONENT in Components){
    
        ## GENERATING ONSET ----
                
        # If we're working with a test component
        if (COMPONENT == "Test"){
          
          #Specifying a few misceallenous details specific to each task
          TaskFiles <- 519
          OnsetBuffer <- 17
          Folder <- TASK
        }
        
        # If we're working with a control component
        if (COMPONENT == "Control"){
          
          #Specifying a few misceallenous details specific to each task
          TaskFiles <- 0
          OnsetBuffer <- 0
          Folder <- "7_task-3"
        }
       
        # Calculate how many files they have in their raw directory
        # If we find any files at this path, proceed ...
        if (length(list.files(paste0(RawDir, "sub-", PID, "/",  Folder, "/DICOM/"))) != 0){
          nFiles <- length(list.files(paste0(RawDir, "sub-", PID, "/",  Folder, "/DICOM/")))
        }

        # ... but if that directory doesn't work, try this other one.
        if (length(list.files(paste0(RawDir, "sub-", PID, "/scans/",  Folder, "/DICOM/"))) != 0){
          nFiles <- length(list.files(paste0(RawDir, "sub-", PID, "/scans/",  Folder, "/DICOM/")))
        }
        
        # ... and if that didn't work, try this
        if ((length(list.files(paste0(RawDir, "sub-", PID, "/",  Folder, "/DICOM/"))) == 0) & 
            (length(list.files(paste0(RawDir, "sub-", PID, "/scans/",  Folder, "/DICOM/"))) == 0) &
            (length(list.files(paste0(RawDir, "sub-", PID, "/",  Folder, "/"))) != 0)){
          nFiles <- length(list.files(paste0(RawDir, "sub-", PID, "/",  Folder, "/")))
        }

        # If the participant's uncertainty Task has 759 files 
        if (nFiles == 240 + TaskFiles){
          # Create an onset sequence that removes the first 90 and last 90 seconds 
          onset <- seq((90 + OnsetBuffer), ((nFiles * TR) - 90 - TR), TR)
        }

        # If the participant's uncertainty Task has 729 files
        if (nFiles == 210 + TaskFiles){
          # Create an onset sequence that removes the first 60 and last 60 seconds
          onset <- seq((60 + OnsetBuffer), (nFiles * TR) - 60 - TR, TR)
        } 

        ## GENERATING DURATION ----
        # Create a duration sequence equal to the TR across the board
        duration <- rep(TR, length(onset))
        
        ## GENERATING PARAMETRIC MODULATOR ----
        # if we want to create a parametric modulator
        if (ParaMod == TRUE){
        
          # Import the dataframe containing this participants behavioral correlate
          behav_file <- list.files(path = BehavDir,
                                   full.names = F,
                                   pattern = paste0("^certainty_neuro_SR-", PID, ".*\\.csv$"))
          
          # Or just set it to '1'
          if (is_empty(behav_file) | ((str_detect(behav_file, "condB") & str_detect(TASK, "task-1")) | (str_detect(behav_file, "condA") & str_detect(TASK, "task-2")))){
            paramod <- rep(1, length(onset))
            next
          }
          
          # If the behavioral correlate file exists
          if (!is_empty(behav_file)){
            
            # and if this is the run that participants actually gave ratings for
            if ((str_detect(behav_file, "condB") & str_detect(TASK, "task-2")) | (str_detect(behav_file, "condA") & str_detect(TASK, "task-1"))){
  
              # Turn the time course in that behavioral file into parametrid modulator
              paramod <- rucleaner(file = behav_file,
                                   dir = BehavDir,
                                   unit_secs = TR,
                                   shave_secs = ShaveLength) %>%
                
                # Remove observations for the control video
                subset(str_detect(.$Video, COMPONENT), select = (CertRate)) %>%
                
                # Get the absolute value of certainty 
                abs() %>% 
                
                # And reformat the values to an array of numeric values
                unlist %>%
                as.numeric()
          
              # If we want to detrend the data ...
              if (Detrend == TRUE){
                
                # Use the difference function to basically subtract from any given datapoint the value that preceded it
                paramod <- c(0, diff(paramod))
              }
              
              # We may later transform this array and non-zero values will be hard to distinguish from zero values. 
              # Therefore, we are identifying a zero value now for reference
              zero_point <- which(paramod == 0)[1]
              
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
              
              # Any data points less than threshold will be converted to a value equal to zero
              if (Threshold > 0){
                
                # First I'm removing those timepoints from the inflections variable, if they exist there
                for (CENSORED in which(abs(paramod) < Threshold)){
                  if (any(inflections == CENSORED)){
                    inflections <- inflections[-which(inflections == CENSORED)]
                  }
                }
                
                # Now we'll actually censor those timepoints
                paramod[which(abs(paramod) < Threshold)] <- paramod[zero_point]
              }
              
              # Our zero value may be vastly different from what it was now, so let's define it more concretely
              zero_value <- paramod[zero_point]
              
              # If we want to buffer inflections by some amount
              if (BufferBefore > 0 | BufferAfter > 0){
                
                # and if we want to smooth across contiguous time points which might otherwise compete
                if (Smoothing == TRUE){
                  
                  # Specifying our starting cluster
                  cluster_num <- 1  
                  
                  # Creating an array to store the clusters each inflection belongs to
                  clusters <- rep(NA, length(inflections))
                  
                  # Iterate through each of the identified inflection points
                  for (CHANGE in inflections){
                    
                    # If a given inflection point stands alone or has neighbors
                    if (length(which(inflections > CHANGE - BufferBefore/TR & inflections < CHANGE + BufferAfter/TR)) >= 1){
                      
                      # Label those contiguous or isolated points as a cluster  
                      clusters[which(inflections > CHANGE - BufferBefore/TR & inflections < CHANGE + BufferAfter/TR)] <- paste0("Cluster", cluster_num)
                    }
                    
                    # If we find an inflection point which doesn't have any contiguous points
                    if (length(which(inflections > CHANGE & inflections < CHANGE + BufferAfter/TR)) == 0){
                      
                      # Increment the cluster number by 1
                      cluster_num <- cluster_num + 1
                    }
                  }
                  
                  # Cleaning Our Space
                  rm(cluster_num)
                  
                  # Iterating through each of the clusters
                  for (CLUSTER in unique(clusters)){
                    
                    # And smoothing values within cluster by averaging across all values
                    paramod[inflections[which(clusters == CLUSTER)]] <- mean(paramod[inflections[which(clusters == CLUSTER)]])
                  }
                  
                  # Iterating through each of those inflections again
                  for (INFLECTION in inflections){
                    
                    # Specifying our target region
                    target <- (INFLECTION - BufferBefore/TR):(INFLECTION + BufferAfter/TR)
                    
                    # Copying the value of each inflection point to the buffered time points before and afterwards
                    paramod[target[target > 0 & target < length(paramod)]] <- paramod[INFLECTION]
                  }
                  
                  # Cleaning Our Space
                  rm(target, clusters)
                }
                
                # and if we do not want to smooth across contiguous time points
                if (Smoothing == FALSE){
                  
                  # Creating a temporary array so that adding value to overlapping inflections/buffer zones doesn't compound their effects
                  paramod_temp <- paramod
                  
                  # Iterating through each of those inflections again
                  for (INFLECTION in inflections){
                    
                    # Iterating through each of the time points within the buffer zone for that inflection
                    for (BUFFERPOINT in (INFLECTION - BufferBefore/TR):(INFLECTION + BufferAfter/TR)){
                    
                      # Excluding iterations where inflection and bufferpoint are synonymous
                      if (BUFFERPOINT != INFLECTION & BUFFERPOINT > 0 & BUFFERPOINT < length(paramod)){
                      
                        # Adding the value of each inflection point to the buffered time points before and afterwards
                        paramod_temp[BUFFERPOINT] <- paramod_temp[BUFFERPOINT] + (paramod[INFLECTION] - paramod[zero_point])  
                      }
                    }
                  }
                    
                  # Cleaning Our Space
                  paramod <- paramod_temp
                  rm(paramod_temp)
                }
                
                # Before we finish up, we may have just paved over our zero point, so let's redefine it to be sure
                zero_point <- which(paramod == zero_value)[1]
                
                # Rescaling our array through demeaning, standardizing, or both (or neither!)
                paramod <- as.numeric(scale(paramod, scale = ZScore, center = Demean))
                
                # We likely just redefined our zero value again, though, so one more time with the new zero point
                zero_value <- paramod[zero_point]
              }
    
              # Offset ratings by a certain number of trials 
              if (OffsetLength != 0){
                paramod <- paramod[(1 + OffsetLength):(length(paramod) + OffsetLength)]
                paramod[is.na(paramod)] <- 0 
              }
          }
          }
        }
        
        # If we just want a standard parametric modulator just set it to '1'
        if (ParaMod == FALSE){
          paramod <- rep(1, length(onset))
        }
        
        if (Method == "CPA" | Method == "Inflections"){
          
          if ((str_detect(behav_file, "condB") & str_detect(TASK, "task-2")) | (str_detect(behav_file, "condA") & str_detect(TASK, "task-1"))){
           
            # Specifying our starting cluster
            cluster_num <- 1  
            
            # Creating an array to store the clusters each collection of points belongs to
            clusters <- rep(NA, length(inflections))
            
            # Iterate through each of the identified inflection points
            for (DATAPOINT in 1:length(paramod)){
              
              # If we're on the first datapoint
              if (DATAPOINT == 1){
                
                # Specify that that observation is in the first cluster
                clusters[DATAPOINT] <- paste0("Cluster", cluster_num)
              }
              
              # If we're past the first observation
              if (DATAPOINT > 1){
                
                # and the current iteration is different than the one that preceded it
                if (paramod[DATAPOINT] != paramod[DATAPOINT - 1]){
                  
                  # Increment the cluster number by 1
                  cluster_num <- cluster_num + 1
                  
                  # Label that observation as belonging to a new cluster
                  clusters[DATAPOINT] <- paste0("Cluster", cluster_num)
                }
                
                # and the current iteration is the same as the one that preceded it
                if (paramod[DATAPOINT] == paramod[DATAPOINT - 1]){
                  
                  # Label that observation as belonging to the same cluster
                  clusters[DATAPOINT] <- paste0("Cluster", cluster_num)
                }
              }
            }
            
            # Cleaning Our Space
            rm(cluster_num)
            
            # Creating empty arrays
            paramod_cluster <- rep(NA, length(unique(clusters)))
            onset_cluster <- rep(NA, length(unique(clusters)))
            duration_cluster <- rep(NA, length(unique(clusters)))
            
            # Iterating through each of the clusters
            for (CLUSTER in 1:length(unique(clusters))){
              
              # Identifying the duration of each cluster by counting how many observations each has and multiplying it by TR
              duration_cluster[CLUSTER] <- length(which(clusters == unique(clusters)[CLUSTER])) * TR
              
              # Identifying the onset of each cluster as the first onset associated with that cluster
              onset_cluster[CLUSTER] <- onset[which(clusters == unique(clusters)[CLUSTER])][1]
              
              # Identifying the parametric modulator of each cluster as the first value associated with that cluster
              paramod_cluster[CLUSTER] <- paramod[which(clusters == unique(clusters)[CLUSTER])][1]
              
            }
            
            # Concatenate onset, duration and parametric modulation into a data frame
            df_temp <- data.frame(onset_cluster, duration_cluster, paramod_cluster) 
            
            # Cleaning Our Space
            rm(paramod, paramod_cluster, onset, onset_cluster, duration, duration_cluster) 
          }
        }
        
        # If we want even length bins of observations
        if (Method == "Bins" | ((str_detect(behav_file, "condB") & str_detect(TASK, "task-1")) | (str_detect(behav_file, "condA") & str_detect(TASK, "task-2")))){
          
          # Create a sequence of onsets from the original onset variable spaced apart according to how large our bins are
          onset_bin <- seq(onset[1], onset[length(onset)], BinLength)
          
          # Make the duration of each onset equal to the bin length
          duration_bin <- rep(BinLength, length(onset_bin))
          
          # Create an empty parametric modulator to use in this for loop
          paramod_bin <- rep(NA, length(onset_bin))
          
          # in which we iterate through each onset
          for (ONSET in 1:length(onset_bin)){
            
            # If we're not on the last iteration
            if (ONSET != length(onset_bin)){
              
              # And calculate the parametric modulator as the average of observations across the bin. 
              paramod_bin[ONSET] <- mean(paramod[which(onset == onset_bin[ONSET]):(which(onset == onset_bin[ONSET + 1]) - 1)])  
            }
            
            # If we're on the last iteration
            if (ONSET == length(onset_bin)){
              
              # And calculate the parametric modulator as the average of observations across the bin. 
              paramod_bin[ONSET] <- mean(paramod[which(onset == onset_bin[ONSET]):(which(onset == onset_bin[ONSET]) + (BinLength/TR) - 1)])  
            }
          }
          
          # Concatenate onset, duration and parametric modulation into a data frame
          df_temp <- data.frame(onset_bin, duration_bin, paramod_bin)
          
          # Cleaning Our Space
          rm(paramod, paramod_bin, onset, onset_bin, duration, duration_bin)
        }
        
        # If we want to override the values
        if (!is.na(Override) & is.numeric(Override)){
          df_temp[,grep(x = names(df_temp), 
                        pattern = "paramod")] <- Override
        }
        
        # If an onset directory doesn't already exist
        if (!dir.exists(paste0(DerivDir, "sub-", PID, "/","onset"))){
          # Create a new directory in the participant's raw files called "Onset"
          dir.create(paste0(DerivDir, "sub-", PID, "/","onset"))
        }
        
        # Set our working directory to that onset directory
        setwd(paste0(DerivDir, "sub-", PID, "/","onset"))
        
        # If we're looking at a test component
        if (COMPONENT == "Test"){{
    #       if (SeparateFiles == T){
    #         #Tracking which rows denote the start of a new trial 
    #         rows <- seq(1, nrow(df_temp), nrow(df_temp) / )
    #       }
          
    #       if (SeparateFiles == F){
    #         #Tracking which rows denote the start of a new trial 
    #         rows <- 1
    #       }
          
          # Iterate through each row that starts a new trial in the new dataframe
    #       for (TRIAL in 1:length(rows)){   
          }
          
          # If we're working with the first half video ...
          if (TASK == "3_task-1"){{
    #           if (length(rows) != 1){
    #             # Save only the target row (which is a single observation) of our dataframe as a text file with this name
    #             write.table(df_temp[rows[TRIAL]:(rows[TRIAL] + ((nrow(df_temp) / Trial_Num) - 1)),],
    #                         paste0("sub-", PID, "_task-uncertainty_run-1_min-", TRIAL, Suffix ,"_timing.txt"),
    #                         sep = "\t",
    #                         row.names = FALSE,
    #                         col.names = FALSE)
    #           }
              
    #           if (length(rows) == 1){
            }
            
            # Save only the target row (which is a single observation) of our dataframe as a text file with this name
            write.table(df_temp,
                        paste0("sub-", PID, "_task-run-1_ParaMod-", ParaMod ,"_Override-", Override, "_Method-", Method, "_Buffer-", BufferBefore, 
                               "s_Smoothing-", Smoothing, "_Threshold-", Threshold, "sd_Offset-", OffsetLength, "s_", Suffix ,"_timing.txt"),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE)
          }

          # If we're working with the second half video ...
          if (TASK == "5_task-2"){{
    #           if (length(rows) != 1){
                
    #             # Save the target row of our dataframe as a text file with a slightly different name
    #             write.table(df_temp[rows[TRIAL]:(rows[TRIAL] + ((nrow(df_temp) / Trial_Num
    #             ) - 1)),],
    #             paste0("sub-", PID, "_task-uncertainty_run-2_min-", TRIAL , Suffix ,"_timing.txt"),
    #             sep = "\t",
    #             row.names = FALSE,
    #             col.names = FALSE)
    #           }
              
    #           if (length(rows) == 1){
          }
            
            # Save the target row of our dataframe as a text file with a slightly different name
            write.table(df_temp,
                        paste0("sub-", PID, "_task-run-2_ParaMod-", ParaMod ,"_Override-", Override, "_Method-", Method, "_Buffer-", BufferBefore, 
                               "s_Smoothing-", Smoothing, "_Threshold-", Threshold, "sd_Offset-", OffsetLength, "s_", Suffix ,"_timing.txt"),
                        sep = "\t",
                        row.names = FALSE,
                        col.names = FALSE)
          }
        }
        
        # If we're looking at a test component
        if (COMPONENT == "Control"){
          
          # Save the target row of our dataframe as a text file with a slightly different name
          write.table(df_temp,
                      paste0("sub-", PID, "_task-control_ParaMod-", ParaMod ,"_Override-", Override, "_Method-", Method, "_Buffer-", BufferBefore, 
                               "s_Smoothing-", Smoothing, "_Threshold-", Threshold, "sd_Offset-", OffsetLength, "s_", Suffix ,"_timing.txt"),
                      sep = "\t",
                      row.names = FALSE,
                      col.names = FALSE)
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
                      paste0("sub-", PID, "_task-CB_", Suffix ,"_timing.txt"),
                      sep = "\t",
                      row.names = FALSE,
                      col.names = FALSE)
        }
        
        if (ShavedFile == TRUE & ShaveLength > 0){
          
          # Writing an onset file for the shaved data
          write.table(data.frame(x=60, 
                                 y=ShaveLength,
                                 z=1),
                      paste0("sub-", PID, "_shaved_", Suffix ,"_timing.txt"),
                      sep = "\t",
                      row.names = FALSE,
                      col.names = FALSE)
        }
        
        # Cleaning Space
        rm(df_temp, nFiles)
        
        # If we want to take the extra step to use ConditionSorter
        if (UseConditionSorter == TRUE){
  
          # Source the function
          source("https://github.com/wj-mitchell/neuRotools/blob/main/ConditionSorter.R?raw=TRUE", local = T)
  
          # We can't have two arguments of the same name with nested functions, so I'm creating a temporary one
          suffix <- paste0("ParaMod-", ParaMod ,"_Override-", Override, "_Method-", Method, "_Buffer-", BufferBefore, 
                               "s_Smoothing-", Smoothing, "_Threshold-", Threshold, "sd_Offset-", OffsetLength, "s_", Suffix)
  
          # Use Condition Sorter
          ConditionSorter(PIDs = PID,
                          Tasks = TASK,
                          ZeroValue = zero_value,
                          Suffix = suffix,
                          Components = COMPONENT)
        }
      }
    }
  }
}
