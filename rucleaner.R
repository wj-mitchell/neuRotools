## rucleaner.R | v2024.05.06
  
  # This function will reorganize any output file from the Regulating Uncertainty
  # neuro project in a format more conducive to analysis. It will automatically identify 
  # the source script the data came from and change the analysis as needed. 
  rucleaner <- function(file, # String containing the filename that must be cleaned
                      dir, # String containing the working directory the file is contained within
                      unit_secs = 2, # A numeric value denoting the interval of time, in seconds, that output should be displayed in.
                      # Data was collected frame-by-frame. A value of NA would produce unaveraged, raw output. A value of 1
                      # would yield data averaged on a second by second basis. A value of 60 would yield data averaged on a 
                      # minute by minute basis, etc.
                      shave_secs = 17) # A numeric value denoting the interval of time, in seconds, that should be ignored/removed from the 
  # beginning of data collection. So, for example, the first few seconds of naturalistic stimuli, if not
  # captured using a checkerboard first, should probably be ignored in fMRI research. 
{
  # Setup ----
  # Loading packages
  if (require("pacman") == FALSE){
    install.packages("pacman")
  }
  
  ## Package Loading ----
  pacman::p_load(assertthat, tidyverse)
  
  ## Options ----
  options(scipen=100)
  options(digits=3)
  options(contrasts = c("contr.helmert", "contr.poly"))
  
  ## Errors ----
  if (!is.string(file) | !is.string(dir)){
    stop(paste("Both file and dir must be entered as strings. Please update file and dir to comply.", sep = " "))
  }
  if (!is.na(unit_secs)){
    if (!is.numeric(unit_secs) | unit_secs > 300){
      stop(paste("unit_secs must be set to either a numeric value below 300 or NA. You have entered:", unit_secs, sep = " "))
    }
  }
  if (!is.numeric(shave_secs)){
    stop(paste("shave_secs must be set to a numeric value; either 0, if you do not wish to remove any observations from the beginning of your data, or a value greater than 0. You have entered:", shave_secs, sep = " "))
  }
  if (!str_detect(file,"\\.csv")){
    stop(paste("Submitted files must be in .csv format. The file you listed [", file, "] does not have a .csv extension. Please add .csv to it.", sep = " "))
  }
  if (!file.exists(paste0(dir,file))){
    stop(paste("A file named", file, "could not be found in the path", dir, "please check your filenames and file paths.", sep = " "))
  }
  if (!str_detect(string = file,
                  pattern = "(_cont_|_pract_|_task_|_full_)")){
    stop(paste("Based upon the name,", file, "does not appear to be data from the practice, control, or experimental task. Please double 
               check your filenames and file paths and try again.", sep = " "))
  }
  
  ## Setting Directory ----
  setwd(dir)
  
  ## Loading data ----
  df <- read.csv(file = file, 
                 na.strings = c("", "NA", " "),  
                 header = T,
                 sep=",", 
                 dec = ".", 
                 stringsAsFactors = F)
  
  ## Creating a Space Saving Variables ----
  if (str_detect(string = file,
                 pattern = "_full_")){
    Task <- ((str_detect(df$Video,
                         "FirstHalf\\.mp4$|FirstHalf_comp\\.mp4$") | 
                str_detect(df$Video,
                           "LastHalf\\.mp4$|LastHalf_comp\\.mp4$")) & !is.na(df$CertRate))
    
    Cont <- (str_detect(df$Video,
                        "Control_*") & !is.na(df$CertRate))
  }
  
  ## Creating a string separating function ----
  strsep <- function(source,
                     output,
                     cond,
                     pattern){
    if (output == "status"){
      array <- unlist(strsplit(x = source[cond],
                               split = "', '", 
                               fixed = T))
      
      array[str_detect(string = array, 
                       pattern = ".'.")] <-str_replace_all(string = array[str_detect(string = array, 
                                                                                     pattern = ".'.")],
                                                           pattern = pattern, 
                                                           replacement = "")
      return(array)
    }
    
    if (output == "rating"){
      array <- unlist(strsplit(x = source[cond],
                               split = ", ", 
                               fixed = T))
      array <-str_replace_all(string = array,
                              pattern = pattern, 
                              replacement = "")
      return(array)
    }
  }
  
  ## Separating Strings ----
  if (str_detect(string = file,
                 pattern = "_full_")){
    CertStat_task <- strsep(source = df$CertStat,
                            output = "status",
                            cond = Task,
                            pattern = "[^[:alnum:]]")
    
    CertRate_task <- strsep(source = df$CertRate,
                            output = "rating",
                            cond = Task,
                            pattern = "'|^\\[|\\]$")
    
    CertStat_cont <- strsep(source = df$CertStat,
                            output = "status",
                            cond = Cont,
                            pattern = "[^[:alnum:]]")
    
    CertRate_cont <- strsep(source = df$CertRate,
                            output = "rating",
                            cond = Cont,
                            pattern = "'|^\\[|\\]$")
  }
  
  if (!str_detect(string = file,
                  pattern = "_full_")){
    CertStat <- strsep(source = df$CertStat,
                       output = "status",
                       cond = !is.na(df$Certainty.Status),
                       pattern = "[^[:alnum:]]")
    
    CertRate <- strsep(source = df$CertRate,
                       output = "rating",
                       cond = !is.na(df$Certainty.Status),
                       pattern = "'|^\\[|\\]$")
  }
  
  ## Creating a longform dataframe ----
  if (str_detect(string = file,
                 pattern = "_full_")){
    rows <- 1:(length(CertStat_task) + length(CertStat_cont))
  } 
  if (!str_detect(string = file,
                  pattern = "_full_")){
    rows <- 1:length(CertStat)
  }
  cols <- c("PID", "Condition", "Video", "Frame", "Seconds", "CertRate", "CertStat", 
            "Time_Video", "Time_Overall", "FrameRate", "Date")
  df_long <- data.frame(matrix(NA, 
                               nrow = length(rows), 
                               ncol = length(cols), 
                               dimnames = list(rows, cols)))
  
  ## Copying data from original dataframe to long form ----
  ### Variables that just copy the value in the same row as video and duplicate it many times ...
  df_long$PID <- rep(df$Participant.[1], length(rows))
  
  ### Correcting PID Errors
  if (any(df_long$PID == "SR-6977") & str_detect(file, pattern = "SR-6799")){
    df_long$PID[df_long$PID == "SR-6977"] <- "SR-6799"
  }
  
  if (str_detect(string = file, pattern = "_cond._")){
    df_long$Condition <- rep(str_extract(string = str_extract(string = file, pattern = "cond."), pattern =".$"), length(rows))
  }
  if (!str_detect(string = file, pattern = "_cond._")){
    df_long$Condition <- NA
  }
  if (!str_detect(string = file, pattern = "_full_")){
    df_long$Video <- rep(df$Video[!is.na(df$Video)], length(rows))
    df_long$FrameRate <- rep(df$frameRate[!is.na(df$Video)], length(rows))
  }
  if (str_detect(string = file, pattern = "_full_")){
    df_long$Video <- c(rep(df$Video[Task], length(CertRate_task)),
                       rep(df$Video[Cont], length(CertRate_cont)))
    df_long$FrameRate <- c(rep(df$frameRate[Task], length(CertRate_task)),
                           rep(df$frameRate[Cont], length(CertRate_cont)))
  }
  df_long$Date <- rep(df$date[1], length(rows))
  df_long$Time_Overall <- rep(df$Offset[length(rownames(df))], length(rows))
  
  ### Variables that just copy pre-existing vectors ...
  if (!str_detect(string = file, pattern = "_full_")){
    df_long$Frame <- rows
    df_long$CertRate <- CertRate
    df_long$CertStat <- CertStat
  }
  
  if (str_detect(string = file, pattern = "_full_")){
    df_long$Frame <- c(1:length(CertRate_task), 1:length(CertRate_cont))
    df_long$CertRate <- c(CertRate_task, CertRate_cont)
    df_long$CertStat <- c(CertStat_task, CertStat_cont)
  }
  
  ### Variables that demand some calculations ...
  if (str_detect(string = file, pattern = "_full_")){
    df_long$Time_Video <- c(rep(df$Offset[Task] - df$Onset[Task], length(CertRate_task)),
                            rep(df$Offset[Cont] - df$Onset[Cont], length(CertRate_cont)))
    df_long$Seconds <- c((df_long$Frame[1:length(CertRate_task)] / max(df_long$Frame[1:length(CertRate_task)])) * df_long$Time_Video[1:length(CertRate_task)],
                         (df_long$Frame[(length(CertRate_task) + 1):length(rownames(df_long))] / max(df_long$Frame[(length(CertRate_task) + 1):length(rownames(df_long))])) * df_long$Time_Video[(length(CertRate_task) + 1):length(rownames(df_long))])
  }
  if (!str_detect(string = file, pattern = "_full_")){
    df_long$Time_Video <- rep(df$isi.Onset[(1:nrow(df))[!is.na(df$Video)] + 1] - df$isi.Offset[(1:nrow(df))[!is.na(df$Video)]], length(rows))
    df_long$Seconds <- (df_long$Frame / max(df_long$Frame)) * df_long$Time_Video
  }
  
  ### Cleaning space ...
  if (str_detect(string = file, pattern = "_full_")){
    rm(CertRate_cont, CertStat_cont, CertStat_task, CertRate_task, rows, cols, Task, Cont, df)
  }
  if (!str_detect(string = file, pattern = "_full_")){
    rm(CertRate, CertStat, CertStat, CertRate, rows, cols, Task, Cont, df)
  }
  
  # Removing observations captured before the shave cutoff ----
  if (shave_secs > 0){
    df_long <- df_long[-which(df_long$Seconds < shave_secs & 
                                !str_detect(df_long$Video, "Control")),]
  }
  
  if (!is.na(unit_secs)){
    # Averaging rating values per second ----
    if (str_detect(string = file, pattern = "_full_")){
      Task_len <- floor(max(df_long$Time_Video[df_long$Video == unique(df_long$Video)[1]])) - shave_secs
      Task_rows <- df_long$Video == unique(df_long$Video)[1]
      Cont_len <- floor(max(df_long$Time_Video[df_long$Video == unique(df_long$Video)[2]]))
      Cont_rows <- df_long$Video == unique(df_long$Video)[2]
    }
    
    ## Creating an average dataframe ----
    if (str_detect(string = file, pattern = "_full_")){
      rows <- c(1:ceiling(Task_len / unit_secs),
                1:ceiling(Cont_len / unit_secs))
    }
    if (!str_detect(string = file, pattern = "_full_")){
      rows <- 1:(max(df_long$Time_Video) / unit_secs)
    }
    cols <- colnames(df_long)[-which(colnames(df_long) == "Seconds")]
    df_avg <- data.frame(matrix(NA,
                                nrow = length(rows),
                                ncol = length(cols),
                                dimnames = list(rows, cols)))
    
    ## Copying data from original dataframe to long form ----
    ### Variables that just copy the value in the same row as video and duplicate it many times ...
    if (str_detect(string = file, pattern = "_full_")){
      df_avg$PID <- rep(df_long$PID[1], length(rows))
      df_avg$Condition <- rep(df_long$Condition[1], length(rows))
      df_avg$Video <- c(rep(unique(df_long$Video)[1], ceiling(Task_len / unit_secs)),
                        rep(unique(df_long$Video)[2], ceiling(Cont_len / unit_secs)))
      df_avg$Time_Video <- c(rep(Task_len, ceiling(Task_len / unit_secs)),
                             rep(Cont_len, ceiling(Cont_len / unit_secs)))
    }
    
    if (!str_detect(string = file, pattern = "_full_")){
      df_avg$PID <- rep(df_long$PID[1], length(rows))
      df_avg$Condition <- rep(df_long$Condition[1], length(rows))
      df_avg$Video <- rep(df_long$Video[1], length(rows))
      df_avg$Time_Video <- rep(df_long$Time_Video[1], length(rows))
    }
    
    df_avg$Time_Overall <- rep(df_long$Time_Overall[1], length(rows))
    df_avg$Date <- rep(df_long$Date[1], length(rows))
    
    ### Variables that just copy pre-existing vectors ...
    for (h in unique(df_avg$Video)){
      target_rows <- which(df_avg$Video == h)
      for (i in 1:length(target_rows)){
        if (!str_detect(h, "Control")) {
          df_avg$SecondStart[target_rows[i]] <- (i * unit_secs) - unit_secs + 0.0001 + shave_secs
          df_avg$SecondEnd[target_rows[i]] <- i * unit_secs + shave_secs
        }
        if (str_detect(h, "Control")) {
          df_avg$SecondStart[target_rows[i]] <- (i * unit_secs) - unit_secs + 0.0001
          df_avg$SecondEnd[target_rows[i]] <- i * unit_secs
        }
        if (df_avg$SecondEnd[target_rows[i]] > max(df_long$Seconds[df_long$Video == h])){
          df_avg$SecondEnd[target_rows[i]] <- max(df_long$Seconds[df_long$Video == h])
        }
      }
    }
    
    
    ### Variables that demand some calculations ...
    if (!str_detect(string = file, pattern = "_full_")){
      df_avg$FrameRate <- rep(max(df_long$Frame)/df_long$Time_Video[1], length(rows))
      for (i in 1:nrow(df_avg)){
        df_avg$Frame[i] <- paste(min(df_long$Frame[df_avg$SecondStart[i] & 
                                                     df_long$Seconds <= df_avg$SecondEnd[i]]),
                                 max(df_long$Frame[df_avg$SecondStart[i] & 
                                                     df_long$Seconds <= df_avg$SecondEnd[i]]),
                                 sep = " - ")
        df_avg$CertRate[i] <- mean(as.numeric(df_long$CertRate[df_avg$SecondStart[i] & 
                                                                 df_long$Seconds <= df_avg$SecondEnd[i]]))
        df_avg$CertRateVar[i] <- var(as.numeric(df_long$CertRate[df_avg$SecondStart[i] & 
                                                                   df_long$Seconds <= df_avg$SecondEnd[i]]))
        if (str_detect(string = file,
                       pattern = "_cont_")){
          if (df_avg$CertRate[i] < 0){
            df_avg$CertStat[i] <- "Left"
          }
          if (df_avg$CertRate[i] > 0){
            df_avg$CertStat[i] <- "Right"
          }
          if (df_avg$CertRate[i] == 0){
            df_avg$CertStat[i] <- "Neither"
          }
        }
        if (str_detect(string = file,
                       pattern = "_pract_|_task_")){
          if (df_avg$CertRate[i] < 0){
            df_avg$CertStat[i] <- "Guilty"
          }
          if (df_avg$CertRate[i] > 0){
            df_avg$CertStat[i] <- "Innocent"
          }
          if (df_avg$CertRate[i] == 0){
            df_avg$CertStat[i] <- "Neutral"
          }
        }
      }
    }
    
    if (str_detect(string = file, pattern = "_full_")){
      df_avg$FrameRate <- c(rep(max(df_long$Frame[Task_rows]) / ceiling(Task_len / unit_secs), ceiling(Task_len / unit_secs)),
                            rep(max(df_long$Frame[Cont_rows]) / ceiling(Cont_len / unit_secs), ceiling(Cont_len / unit_secs)))
      for (h in unique(df_avg$Video)){
        for (i in which(df_avg$Video == h)){
          df_avg$Frame[i] <- paste(min(df_long$Frame[df_long$Seconds >= df_avg$SecondStart[i] & 
                                                       df_long$Seconds <= df_avg$SecondEnd[i]  & 
                                                       df_long$Video == h]),
                                   max(df_long$Frame[df_long$Seconds >= df_avg$SecondStart[i] & 
                                                       df_long$Seconds <= df_avg$SecondEnd[i]  & 
                                                       df_long$Video == h]),
                                   sep = " - ")
          df_avg$CertRate[i] <- mean(as.numeric(df_long$CertRate[df_long$Seconds >= df_avg$SecondStart[i] & 
                                                                   df_long$Seconds <= df_avg$SecondEnd[i]  & 
                                                                   df_long$Video == h]))
          df_avg$CertRateVar[i] <- var(as.numeric(df_long$CertRate[df_long$Seconds >= df_avg$SecondStart[i] & 
                                                                     df_long$Seconds <= df_avg$SecondEnd[i]  & 
                                                                     df_long$Video == h]))
          if (h == unique(df_avg$Video)[1]){
            if (df_avg$CertRate[i] < 0){
              df_avg$CertStat[i] <- "Guilty"
            }
            if (df_avg$CertRate[i] > 0){
              df_avg$CertStat[i] <- "Innocent"
            }
            if (df_avg$CertRate[i] == 0){
              df_avg$CertStat[i] <- "Neutral"
            }
          }
          if (h == unique(df_avg$Video)[2]){
            if (df_avg$CertRate[i] < 0){
              df_avg$CertStat[i] <- "Left"
            }
            if (df_avg$CertRate[i] > 0){
              df_avg$CertStat[i] <- "Right"
            }
            if (df_avg$CertRate[i] == 0){
              df_avg$CertStat[i] <- "Neither"
            }
          }
        }
      }
    }
  }
  
  if (!is.na(unit_secs)){
    ### Cleaning space ...
    rm(rows, cols, h, i, Task_len, Task_rows, Cont_len, Cont_rows, strsep, target_rows) 
  }
  
  ## Output ----
  if (!is.na(unit_secs)){
    return(df_avg)
  }
  if (is.na(unit_secs)){
    return(df_long)
  }
}
