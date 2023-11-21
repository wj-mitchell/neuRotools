lmer_extracter <- function(Old_Model = NA, # The base or null model to compare against
                           New_Model, # The new or experimental model to compare 
                           Old_Model_Name = NA, # What you'd like to name the null model
                           New_Model_Name,
                           ICC = T){ # What you'd like to name the new model
  
  # If the pacman package manager is not currently installed on this system, install it.
  if (require("pacman") == FALSE){
    install.packages("pacman")
  }
  
  # Loading in my packages with my pacman manager
  pacman::p_load(broom.mixed,
                 lme4, 
                 lmerTest,
                 performance, 
                 tidyverse)
  
  # ----- LMER MODEL SUMMARY -----
  
  # Pulling the relevant fixed effect statistics from the Compare Model
  results <- tidy(New_Model) %>%
    subset(.$effect == "fixed",
           select = c("term", "estimate", "std.error", "statistic", "df", "p.value"))
  
  # Noting what model these results came from
  results$Model <- New_Model_Name
  
  # Reordering Variables
  results <- subset(results, select = c("Model", "term", "estimate", "std.error", "statistic", "df", "p.value"))
  
  # ----- LMER MODEL ICC -----
  if (ICC == T){
    
    # Calculating the ICC of the New Model
    results$ICC_adj <- as.numeric(icc(New_Model)[1])
  }
  
  # ----- LMER MODEL COMPARISON -----
  
  if (!is.na(Old_Model)){
    
    # Conducting an ANOVA to compare models
    anova <- anova(Old_Model,New_Model)
    
    # Extracting the first and second rows of the ANOVA separately
    df_Old_Model <- anova[1,] %>% select(c("npar", "AIC", "BIC", "logLik", "deviance")) 
    df_New_Model <- anova[2,] %>% select(c("npar", "AIC", "BIC", "logLik", "deviance")) 
    
    # Fitting the model names to the column headers so they don't get confused
    col <- c("ModelComparison",
             paste("Base", names(df_Old_Model), sep = "_"),
             paste("Compare", names(df_New_Model), sep = "_"),
             "Chisq", "Df","Pr(>Chisq)")
    
    # Constructing the new dataframe
    compare <- c(paste(Old_Model_Name, New_Model_Name, sep = "_"),
                 unlist(df_Old_Model), 
                 unlist(df_New_Model),
                 anova[2,which(names(anova) == "Chisq" |
                                 names(anova) == "Df" |
                                 names(anova) == "Pr(>Chisq)" )]) %>%
      t() %>%
      as.data.frame()
    
    # Adding new column names
    names(compare) <- col
  }
  
  # ----- MERGING SUMMARY AND COMPARISON -----
  
  if (is.na(Old_Model)){
    return(results)
  }
  
  if (!is.na(Old_Model)){
    return(as.data.frame(cbind(results,compare)))
  }
}
