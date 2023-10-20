## cormat_long

cormat_long <- function(data = .){
  
  df <- data %>%
        cor %>%
        as.data.frame %>%
        rownames_to_column(var = 'var1') %>%
        gather(var2, value, -var1) %>%
        mutate(var_order = paste(var1, var2) %>%
                 strsplit(split = ' ') %>%
                 map_chr( ~ sort(.x) %>% 
                            paste(collapse = ' '))) %>%
        mutate(count = 1) %>%
        group_by(var_order) %>%
        mutate(cumsum = cumsum(count)) %>%
        filter(cumsum != 2) %>%
        ungroup %>%
        select(-var_order, -count, -cumsum) %>%
        subset(.$var1 != .$var2)
  
  return(df)
}