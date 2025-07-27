#' Get Neighbourhood
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param targets A file of targets with dtxsid identifiers and smiles 
#' @param source_fps Source fingerprint set
#' @param target_fps Target fingerprint set
#' @return A list of analogues
#' @export
source_lst <- function(targets,source_fps, target_fps, k){
  target_ids <- targets %>% distinct(dtxsid) %>% pull()
  all_analogues <- list()
  for (idx in target_ids) {
    
    all_analogues[[idx]] <- get_analogues(source_fps, target_fps[[idx]],k)
  }
  return(all_analogues)
}





#' Get Neighbourhood List
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param all_analogues List of dataframes for each target and their associated pairwise similarities and dtxsid identifiers
#' @return A vector of source analogues
#' @export
source_array <- function(all_analogues){
  analogue_lst <- list()
  for (idx in target_ids){
    analogue_lst[[idx]] <- all_analogues[[idx]]$dtxsid
  }
  return(analogue_lst)
}


#' Helper function to calculate similarity weighted activity
#' @import tidyverse
#' @import tibble
#' @import readxl
#' @import readr
#' @import dplyr
#' @import rcdk
#' @param id dtxsid
#' @param df The df with activity and similarities
#' @param outcome endpoint of interest
#' @return A tibble with predicted value
#' @export
wtavg <- function(id, df,outcome_col) {
  df <- df %>% filter(!is.na(sim), !is.na(.data[[outcome_col]]),dtxsid != id)
  
  # Calculate weighted average
  value <- sum(df[[outcome_col]] * df$sim, na.rm = TRUE) / sum(df$sim, na.rm = TRUE)
  
  tibble(dtxsid = id, pred_act = value )
}

add_scale <- function (df) {
  normalize <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  df %>% group_by(property) %>%reframe(id = dtxsid,                              
                                       property = property,                   
                                       role = role,                           
                                       property_value = property_value,
                                       property_scaled = normalize(property_value)) 
  
  
  
}

#' Helper function to calculate similarity weighted activity
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param df The df with physchem properties for target and analogues
#' @export

prep_df <- function(df){
  
  test <- pivot_longer(df, cols = c(XLogP,MW, nHBAcc, nHBDon  ), names_to="property", values_to="property_value")
return(test)

}

#' Helper function to compute the SD threshold
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param df A melted df with pchem properties derived from the summary_analysis
#' @param n The sd threshold where 3SD is the default
#' @return A df with source analogues falling outside of the SD threshold
#' @export
filter_sd <- function(df, n=3) {
  # Separate analogues and target data
  analogues <- df %>% filter(role == 'analogue')
  target <- df %>% filter(role == 'target')
  
  # Calculate the mean for each target's property
  target_stats <- target %>%
    group_by(property) %>%
    summarise(target_mean = mean(property_value, na.rm = TRUE))
  
  # Calculate the standard deviation for analogues per property
  analogue_stats <- analogues %>%
    group_by(property) %>%
    summarise(target_std = sd(property_value, na.rm = TRUE))
  
  # Join the target mean and analogue std together
  analogues <- analogues %>%
    left_join(target_stats, by = "property") %>%
    left_join(analogue_stats, by = "property") %>%
    mutate(
      # Check if the analogue is outside n standard deviations of the target mean
      outside_sd = !(property_value >= (target_mean - n * target_std) & property_value <= (target_mean + n * target_std))
    )
  
  # Return only the analogues that are outside the sd range
  analogues %>% filter(outside_sd == TRUE)
  
}
