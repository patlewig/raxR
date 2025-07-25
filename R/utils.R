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
#' @import readxl
#' @import readr
#' @import rcdk
#' @param id dtxsid
#' @param df The df with activity and similarities
#' @return A tibble with predicted value
#' @export
wtavg <- function(id, df) {
  df <- df %>% filter(!is.na(sim), !is.na(pPOD))
  
  # Calculate weighted average
  value <- sum(df$pPOD * df$sim) / sum(df$sim)
  
  tibble(dtxsid = id, pred_act = value )
}

