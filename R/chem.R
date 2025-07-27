#' Preprocess csv file of substances and their identifiers
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param substance_path A file path and csv file name for a set of substances.
#' @return A dataframe of dtxsid identifiers and smiles
#' @export


process_substances <- function(substance_path) {
  substances <- read.csv(substance_path)
  
  # Normalize column names for flexible matching
  names(substances) <- tolower(names(substances))
  
  # Identify SMILES column
  smiles_col <- names(substances)[grepl("smiles", names(substances))]
  if (length(smiles_col) == 0) {
    stop("No SMILES column found. Expected a column like 'Structure_SMILES'.")
  }
  
  # Identify ID column
  id_col <- names(substances)[grepl("dtxsid|id", names(substances))]
  if (length(id_col) == 0) {
    stop("No ID column found. Expected a column like 'DTXSID'.")
  }
  
  # Rename for consistency
  substances <- substances %>%
    rename(
      smiles = all_of(smiles_col[1]),
      dtxsid = all_of(id_col[1])
    ) %>%
    select(dtxsid, smiles)
  substances <- substances %>% distinct(dtxsid, .keep_all = TRUE)
  return(substances)
}

#' Generate Morgan fingerprints for a set of substances characterised by their rcdk molecular objects
#' @import rcdk
#' @param smiles smiles
#' @return A rcdk fingerprint object
#' @export
generate_fp <- function(smiles, dtxsid){
  substance_mol <- parse.smiles(as.character(smiles))
  names(substance_mol) <- as.character(dtxsid)
  single_fps <- lapply(substance_mol, get.fingerprint, type = "circular")
  
  tibble(
    dtxsid = names(single_fps),
    smiles = smiles,
    fingerprint = single_fps
  )
  
}




#' Generate Morgan fingerprints for a set of substances characterised by their rcdk molecular objects
#' @import rcdk
#' @param substances A dataframe of dtxsid and smiles
#' @return A rcdk fingerprint object
#' @export
generate_fingerprints <- function(substances){
  substances_mol <- parse.smiles(as.character(substances$smiles))
  substances_mol <- substances_mol[!sapply(substances_mol, is.null)]
  substances <- substances %>%
    mutate(is_parsed = smiles %in% names(substances_mol))
  
  substances <- substances %>% filter(is_parsed)
  
  names(substances_mol) <- substances$dtxsid
  
  substances_fps <- lapply(substances_mol, get.fingerprint, type = "circular")
  
  return(substances_fps)
  
}

#' Identify candidate source analogues
#' @import rcdk
#' @param sfp Source fingerprint object
#' @param tfp Target fingerprint object
#' @param k Number of source analogues to return
#' @return A dataframe containing the source analogues and their pairwise Jaccard similarities
#' @export
get_analogues <- function(sfp, tfp, k){
  neighbourhood <- data.frame(sim=do.call(rbind, lapply(sfp,
                                               fingerprint::distance,
                                               fp2=tfp, method='tanimoto')))
  neighbourhood <- neighbourhood %>% arrange(desc(sim)) %>% head(k)
  
  neighbourhood <- neighbourhood %>%
    rownames_to_column() %>% mutate(dtxsid = rowname) %>% select(-c(rowname))
  return(neighbourhood)
}
