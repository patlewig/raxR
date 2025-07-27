#' Generate mol objects from a dataframe of substances and their smiles
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param substances A dataframe of dtxsid identifiers and smiles.
#' @return A rcdk molecular object
#' @export
generate_mol <- function(substances){
  substances_mol <- parse.smiles(as.character(substances$smiles))
  substances_mol <- substances_mol[!sapply(substances_mol, is.null)]
  substances <- substances %>%
    mutate(is_parsed = smiles %in% names(substances_mol))
  
  substances <- substances %>% filter(is_parsed)
  
  names(substances_mol) <- substances$dtxsid
  
  return(substances_mol)
}

#' Generate mol objects from a dataframe of substances and their smiles
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param substances A dataframe of dtxsid identifiers and smiles.
#' @return A rcdk molecular object
#' @export
generate_mol_obj <- function(smiles, dtxsid){
  substance_mol <- parse.smiles(as.character(smiles))
  names(substance_mol) <- as.character(dtxsid)

  return(substance_mol)
}


#' Generate physchem parameters from a set of molecular objects
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param substances_mol A set of molecular objects.
#' @return A dataframe of physchem parameters for the set of substances tagged by their dtxsid identifier
#' @export 
generate_physchem <- function(substances_mol){
  descNames <- c(
    "org.openscience.cdk.qsar.descriptors.molecular.XLogPDescriptor"  ,
    "org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor"   ,
    "org.openscience.cdk.qsar.descriptors.molecular.HBondAcceptorCountDescriptor",
    "org.openscience.cdk.qsar.descriptors.molecular.HBondDonorCountDescriptor")
  descs <- eval.desc(substances_mol, descNames) 
  
  descs <- descs %>%
    rownames_to_column() 
  
  descs <- descs %>% mutate(dtxsid = rowname) %>% select(-c(rowname))
  
  return(descs)
}

#' Generate physchem parameters from a set of molecular objects
#' @import tidyverse
#' @import readxl
#' @import readr
#' @import rcdk
#' @param source_pc Source analogue pchem.
#' @param source_analogues Source analogue identifiers
#' @param target_pc
#' @return A dataframe of physchem parameters for the set of source analogues for a specific target tagged by their dtxsid identifier and role
#' @export 

create_pchem <- function(source_pc, source_analogues, target_pc) {
  source_pc <- source_pc %>% filter(dtxsid %in% source_analogues$dtxsid) %>% mutate(role = 'analogue')
  target_pc <- target_pc  %>% mutate(role = 'target')
  
  return(bind_rows(target_pc, source_pc))
  
}


