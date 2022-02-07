#all functions use the tidyverse
require(tidyverse)

split_into_2_seq <- function(data, col){
  #'
  #' Split a pLink Peptide column 
  #' 
  #' @description splits a Petide column into 4 columns, 2 for each peptide
  #' 
  #' @param data a data frame (pLink output table)
  #' @param col the column name of the sequences - 'Peptide'
  #' 
  #' @return The input data frame (data) with 4 new columns: PepPos_A, PepPos_B, Sequence_A, Sequence_B
  #' @note PepPos_# is the crosslinked position within the peptide 
  #' @note Sequence_# is the sequence with square brackets around the crosslinked position (matching XlinkX notation)
  
  
  num_2_bracks <- function(vec){
    pos <- parse_number(vec)
    pep <- str_extract(vec,"\\w+")
    paste0(str_sub(pep,1,pos-1),"[",str_sub(pep,pos,pos),"]",str_sub(pep,pos+1,-1)) #this function does the  ABC(3) to AB[C] conversion
  }
  data %>% separate(col = col, c("Sequence_A", "Sequence_B"), sep = "-") %>%
    mutate(PepPos_A = parse_number(Sequence_A), #takes number in the string in between brackets
           PepPos_B = parse_number(Sequence_B), #num peptide 2
           Sequence_A = num_2_bracks(Sequence_A), #puts square brackets around number position. ABC(3) to AB[C]
           Sequence_B = num_2_bracks(Sequence_B)
           )
}

split_into_2_acc <- function(data, col){
  #'
  #' Split a pLink Proteins column 
  #' 
  #' @description Splits a pLink Proteins column into 6 columns, 3 for each protein.
  #' @description Multiple protein pairs are sorted so intra-links appear first.
  #' @description Discards crosslinks involving contaminants 
  #' 
  #' @param data a data frame (pLink output table)
  #' @param col the column name of the proteins - 'Proteins'
  #' 
  #' @return The input data frame (data) with 6 new columns: 
  #' Position_A, Position_B, Name_A, Name_B, Accession_A, Accession_B
  #' 
  #' @note Discards entries with a contaminant - does not start with a "sp|" 
  #' @note When there is more than one entry in the Proteins column (separated by a "/"), 
  #' values are separated by ";" where intra-links appear first.
  #' @note Position_# is the crosslinked position within the protein
  #' @note Accession_# and Name_# are the database header: sp|accession|name
  
  
  
  intra_to_top <- function(df){
  # helper function to bring an intra option as first
  if(nrow(df) >1){
  df %>% mutate(type = if_else(Name_A == Name_B,0,1)) %>%
    arrange(type) %>%
    select(-type)
  }
    else{df}
  }

  data %>% bind_cols(data[[col]] %>% # add the new columns to the data
                                     str_split("/") %>% # split to the different options
                                     map(function(x){x[-length(x)]}) %>% # discard the empty last entry (all end with a "/")
                                     map(tibble) %>%
                                     map(rename_with,~ col) %>% # rename any column name to col
                                     map(filter, str_detect(!!as.symbol(col), "sp\\|[^/]+-sp\\|")) %>% # discard the entries with a contaminant. see next row.
                                     # The !!as.symbol(col) is so filter will look for the string stored within col and not for a column named "col"
                                     map(parse_pLink_Proteins, col) %>% # parse into different columns
                                     map(intra_to_top) %>% # if there is an intra option bring it to the top
                                     map(summarise, across(.fns = paste, collapse = ";")) %>% # collapse all options with ";"
                                     bind_rows() # merge the list into one data frame
                     )
}

parse_pLink_Proteins <- function(data, col){
  acc_expr <- "[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
  data %>% separate(col = col, into = c("Accession_A", "Accession_B"), sep = "\\)-") %>%
    mutate(Position_A = str_extract(Accession_A,"\\(\\d+") %>% parse_number,
           Position_B = str_extract(Accession_B,"\\(\\d+") %>% parse_number,
           Name_A = str_extract(Accession_A,"(?<=sp\\|\\w{6,10}\\|)\\w+"),
           Name_B = str_extract(Accession_B,"(?<=sp\\|\\w{6,10}\\|)\\w+"))%>%
    mutate(Accession_A = str_extract(Accession_A,acc_expr),
           Accession_B = str_extract(Accession_B,acc_expr))
}

split_title <- function(data, col){
  #' Splits the title on a pLink
  #' 
  #' @description Splits the title column into Spectrum_File, Spectrum_ID, Precursor_ID, scans
  #' 
  #' @param data a data frame - pLink results
  #' @param col the column name of the sequences - 'Title'
  #' 
  #' @return The input data frame (data) with 4 new columns: 
  #' Spectrum_File, Spectrum_ID, Precursor_ID, scans
  data %>% separate(col = col,
                    into = c("Spectrum_File", "Spectrum_ID", "Precursor_ID","scans"),
                    sep = ";") %>%
    mutate(Spectrum_File = str_extract(Spectrum_File,'(?<=").+(?=")'),
           Spectrum_ID = parse_number(Spectrum_ID),
           Precursor_ID = parse_number(Precursor_ID),
           scans = parse_number(scans))
}


CSMs_2_crossID <- function(CSM_data){
  #' Prepare a spectra table for crossID
  #' 
  #' Takes a processed pLink spectra table and creates a table with 
  #' columns and column names suitable for CrossID

  CSM_data %>% select(Crosslinker = Linker, Charge, Score, "MzInDa" = `Precursor_Mass_Error(Da)`,
                      "SequenceA" = Sequence_A, "CrosslinkerPositionA" = PepPos_A,
                      "IdentifierA" = Accession_A, "LeadingProteinPositionA" = Position_A,
                      "SequenceB" = Sequence_B, "CrosslinkerPositionB" = PepPos_B,
                      "IdentifierB" = Accession_B, "LeadingProteinPositionB" = Position_B,
                      "AllScans" = scans, "SpectrumFile" = Spectrum_File)
}

Sites_2_crossID <- function(Sites_data){
  #' Prepare a sites table for crossID
  #' 
  #' Takes a processed pLink sites table and creates a table with 
  #' columns and column names suitable for CrossID
  Sites_data %>% select("AccessionA" = Accession_A, "PositionA" = Position_A,
                      "AccessionB" = Accession_B, "PositionB" = Position_B,
                      "NumCsm" = Spectrum_Number, "MaxXlinkScore" = min_score,
                       "Gene_A" = Name_A, "Gene_B" = Name_B)
}

#Renumber the positions based on the offset_table for sites table - multiple proteins in each row
renumber_plink_multival <- function(data, offset_table){
  #' Renumber the positions based on offset_table
  #' 
  #' @description  Renumber the crosslinked position in the sites table - multi-value (separated by ";") 
  #' @description  Was developed to renumber proteins whose signal peptides were removed from the database
  #' 
  #' @param data a data frame - processed pLink sites table (has Accession_A, Position_A, Accession_B, Position_B)
  #' @param offset_table a data frame with two columns: "ID" (matching Accession_#) and "offset" (the offset from the uniprot sequence) 
  #' 
  #' @return The input data frame (data) with corrected Position_# columns
  
  # helper function to split a row into a dataframe so it can be passed to renumber_plink
  expand_row <- function(rowDF, sep = ";"){
    rowDF %>% str_split(sep) %>% # split on the separator ";"
    as.data.frame(col.names = colnames(rowDF)) %>% # convert the list into a dataframe with the original column names
    type_convert(col_types = cols()) # convert the position column types to numeric - using a readr function (and logic)
}
  # helper function to get a list of DFs of Accession_A, Position_A, Accession_B, Position_B
  get_list_of_DFs <- function(DF){
    DF %>% select(Accession_A, Position_A, Accession_B, Position_B)%>%
           mutate(temp = row_number()) %>% # temp variable that is unique for every row
           group_by(temp) %>% nest() %>% # nests every row as its own DF
           ungroup() %>% select(data) %>% # ungroup to be able to remove temp and select only the data column
           unlist(recursive = F) # extract the list of DFs from the column (not do not unlist that list)
  }

  renumbered_cols <- data %>% get_list_of_DFs() %>% # get a list with a single row DF of every row
      map(expand_row, sep = ";") %>% # expands the single row DF into a multi row DF based on the "; separator"
      map(renumber_plink,offset_table) %>% # renumber based on offset table - now we can use the renumber_plink function since each row has one accession
      map(summarise, across(.fns = paste, collapse = ";")) %>% # undo the expand_row - collapse the columns with ";"
      bind_rows() # merge into one DF
      data<-data %>% mutate(Accession_A = renumbered_cols$Accession_A,
                            Position_A = renumbered_cols$Position_A,
                            Accession_B = renumbered_cols$Accession_B,
                            Position_B = renumbered_cols$Position_B)
}

renumber_plink <- function(data, offset_table){
    indA <-match(data$Accession_A, offset_table$ID)#the index of protein A in the offset_table
    if(indA %>% is.na %>% any){ #in case some accetions are not in the offset table
     stop(paste(paste(data$Accession_A[indA %>% is.na] %>% unique,collapse = ","),
                "not found in offset_table"))}
    
  # Add offset 
    data$Position_A <- data$Position_A + offset_table$Offset[indA]
    
    #get the indexes for B
    indB <-match(data$Accession_B, offset_table$ID) #the index of protein A in the offset_table
    if(indB %>% is.na %>% any){ #in case some accessions are not in the offset table
      stop(paste(paste(data$Accession_B[indB %>% is.na] %>% unique,collapse = ","),
                 "not found in offset_table"))}
  
    # Add offset
    data$Position_B <-data$Position_B + offset_table$Offset[indB]
  data
}

read_pLink_sites <- function(file_path){
  #' extracts a sensible table from pLink sites table
  #' 
  #' @description reads the plink sites table and reshape it into a XL table with nested spetra entries
  #' 
  #' @param file_path a path to the text file with the pLink sites table
  #' 
  #' @return data frame with columns: Protein_Order, Protein, Unique_Peptide_Number, Spectrum_Number, spectra_data
  #' 
  #' @note Protein is lists ("/" separated) of the possible crosslinked proteins
  #' @note spectra_data is a nested data frame which holds the data about the CSMs supporting the crosslink
  
  #read as text lines
  data_raw <- read_lines(file_path)
  # get the XL data - starts with a non blank, so not a comma "[^,]"
  XL_data <- read_csv(data_raw[str_detect(data_raw,"^[^,]")])
  # merge the SameSet and SubSet with their parent row
  XL_data$group_id <- XL_data$Protein_Order # create a group for every crosslink
  for (i in 1:nrow(XL_data)) { 
    if (str_detect(XL_data$group_id[i],"^S.+")) { # replace the SameSet and SubSet with the protein order above it
      XL_data$group_id[i] <- XL_data$group_id[i-1]
    }
  }
  XL_data <- XL_data %>% # map one row for each group
    group_by(group_id) %>% nest() %>% 
    mutate(data = map(data, function(x){
      x$Protein[1] = paste0(c(x$Protein,""),collapse = "/") # collapse with / and add a / at the end - matches the spectra table synthax
      x[1,]})# keep only the first row that now has the Protein data of all rows
    ) %>% unnest(cols = data) %>% ungroup() %>% select(-group_id)

  # get the spectra data - starts with a blank so a comma ","
  spectra_data <- read_csv(c(str_sub(data_raw[str_detect(data_raw,"^,")],2,-1)))
  #add the Protein_Order column, based in the Spectrum number column
  spectra_data <- spectra_data %>% mutate(Protein_Order = rep(XL_data$Protein_Order, times = XL_data$Spectrum_Number)) %>%
    nest(spectra_data = -Protein_Order)
  #add the nested spectra data to the XL table
   XL_data <- full_join(XL_data, spectra_data, by = "Protein_Order")
}

## functions to extract data from the nested spectra data - produced by read_pLink_sites()

extract_spec_files_from_nested_spec_data <- function(spectra_data,exp_pattern){
  #' Extract experiment names from split Title columns of spectra_data
  #' 
  #' @description Extract the experiment names based on exp_pattern from the Spectrum_file column in spectra_data
  #' @description spectra_data is generated as a nested data frame by read_pLink_sites(). 
  #' Spectrum_file is generated by split_title() when splitting the title column
  #' 
  #' @param spectra_data the nested data frame by read_pLink_sites()
  #' @param exp_pattern a regular expression to extract experiment name from file names
  #' 
  #' @return The extracted regular expression for all the spectra, separated by ";"
  spectra_data$Spectrum_File %>%
  str_extract(exp_pattern) %>%
  unique() %>%
  paste(collapse = ";")
}

extract_ID_from_nested_spec_data <- function(spectra_data, col){
  #' Extract a value from spectra_data
  #' 
  #' @description Extract a value from a column in spectra_data
  #' @description spectra_data is generated as a nested data frame by read_pLink_sites(). 
  #' 
  #' @param spectra_data the nested data frame by read_pLink_sites()
  #' @param col the column to extract from
  #' 
  #' @return The extracted values for all the spectra, separated by ";"
  spectra_data[,col] %>%
  unlist() %>%
  paste(collapse = ";")
}

extract_score_from_nested_spec_data <- function(spectra_data, col){
  #' Extract the minimal score from spectra_data
  #' 
  #' @description Extract the minimal value from a column in spectra_data
  #' @description spectra_data is generated as a nested data frame by read_pLink_sites(). 
  #' 
  #' @param spectra_data the nested data frame by read_pLink_sites()
  #' @param col the column to extract from
  #' 
  #' @return The minimal value of col
  spectra_data[,col] %>%
  unlist() %>%
  min()
}



