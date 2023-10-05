library(tidyr)
library(optparse)
suppressPackageStartupMessages(library(dplyr))

# Define the command-line arguments
option_list <- list(
  make_option(c("--input_res"), type = "character", help = "Input file path (e.g., resmarker_table.txt)" , default ="filtered_results/resmarker_table_global_max_0_filtered.csv"),
  make_option(c("--output_res"), type = "character", help = "Output file path for the formatted data (e.g., resmarker_table_old_format.csv)", default ="filtered_results/resmarker_table_old_format.csv"),
  make_option(c("--input_haps"), type = "character", help = "Input file path (e.g., resmarker_microhap_table.txt)", default ="filtered_results/resmarker_microhap_table_global_max_0_filtered.csv"),
  make_option(c("--output_haps"), type = "character",help = "Output file path for the formatted data (e.g., resmarker_microhap_table_old_format.csv)", default ="filtered_results/resmarker_microhap_table_old_format.csv")
)

# Parse the command-line arguments
opt_parser <- OptionParser(usage = "Usage: %prog [options]", option_list = option_list)
opt <- parse_args(opt_parser)


### RESMARKERS----

# Read the input data
if (!is.null(opt$input_res)) {
  resmarkers <- read.csv(opt$input_res, header = TRUE)

  resmarkers <- resmarkers[, !colnames(resmarkers) %in% c("norm.reads.locus", "n.alleles", "locus")]
  
  #sum reads when needed
  # Group rows by all columns except the last one
  #resmarkers <- resmarkers %>%
  #  group_by(across(-Reads)) %>%
  #  # Sum the last column within each group and replace the original values
  #  mutate(Reads = sum(as.numeric(Reads))) %>%
  #  # Remove duplicates
  #  distinct() %>%
  #  # Reset grouping
  #  ungroup()
  
  resmarkers<-as.data.frame(resmarkers)
  
  #processing
  resmarkers$resmarker <- paste(resmarkers$Gene, resmarkers$CodonID, sep = "_") #resmarkers names
  resmarkers$resmarker_sampleID <- paste(resmarkers$SampleID, resmarkers$resmarker, sep = "_") #column to check for multiple markers in a single sample
  resmarkers$contents <- paste(resmarkers$AA, " [", resmarkers$Reads, "]", sep = "")
  
  #check for missing resmarkers and add NAs when needed before aggregating
  #1) make all possibloe combinations of unique(resmarkers$SampleID) and unique(resmarkers$resmarker)
  combinations <- expand.grid(unique(resmarkers$SampleID), unique(resmarkers$resmarker))
  combined_combinations <- apply(combinations, 1, paste, collapse = "_")
  
  #2) compare to resmarkers$resmarker_sampleID. when the combination is not found, add a new line with the respective sampleID and resmarker and fill ebverything else with NA
  missing_indices <- which(!combined_combinations %in% resmarkers$resmarker_sampleID)
  
  for (index in missing_indices) {
    # Extract the elements from combined_combinations
    elements <- combinations[index,]
    
    # Create a new row with NA values
    new_row <- data.frame(
      SampleID = elements$Var1,
      GeneID = NA,
      Gene = NA,
      CodonID = NA,
      RefCodon = NA,
      Codon = NA,
      CodonStart = NA,
      CodonRefAlt = NA,
      RefAA = NA,
      AA = NA,
      AARefAlt = NA,
      Reads = NA,
      resmarker = elements$Var2,
      resmarker_sampleID = paste(elements$Var1, elements$Var2, sep="_"),
      contents = "[NA]"
    )
    
    # Add the new row to resmarkers
    resmarkers <- rbind(resmarkers, new_row)
  }
  
  #3)as.actor
  # Create a factor variable for AARefAlt with the desired order
  resmarkers$AARefAlt <- factor(resmarkers$AARefAlt, levels = c("REF", "ALT"), ordered = TRUE)
  
  # Sort resmarkers by resmarker_sampleID and AARefAlt
  resmarkers <- resmarkers[order(resmarkers$resmarker_sampleID, resmarkers$AARefAlt, decreasing = F), ]
  
  #### for rows that have the same string in resmarker_sampleID, delete the one with the lowest value in Reads (this is for amplicons from the same region)
  # Subset rows with NA in Reads
  na_rows <- resmarkers %>% filter(is.na(Reads))
  
  # Filter rows by group based on non-NA values
  resmarkers <- resmarkers %>%
    group_by(resmarker_sampleID) %>%
    filter(Reads == max(Reads)) %>%
    ungroup()
  
  # Combine the two subsets
  resmarkers <- bind_rows(resmarkers, na_rows)
  
  #aggregate contents by resmarker_sampleID 
  agg_contents <- aggregate(contents ~ resmarker_sampleID, data = resmarkers, FUN = paste, collapse = ", ")
  
  # Split the agg_contents$resmarker_sampleID by the penultimate "_" into SampleID and resmarker columns
  split_cols <- strsplit(agg_contents$resmarker_sampleID, "_")
  agg_contents$SampleID <- sapply(split_cols, function(x) paste(x[1:(length(x) - 2)], collapse = "_"))
  agg_contents$resmarker <- sapply(split_cols, function(x) paste(x[(length(x) - 1):length(x)], collapse = "_"))
  sampleIDs_unique<-unique(agg_contents$SampleID)
  resmarker_unique<-unique(agg_contents$resmarker)
  
  #create df and fill it with the aggregated contets
  new_df <- as.data.frame(matrix(agg_contents$contents, nrow = length(sampleIDs_unique), ncol = length(resmarker_unique), byrow = TRUE))
  
  # Assign row and column names
  rownames(new_df) <- sampleIDs_unique
  colnames(new_df) <- resmarker_unique
  
  # Extract the column names and split into name and number
  col_names <- colnames(new_df)
  col_names_split <- strsplit(col_names, "_")
  col_names_name <- sapply(col_names_split, function(x) x[1])
  col_names_number <- as.integer(sapply(col_names_split, function(x) x[2]))
  
  # Sort the column names by name and then by number
  sorted_col_names <- col_names[order(col_names_name, col_names_number)]
  
  # Reorder the columns of new_df
  new_df <- new_df[, sorted_col_names]
  
  # Save the formatted data to the output file
  write.csv(new_df, opt$output_res, row.names = TRUE)
  cat("Formatted resmarkers data saved to", opt$output_res, "\n")
}


### HAPLOTYPES----

# Read the input data
if (!is.null(opt$input_haps)) {
  resmarkers <- read.table(opt$input_haps, header = TRUE, sep =",")
  
  resmarkers <- resmarkers[, !colnames(resmarkers) %in% c("norm.reads.locus", "n.alleles")]
  
  #sum reads when needed
  # Group rows by all columns except the last one
  resmarkers <- resmarkers %>%
    group_by(across(-Reads)) %>%
    # Sum the last column within each group and replace the original values
    mutate(Reads = sum(as.numeric(Reads))) %>%
    # Remove duplicates
    distinct() %>%
    # Reset grouping
    ungroup()
  
  resmarkers<-as.data.frame(resmarkers)
  
  #processing
  resmarkers$resmarker <- paste(resmarkers$Gene, resmarkers$MicrohapIndex, sep = "_") #micfrohaps names
  resmarkers$resmarker_sampleID <- paste(resmarkers$SampleID, resmarkers$resmarker, sep = "_") #column to check for multiple microhaps in a single sample
  resmarkers$contents <- paste(resmarkers$Microhaplotype, " [", resmarkers$Reads, "]", sep = "")
  
  #check for missing resmarkers and add NAs when needed before aggregating
  #1) make all possibloe combinations of unique(resmarkers$SampleID) and unique(resmarkers$resmarker)
  combinations <- expand.grid(unique(resmarkers$SampleID), unique(resmarkers$resmarker))
  combined_combinations <- apply(combinations, 1, paste, collapse = "_")
  
  #2) compare to resmarkers$resmarker_sampleID. when the combination is not found, add a new line with the respective sampleID and resmarker and fill ebverything else with NA
  missing_indices <- which(!combined_combinations %in% resmarkers$resmarker_sampleID)
  
  for (index in missing_indices) {
    # Extract the elements from combined_combinations
    elements <- combinations[index,]
    
    # Create a new row with NA values
    new_row <- data.frame(
      SampleID = elements$Var1,
      GeneID = NA,
      Gene = NA,
      MicrohapIndex = NA,
      RefMicrohap = NA,
      Microhaplotype = NA,
      MicrohapRefAlt = NA,
      Reads = NA,
      resmarker = elements$Var2,
      resmarker_sampleID = paste(elements$Var1, elements$Var2, sep="_"),
      contents = "[NA]"
    )
    
    # Add the new row to resmarkers
    resmarkers <- rbind(resmarkers, new_row)
  }
  
  #3)as.actor
  # Create a factor variable for MicrohapRefAlt with the desired order
  resmarkers$MicrohapRefAlt <- factor(resmarkers$MicrohapRefAlt, levels = c("REF", "ALT"), ordered = TRUE)
  
  # Sort resmarkers by resmarker_sampleID and MicrohapRefAlt
  resmarkers <- resmarkers[order(resmarkers$resmarker_sampleID, resmarkers$MicrohapRefAlt, decreasing = F), ]
  
  #aggregate contents by resmarker_sampleID 
  agg_contents <- aggregate(contents ~ resmarker_sampleID, data = resmarkers, FUN = paste, collapse = ", ")
  
  # Split the agg_contents$resmarker_sampleID by the penultimate "_" into SampleID and resmarker columns
  split_cols <- strsplit(agg_contents$resmarker_sampleID, "_")
  agg_contents$SampleID <- sapply(split_cols, function(x) paste(x[1:(length(x) - 2)], collapse = "_"))
  agg_contents$resmarker <- sapply(split_cols, function(x) paste(x[(length(x) - 1):length(x)], collapse = "_"))
  sampleIDs_unique<-unique(agg_contents$SampleID)
  resmarker_unique<-unique(agg_contents$resmarker)
  
  #create df and fill it with the aggregated contets
  new_df <- as.data.frame(matrix(agg_contents$contents, nrow = length(sampleIDs_unique), ncol = length(resmarker_unique), byrow = TRUE))
  
  # Assign row and column names
  rownames(new_df) <- sampleIDs_unique
  colnames(new_df) <- resmarker_unique
  
  # Save the formatted data to the output file
  write.csv(new_df, opt$output_haps, row.names = TRUE)
  cat("Formatted haplotype data saved to", opt$output_haps, "\n")
}
