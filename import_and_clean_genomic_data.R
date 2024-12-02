library(dplyr)


# provice folder name containing all *.FILTERED genomic data
main_dir <- "ALL_SEQ_RESULTS_FILTERED"


# import genomic data
filtered_dirs <- list.dirs(main_dir, recursive = FALSE, full.names = TRUE)

list_of_dfs <- list()

for (dir in filtered_dirs) {
  
  print(paste0("Importing: ", dir))
  
  folder_name <- gsub("_RESULTS.*", "", basename(dir))
  
  csv_file <- file.path(dir, "allele_data_global_max_0_filtered.csv")
  
  if (file.exists(csv_file)) {
    
    df <- read.csv(csv_file)
    df$run <- folder_name
    list_of_dfs[[folder_name]] <- df
  }
}


#merge genomic data into one df
merged_dfs <- bind_rows(list_of_dfs)


#create allele column
merged_dfs$allele <- paste0(merged_dfs$locus, "__", merged_dfs$pseudo_cigar)


# identify sampleID with >=50 loci having >= 100 reads (NANNA'S AND SIMONE'S FILTER)
good_sampleID <- merged_dfs %>%
  
  group_by(sampleID, locus) %>%
  summarize(reads = sum(reads)) %>%
  
  # Filter rows where reads >= 100
  filter(reads >= 100) %>%   

  # Filter NIDA where n_loci >= 50
  summarise(n_loci = n_distinct(locus)) %>%
  filter(n_loci >= 50) %>%

  pull(sampleID)


# Keep good quality samples
merged_dfs <- merged_dfs[merged_dfs$sampleID %in% good_sampleID, ]

#remove shit
MAF = 0.01

merged_dfs <- merged_dfs[!grepl("I=", merged_dfs$allele),] #remove alleles with I (insertion)
merged_dfs <- merged_dfs[!grepl("D=", merged_dfs$allele),] #remove alleles with D (deletion)
merged_dfs <- merged_dfs[!merged_dfs$reads < 10,] #remove alleles with low read counts
merged_dfs <- merged_dfs[!merged_dfs$norm.reads.locus < MAF,] #remove alleles with low read counts

#output
write.csv(merged_dfs, "high_quality_genomic_Data.csv", row.names = F)
