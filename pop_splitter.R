

library(dplyr)
library(stringr)


######################################################################
# input data
db <- read.csv("HFS_ANO1.csv", sep = ";") # metadata db
seq <- read.csv("Allele_merged.csv") # sequencing data

# parameters (you need to input the name of the columns of interest here)
pop <- "district.factor" # select population scale: 1) province.factor, 2) district.factor, 3) us.factor
seq_sample_col <- "SampleID" # column with sample ids in seq data
db_sample_col <- "nida_true" # column with sample ids in db data
######################################################################


# remove rows from the db that don't have population metadata 
db <- db %>%
  filter(!is.na(.data[[pop]]), .data[[pop]] != "")


# fix nidas in seq data
process_sample_data <- function(merged_dfs, sampleID_col) {

  col_name <- sampleID_col
  
  merged_dfs <- merged_dfs %>%
    mutate(
      !!col_name := as.character(.data[[col_name]]),
      !!col_name := gsub("N|_S.*$", "", .data[[col_name]]),
      !!col_name := if_else(
        str_detect(.data[[col_name]], "_") | str_detect(.data[[col_name]], "\\."),
        .data[[col_name]], 
        paste0(.data[[col_name]], ".0")
      ),
      !!col_name := gsub("_", ".", .data[[col_name]])
    )
  return(merged_dfs)
}
  
seq <- process_sample_data(seq, seq_sample_col) # correct nidas from seq df
db <- process_sample_data(db, db_sample_col) # correct nidas fron db


# find common nidas and subset both dfs
commonn_samples <- intersect(unique(seq[[seq_sample_col]]), unique(db[[db_sample_col]]))
seq <- seq[seq[[seq_sample_col]] %in% commonn_samples,]
db <- db[db[[db_sample_col]] %in% commonn_samples,]


# separate population data
pop_list <- list()

pop_factors <- unique(db[[pop]])

for (pop_factor in pop_factors) {
  
  pop_samples <- db %>% 
    filter(.data[[pop]] == pop_factor) %>%
    pull(.data[[db_sample_col]])
  
  pop_seq_data <- seq %>%
    filter(.data[[seq_sample_col]] %in% pop_samples)

  pop_list[[as.character(pop_factor)]] <- pop_seq_data
}


# population summary results
summary_df <- data.frame(
  population = names(pop_list),
  n_unique_samples = sapply(pop_list, function(df) length(unique(df[[seq_sample_col]])))
)

summary_df
write.csv(summary_df, paste0("SUMMARY_", pop), row.names = F)


# output population dfs
for (pop_name in names(pop_list)) {
  filename <- paste0(pop_name, "_seq_data.csv")
  write.csv(pop_list[[pop_name]], file = filename, row.names = FALSE)
  cat("Exported:", filename, "\n")
}

