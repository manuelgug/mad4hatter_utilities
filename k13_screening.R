library(purrr)
library(readr)
library(dplyr)

# List all directories in the current working directory
all_dirs <- list.files(path = ".", full.names = TRUE, recursive = FALSE)

# Filter directories containing "*FILTERED"
folder_path <- grep("FILTERED", all_dirs, value = TRUE, fixed = TRUE)

# Function to read CSV files from a folder
read_csv_in_folder <- function(folder_path) {
  csv_files <- list.files(path = folder_path, pattern = "resmarker_table.*filtered.csv", full.names = TRUE)
  df_list <- map(csv_files, read_csv)
  folder_name <- basename(folder_path)
  list(df_name = folder_name, df_list = df_list)
}

# Apply the function to each folder
result_list <- map(folder_path, read_csv_in_folder)

#add run name
for (i in seq_along(result_list)) {
  result_list[[i]]$df_list[[1]] <- result_list[[i]]$df_list[[1]] %>%
    mutate(Run = result_list[[i]]$df_name)
}

# Extract data frames and assign them to named objects
walk(result_list, ~assign(.x$df_name, bind_rows(.x$df_list), envir = .GlobalEnv))


results<- as.data.frame(NULL)

for (run in 1:length(result_list)){
  
  ok <- result_list[[run]]$df_list
  ok <- as.data.frame(ok)
  
  ok <- ok[ok$Gene == "k13" & ok$AARefAlt == "ALT",]
  
  results <- rbind(results, ok)
}

write.csv(results, "k13_screening_all_runs.csv", row.names = F)
