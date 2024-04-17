
pipe_nsym_result <- read.csv("k13_screening_PIPELINE_nsym_mutations_FINAL_17apr2024.csv")


# List all directories in the current working directory
all_dirs <- list.files(path = "../results_v0.1.8_RESMARKERS_FIX/", full.names = TRUE, recursive = FALSE)

# Filter directories containing "*FILTERED"
folder_path <- grep("FILTERED", all_dirs, value = TRUE, fixed = TRUE)

#runs to include in the validation (those with nsym mutations found by the pipeline)
runs_to_include <- unique(pipe_nsym_result$Run)

target_dir <- c("../results_v0.1.8_RESMARKERS_FIX//")

folder_path <- folder_path[folder_path %in% paste0(target_dir, runs_to_include)]


# Function to read CSV files from a folder
read_csv_in_folder <- function(folder_path) {
  csv_files <- list.files(path = folder_path, pattern = "resmarker_table.*filtered.csv", full.names = TRUE)
  df_list <- map(csv_files, read_csv)
  folder_name <- basename(folder_path)
  list(df_name = folder_name, df_list = df_list)
}

# Apply the function to each folder
result_list <- map(folder_path, read_csv_in_folder)

# Add run name
for (i in seq_along(result_list)) {
  print(paste("Processing element", i))
  if (!is.null(result_list[[i]]$df_list) && length(result_list[[i]]$df_list) >= 1) {
    result_list[[i]]$df_list[[1]] <- result_list[[i]]$df_list[[1]] %>%
      mutate(Run = result_list[[i]]$df_name)
  } else {
    print("Error: df_list[[1]] does not exist for this element.")
  }
}

# Extract data frames and assign them to named objects
walk(result_list, ~assign(.x$df_name, bind_rows(.x$df_list), envir = .GlobalEnv))

k13_amps_resmarker_data<- as.data.frame(NULL)

for (run in seq_along(result_list)){
  ok <- result_list[[run]]$df_list
  ok <- as.data.frame(ok)
  #ok <- ok[ok$locus %in% amp_data$amplicon & ok$pseudo_cigar != ".",] #subset k13 amplicons and exclude the ones identical to reference
  
  k13_amps_resmarker_data <- rbind(k13_amps_resmarker_data, ok)
}

k13_amps_resmarker_data <- k13_amps_resmarker_data[k13_amps_resmarker_data$Gene == "k13",]

#VALIDATE. IF MUTATIONS IN PIPE DATA ARE ALSO IN SCREENING DATA, ALL GOOD.
pipe_nsym <- paste0(k13_amps_resmarker_data$SampleID, "__", k13_amps_resmarker_data$CodonID, "__", k13_amps_resmarker_data$RefAA, "__", k13_amps_resmarker_data$AA, "__", k13_amps_resmarker_data$Reads, "__", k13_amps_resmarker_data$Run)

pipe_screening_nsym <- paste0(pipe_nsym_result$sampleID, "__", pipe_nsym_result$non_synonymous_codon, "__", pipe_nsym_result$REF, "__", pipe_nsym_result$ALT, "__", pipe_nsym_result$reads, "__", pipe_nsym_result$Run)

# nsym mutations from the screeining data also found in he pipe run:
pipe_screening_nsym %in% pipe_nsym #ALL TRUE MEANS YES!

#PROOF!
if (length(pipe_screening_nsym) == sum(pipe_screening_nsym %in% pipe_nsym)){
  
  print(paste0("All ", length(pipe_screening_nsym), " nsym muatations from codons already included in the pipeline were found in both the pipeline run and the screening"))
  
}else{
  
  "contact your coke dealer."
  
}
