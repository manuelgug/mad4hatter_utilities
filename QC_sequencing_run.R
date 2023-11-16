# Check if the required packages are installed, if not, install them
if (!requireNamespace("rvest", quietly = TRUE)) {
  install.packages("rvest")
}

if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}

if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr")
}

if (!requireNamespace("optparse", quietly = TRUE)) {
  install.packages("optparse")
}

library(rvest)
library(readxl)
library(tidyr)
library(dplyr)
library(stringr)
library(optparse)

## COMMAND IN TERMINAL: Rscript your_script.R --lab_excel path/to/your/lab_excel.xlsx --qc_html path/to/your/QCplots.html

# Argument parsing
option_list <- list(
  make_option(c("--lab_excel", "-l"), type = "character", default = "POOP", 
              help = "Path to lab_excel file"),
  make_option(c("--qc_html", "-q"), type = "character", default = "PEEP", 
              help = "Path to qc_html file")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Read lab_excel
lab_excel <- read_excel(opt$lab_excel, sheet = "Balancing")

#extract run name
run_name <- lab_excel$RUN[1]

#format lab_excel
colnames(lab_excel)[4]<-"SampleID"
lab_excel <- lab_excel[,c("SampleID", "Well", "Parasitemia", "Group")]

lab_excel <- lab_excel %>%
  mutate(SampleID = paste0("N", gsub(" ", "_", gsub("\\.", "_", SampleID))))

lab_excel$SampleID <- gsub("NNA", NA, lab_excel$SampleID)

lab_excel<- lab_excel %>%
  separate(Well, into = c("row_Well", "column_Well"), sep = "(?<=\\D)(?=\\d)", remove = FALSE)


# Extract tables from HTML
qc_html <- read_html(opt$qc_html)
tables <- html_table(html_nodes(qc_html, "table"))

#cutadapt table
cutadapt_ss <- tables[[2]]
cutadapt_ss <- cutadapt_ss[,-1]
colnames(cutadapt_ss)<- c("SampleID", "Input", "No_Dimers", "Amplicons")

cutadapt_ss <- cutadapt_ss %>%
  mutate(SampleID = ifelse(str_detect(SampleID, "_"), SampleID, paste0(SampleID, "_0")))

cutadapt_ss$SampleID <- gsub("Nn", "NN", cutadapt_ss$SampleID)

# >100 reads table
hundred_or_more <- tables[[6]]

hundred_or_more <- hundred_or_more %>%
  mutate(SampleID = ifelse(str_detect(SampleID, "_"), SampleID, paste0(SampleID, "_0")))

hundred_or_more$SampleID <- gsub("Nn", "NN", hundred_or_more$SampleID)

#write.table(lab_excel, "lab_excel_OUT.csv", row.names = F, sep =",")
#write.table(cutadapt_ss, "cutadapt_ss_OUT.csv", row.names = F, sep =",")
#write.table(hundred_or_more, "hundred_or_more_OUT.csv", row.names = F, sep =",")
#write.table(run_name, "library_OUT.csv", row.names = F, col.names = F, sep=",")

merged_data_FINAL <- merge(merge(hundred_or_more, cutadapt_ss, by = "SampleID", all.x =TRUE), lab_excel, by = "SampleID", all.x =TRUE)

write.table(merged_data_FINAL, paste0("final_lab_df_",run_name, ".csv") , row.names = F, sep =",")

cat("⣿⣟⡿⣟⣿⣻⢿⣻⣟⡿⣟⣿⣻⢿⣻⣟⡿⣟⣿⣻⢿⣻⣟⡿⣟⣿⣻⢿⣻⢿\n",
    "⣿⣞⡿⣯⣷⢿⣻⣽⡾⣟⣯⣷⢿⣻⣽⡾⣟⣯⣷⢿⣻⣽⡾⣟⣯⣷⢿⣻⣯⢿\n",
    "⣿⣞⣿⣽⢾⣟⣯⡷⣿⣻⣽⢾⣟⣯⢿⣯⣍⡙⠺⣿⢯⣷⣟⡿⣽⣾⣻⢷⣻⣿\n",
    "⣿⣞⡿⣞⣯⡿⣞⣿⣳⡿⣽⣻⣾⣯⢿⣞⣿⣽⣦⡀⠉⠛⣾⢿⣽⣞⣿⣻⣽⣾\n",
    "⣿⣞⣿⣻⣽⣻⢯⣷⡿⠙⠁⠀⠀⢨⣽⣟⣾⣽⢾⡿⣦⡀⠀⠙⣿⣾⣳⢿⣳⣿\n",
    "⣿⣾⣹⣷⢿⣹⣿⠏⠀⠀⠀⢀⣾⣿⢿⣾⣹⣾⣿⣹⣿⣷⠀ ⠀⠈⣷⣿⢿⣹⣾\n",
    "⣿⣞⣷⣟⡿⣝⠁⠀⠀⢀⣄⠀⠙⠯⣿⣞⣯⣷⢯⣷⣟⣾⣷⠀⠀ ⠘⣯⣿⢯⣿\n",
    "⣿⣞⡿⣞⣿⣻⣷⣦⣶⣿⡿⣷⣆⡀⠈⠻⣷⣯⢿⣳⣯⡷⣿⠀⠀ ⠀⣹⣯⣿⢾\n",
    "⣿⣞⡿⣯⣷⣟⣷⣻⣽⣾⣻⣽⣻⣷⣄⠀⠈⠙⢿⣽⣳⣿⣻⠀⠀⠀ ⢼⣟⣾⢿\n",
    "⣿⣞⣿⣽⡾⣽⡾⣯⢷⣯⡷⣟⣷⢯⣿⣷⣤⡀⠈⠙⢷⣯⠏⠀⠀⠀⣾⡿⣽⣻\n",
    "⣿⢾⣽⣞⣿⣝⠛⠃⣄⠉⠙⠻⢽⢿⣞⣷⣻⢿⣦⡀⠀⠁⠀⠀⠀ ⣼⣿⣻⣽⢿\n",
    "⣿⢯⣷⡿⠚⠉⢰⣤⡿⣷⣤⣀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀     ⢺⡿⣷⣻\n",
    "⣿⢯⡇⠁⠀⣸⣿⢯⣿⡽⣯⡿⣿⢷⣶⣤⣤⣤⣤⣤⣶⣾⣦⣀⠀⠀⢽⣿⡽⣿\n",
    "⣿⣻⢷⣶⣾⣟⣯⣿⣞⣿⣳⡿⣯⡿⣽⢯⡿⣽⣻⣽⣻⢾⡽⣿⣳⣶⡿⣯⢿⣽\n",
    "⣿⡽⣟⣾⣳⣯⣟⣾⣽⢾⣯⣟⣷⢿⣻⣯⢿⣻⡽⣷⣻⣟⣿⣳⣿⣳⢿⣻⣯⢿\n",
    "⣿⡽⣿⣽⣳⣯⢿⡾⠽⣟⣾⣽⠾⢟⡿⣞⡿⠯⣿⡽⠷⠻⣾⣻⣾⡽⣟⣷⣻⢿")


