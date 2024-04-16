library(purrr)
library(readr)
library(dplyr)
library(Biostrings)
library(tidyr)
library(DECIPHER)

### PLAN ###
# 1) list k13 amplicons
# 2) input ref gene in codon format + codon table
# 3) extract k13 amplicons from filtered allele data
# 4) align sequences to k13 gene reference
# 5) translate aligned amplicons (reveverse compliment)
# 6) identify all non-synonymous mutations on k13
# 7) final k13 allele_data formatting


# 1) list k13 amplicons
amp_data0 <- read.csv("resistance_markers_amplicon_v4.txt", sep ="\t")
amp_data0 <- amp_data0 %>% 
  filter(grepl("k13", V5))
amp_data <- amp_data0[,c(-1:-6,-16)]
amp_data <- distinct(amp_data)


# 2) input ref gene in codon format + codon table
codon_table <- read.table("codontable.txt")
k13_fasta <- readDNAStringSet("Pf3d7_k13_ref_gene.fasta")
k13_prot_rev_comp <- translate(reverseComplement(k13_fasta)) # 3'5' Frame 1
k13_fasta_rev_comp <- reverseComplement(k13_fasta) # 3'5' Frame 1

#just to be sure the ref is good
k13_fasta_PLASMODB_rev_comp <- readDNAStringSet("Pf3d7_k13_ref_gene_PLASMODB.fasta")
refs_k13_alignment<-AlignSeqs(DNAStringSet(c(k13_fasta_rev_comp, k13_fasta_PLASMODB_rev_comp)))

writeXStringSet(refs_k13_alignment, "refs_k13_alignment.fasta", format = "fasta")


# 3) extract k13 amplicons from filtered allele data
# List all directories in the current working directory
all_dirs <- list.files(path = "../results_v0.1.8_RESMARKERS_FIX/", full.names = TRUE, recursive = FALSE)

# Filter directories containing "*FILTERED"
folder_path <- grep("FILTERED", all_dirs, value = TRUE, fixed = TRUE)

# Function to read CSV files from a folder
read_csv_in_folder <- function(folder_path) {
  csv_files <- list.files(path = folder_path, pattern = "allele_data.*filtered.csv", full.names = TRUE)
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

k13_amps_allele_data<- as.data.frame(NULL)

for (run in seq_along(result_list)){
  ok <- result_list[[run]]$df_list
  ok <- as.data.frame(ok)
  ok <- ok[ok$locus %in% amp_data$amplicon & ok$pseudo_cigar != ".",] #subset k13 amplicons and exclude the ones identical to reference
  
  k13_amps_allele_data <- rbind(k13_amps_allele_data, ok)
}


# 4) align sequences to k13 gene reference
#to avoid doing excessive amounts of alignment, perform only on unique alleles (goes from 3000+ to just 86!), fill table later
unique_alleles <- k13_amps_allele_data[c("locus", "asv")]
unique_alleles <- distinct(unique_alleles)
dim(unique_alleles)

# Loop through each sequence in unique_alleles$asv
amp_rev_comps<-c()

for (seq in unique_alleles$asv) {
  ok <- DNAString(seq)
  rev_comp_seq <- as.character(reverseComplement(ok)) # Reverse complement the sequence
  amp_rev_comps <- c(amp_rev_comps, rev_comp_seq)
}

unique_alleles$amp_reverse_compleemntary <- amp_rev_comps

#get dna alignments
aligned_amp_rev_comp<- c()
aligned_amp_rev_comp_aln<- c()

for (amprevcomp in unique_alleles$amp_reverse_compleemntary){
  
  # Create example DNA sequences as a DNAStringSet
  sequences <- DNAStringSet(c(k13_fasta_rev_comp, amprevcomp))
  
  # Perform pairwise sequence alignment
  alignment <- AlignSeqs(sequences)
  
  aligned_amp_rev_comp <- c(aligned_amp_rev_comp, as.character(alignment[2]))
}

unique_alleles$aligned_amp_rev_comp <- aligned_amp_rev_comp

#output full alignment just because
full_alignment_concat <- DNAStringSet(c(k13_fasta_rev_comp, aligned_amp_rev_comp))
names(full_alignment_concat)[2:251]<- unique_alleles$locus

writeXStringSet(full_alignment_concat, "full_aligment_amplicons.fasta", format = "fasta")


# 5) translate aligned amplicons (reveverse compliment)

getCodons <- function(myAln) {
  seqs <- as.character(myAln)
  len <- width(myAln)[1]
  starts <- seq(from=1, to=len, by=3)
  ends <- starts + 2
  myViews <- lapply(myAln, function(x) { 
    Views(x, starts, ends)
  })
  myCodons <- lapply(myViews, function(x) {
    as.character(DNAStringSet(x))
  })
  myCodons
}

## translateCodons - takes a character vector of codons as input, outputs the corresponding amino acids
translateCodons <- function(myCodons, unknownCodonTranslatesTo="-") {
  ## make new genetic code
  gapCodon <- "-"
  names(gapCodon) <- "---"
  my_GENETIC_CODE <- c(GENETIC_CODE, gapCodon)
  
  ## translate the codons
  pep <- my_GENETIC_CODE[myCodons]
  
  ## check for codons that were not possible to translate, e.g. frameshift codons
  if (sum(is.na(pep))>0) {
    cat("\nwarning - there were codons I could not translate. Using this character", unknownCodonTranslatesTo, "\n\n")
    pep[ which(is.na(pep)) ] <- unknownCodonTranslatesTo
  }
  
  ## prep for output
  pep <- paste(pep, collapse="")
  return(pep)
}

## wrap the getCodons and translateCodons functions together into one:
translateGappedAln <- function(myAln, unknownCodonTranslatesTo="-") {
  myCodons <- getCodons(myAln)
  myAAaln <- AAStringSet(unlist(lapply(myCodons, translateCodons, unknownCodonTranslatesTo=unknownCodonTranslatesTo)))
  return(myAAaln)
}

#translate
aas <- translateGappedAln(unique_alleles$aligned_amp_rev_comp, unknownCodonTranslatesTo="X")
unique_alleles$translated_aligned_amp_rev_comp <- as.character(aas)


# 6) identify all non-synonymous mutations on k13
aa_alignment <- c(k13_prot_rev_comp, aas)

aa_matrix <- matrix("", nrow = length(aa_alignment), ncol = width(aa_alignment[1]))
dim(aa_matrix) #check: good

# Fill the matrix character by character
for (i in 1:length(aa_alignment)) {
  current_seq <- as.character(aa_alignment[i])

  for (j in 1:nchar(current_seq)) {
    aa_matrix[i, j] <- substr(current_seq, j, j)
  }
}

aa_df <- as.data.frame(aa_matrix)
colnames(aa_df) <- 1:ncol(aa_df)

# Compare each sequence with the reference
reference_seq_df <- aa_df[1,]
aa_df <- aa_df[-1,]
rownames(aa_df)<- c()

nsym <- data.frame(rowname = character(), non_synonymous_codon = character(), stringsAsFactors = FALSE)

for (i in 1:nrow(aa_df)) {
  current_seq <- as.character(aa_df[i, ])
  differences <- which(current_seq != reference_seq_df & (current_seq != "-" & current_seq != "X"))
  
  if (length(differences) > 0) {
    nsym <- rbind(nsym, data.frame(rowname = rownames(aa_df)[i], non_synonymous_codon = differences))
  } else {
    nsym <- rbind(nsym, data.frame(rowname = rownames(aa_df)[i], non_synonymous_codon = NA))
  }
}

nsym$ALT <- character(nrow(nsym))
nsym$REF<- character(nrow(nsym))

# Populate 'ALT' column based on the coordinates in nsym
for (i in 1:nrow(nsym)) {
  if (is.na(nsym$non_synonymous_codon[i])) {
    nsym$ALT[i] <- NA
  } else {
    row_num <- as.numeric(nsym$rowname[i])
    col_num <- as.numeric(nsym$non_synonymous_codon[i])
    nsym$ALT[i] <- aa_df[row_num, col_num]
  }
}

# Populate 'REF' column based on the coordinates in nsym
for (i in 1:nrow(nsym)) {
  if (!is.na(nsym$non_synonymous_codon[i])) {
    row_num <- 1
    col_num <- as.numeric(nsym$non_synonymous_codon[i])
    nsym$REF[i] <- reference_seq_df[row_num, col_num]
  } else {

  }
}

unique_alleles$rowname<- 1:length(rownames(unique_alleles))
unique_alleles_complete <- merge(nsym, unique_alleles, by = "rowname", all.x = TRUE)
unique_alleles_complete <- unique_alleles_complete[complete.cases(unique_alleles_complete$ALT), ] #removed NA rows (don't have a nsym mutation so not needed)

# 7) final k13 allele_data formatting

# Merge based on "locus" and "asv"
merged_data <- merge(k13_amps_allele_data, unique_alleles_complete, by = c("locus", "asv"))

merged_data<- merged_data[,c(-11,-15:-18)]

# #remove codons already know to be relevant for resistance
# known_resistance_codons <- as.numeric(sub("^[^-]*-[^-]*-(.*)$", "\\1", amp_data0$V5))
# known_resistance_codons <- as.vector(unique(known_resistance_codons))
# 
# FINAL_TABLE <- merged_data[!(merged_data$non_synonymous_codon %in% known_resistance_codons), ]
# 
# FINAL_TABLE_sorted <- FINAL_TABLE[order(FINAL_TABLE$Run, FINAL_TABLE$locus, FINAL_TABLE$sampleID, FINAL_TABLE$pseudo_cigar), ]
# colnames(FINAL_TABLE_sorted)[4]<-"asv"

# Remove specified columns
colnames_to_remove <- c("amp_reverse_compleemntary", "aligned_amp_rev_comp", "translated_aligned_amp_rev_comp", "rowname", "Category", "allele")
merged_data <- merged_data[, !(names(merged_data) %in% colnames_to_remove)]

write.csv(merged_data, "k13_screening_new_nsym_mutations_FINAL_NEW16apr2024.csv", row.names = F)
