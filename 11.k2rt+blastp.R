# Load required libraries
library(dplyr)

# Read the TSV files
blastp_data <- read.csv("blastp_filtered.tsv", sep = "\t", stringsAsFactors = FALSE)
huhu_data <- read.csv("Huhu_k2rt_filtered.tsv", sep = "\t", stringsAsFactors = FALSE)

# Step 1: Clean the 'qseqid' column by removing ".p*" patterns
blastp_data$qseqid <- gsub("\\.p[0-9]+", "", blastp_data$qseqid)

# Step 2: Merge the datasets on qseqid and Component_ID
merged_data <- huhu_data %>%
  inner_join(blastp_data %>% select(qseqid, sseqid, ssciname, scomname, sskingdoms), 
             by = c("Component_ID" = "qseqid"))

# Step 3: Export the merged data to a CSV file
write.csv(merged_data, "k2rt_blastp_data.csv", row.names = FALSE)

# Informative message
cat("The merged data has been saved to:", "k2rt_blastp_data.csv", "\n")

data <- read.csv("k2rt_blastp_data.csv")

# Clean the 'sseqid' column by removing 'sp|' from the start and '|' from the end
data <- data %>%
  mutate(sseqid_clean = gsub("^sp\\|", "", sseqid),
         sseqid_clean = gsub("\\|$", "", sseqid_clean))

# Save the cleaned data to a new CSV file
write.csv(data, "cleaned_k2rt_blastp_data.csv", row.names = FALSE)


# Load the dataset
file_path <- "cleaned_k2rt_blastp_data.csv"
data <- read.csv(file_path)

# Remove the . followed by any number from the 'sseqid' column
data$sseqid <- gsub("\\.\\d+", "", data$sseqid)
data <- data %>%
  mutate(sseqid_clean = gsub("^sp\\|", "", sseqid),
         sseqid_clean = gsub("\\|$", "", sseqid_clean))

# Save the cleaned dataset
write.csv(data, "k2rt_blastp.csv", row.names = FALSE)

# Verify the first few rows of the cleaned 'sseqid' column
head(data$sseqid)

# Load the datasets
blastp_file <- "k2rt_blastp.csv"
uniprot_file <- "UniprotMapping.tsv"

# Read the BLASTP and UniProt data
blastp_data <- read.csv(blastp_file)
uniprot_data <- read.delim(uniprot_file, sep = "\t")
colnames(uniprot_data)

# Merge the two datasets based on 'sseqid_clean' from blastp_data and 'Entry' from uniprot_data
merged_data <- blastp_data %>%
  left_join(uniprot_data %>%
              select(`Entry`, `Entry.Name`, `Protein.names`, `Gene.Names`, `Length`),
            by = c("sseqid_clean" = "Entry"))

# Save the merged dataset to a new CSV file
write.csv(merged_data, "merged_k2rt_blastp_uniprot.csv", row.names = FALSE)

# Print the first few rows of the merged dataset to verify
head(merged_data)

# Set the path to the input file (adjust this as needed)
input_file <- "merged_k2rt_blastp_uniprot.csv"

# Read the CSV file
data <- read.csv(input_file)

# Ensure the Is_not_synonymous column is treated as a logical value
# Convert to lower case in case of "True"/"False" strings
data$Is_not_synonymous <- tolower(as.character(data$Is_not_synonymous))

# Split the data into synonymous and non-synonymous mutations
synonymous <- data %>% filter(Is_not_synonymous == "false")
non_synonymous <- data %>% filter(Is_not_synonymous == "true")

# Function to remove the species suffix (_****) from the Entry.Name column
extract_gene <- function(Entry.Name) {
  # Split the entry name by the underscore and take the first part (gene name)
  return(sub("_.*$", "", Entry.Name))
}

# Load the CSV files
non_synonymous_df <- read.csv("non_synonymous_mutations.csv")
synonymous_df <- read.csv("synonymous_mutations.csv")

# Create the "Gene" column in the non-synonymous dataset
non_synonymous_df <- non_synonymous_df %>%
  mutate(Gene = sapply(Entry.Name, extract_gene))

# Create the "Gene" column in the synonymous dataset
synonymous_df <- synonymous_df %>%
  mutate(Gene = sapply(Entry.Name, extract_gene))

# Write the data to separate CSV files
write.csv(synonymous_df, "synonymous_mutations.csv", row.names = FALSE)
write.csv(non_synonymous_df, "non_synonymous_mutations.csv", row.names = FALSE)

# Output message
cat("Synonymous and non-synonymous mutation files have been created.\n")