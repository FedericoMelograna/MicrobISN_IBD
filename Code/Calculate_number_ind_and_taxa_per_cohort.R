library(readr)
library(dplyr)

# Set the base directory
base_dir <- "---"

result_df <- data.frame(Subfolder = character(), NumColumns = integer(), stringsAsFactors = FALSE)
result_df_ind <- data.frame(Subfolder = character(), NumInd = integer(), stringsAsFactors = FALSE)

# Get all folders that start with "group_"
group_folders <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE) %>%
  grep("group_", ., value = TRUE)


w0_base_dir <- "---"


group_folders2 <- list.dirs(w0_base_dir, full.names = TRUE, recursive = FALSE) %>%
  grep("group_", ., value = TRUE)


global_gf = c(group_folders, group_folders2)

# Iterate over each group folder
for (folder in global_gf) {
  
  # Construct the path to the Data subfolder
  data_folder <- file.path(folder, "Data")
  
  # Check if the Data subfolder exists
  if (dir.exists(data_folder)) {
    
    # List all files starting with "binary_result_group_"
    tsv_files <- list.files(data_folder, pattern = "^binary_result_group_.*\\.tsv$", full.names = TRUE)
    
    # Load each TSV file and get the number of columns
    for (tsv_file in tsv_files) {
      data <- read_delim(tsv_file)
      num_columns <- ncol(data)
      
      # Append the subfolder name and the number of columns to the result data frame
      result_df <- rbind(result_df, data.frame(Subfolder = folder, NumColumns = num_columns))
    }
    
    tsv_files <- list.files(data_folder, pattern = "^Group_selection.*\\.tsv$", full.names = TRUE)
    
    # Load each TSV file and get the number of columns
    for (tsv_file in tsv_files) {
      data <- read_delim(tsv_file)
      num_ind <- nrow(data)
      
      # Append the subfolder name and the number of columns to the result data frame
      result_df_ind <- rbind(result_df_ind, data.frame(Subfolder = folder, NumInd = num_ind))
    }
    
  } else {
    print(paste("Data folder does not exist in:", folder))
  }
}

# Print the final result
print(result_df)


merged_ncol_nind = merge(result_df, result_df_ind, by = "Subfolder")

setwd(base_dir)


xlsx::write.xlsx(merged_ncol_nind, file = "Number_individual_and_taxa.xlsx")
