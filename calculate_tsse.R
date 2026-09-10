# workflow/scripts/calculate_tsse.R

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Error: Missing arguments. Usage: Rscript calculate_tsse.R <output_txt> <input_tab1> <input_tab2> ...")
}

# The first argument is our designated output summary file
output_file <- args[1]

# All remaining arguments are the paths to the individual sample .tab files
tab_files <- args[2:length(args)]

# Create a data frame to hold the final results
results <- data.frame(Sample = character(), TSSE_Score = numeric(), stringsAsFactors = FALSE)

# Loop through each tab file and calculate the ENCODE TSSE score
for (file_path in tab_files) {
  
  # 1. Extract the sample name from the file name
  file_name <- basename(file_path)
  sample_name <- gsub("_tss_profile.tab", "", file_name)
  
  # 2. Read deepTools profile tab file (skips metadata headers)
  # deepTools has a descriptive text row and column names; skipping 2 rows gets straight to the values
  data <- read.table(file_path, sep = "\t", skip = 2, header = FALSE)
  
  # Convert the row into a clean numeric vector (skipping first 2 columns which hold metadata labels)
  signal <- as.numeric(data[1, 3:ncol(data)])
  
  # 3. Calculate Background according to ENCODE standards
  # Our computeMatrix setup uses 1bp bins from -2000 to +2000 (4001 total positions)
  # Flanks are the first 100bp (-2000 to -1901) and the last 100bp (+1901 to +2000)
  left_flank  <- signal[1:100]
  right_flank <- signal[(length(signal) - 99):length(signal)]
  bg_mean     <- mean(c(left_flank, right_flank))
  
  # Prevent division by zero if background is completely flat/empty
  if (is.na(bg_mean) || bg_mean == 0) {
    bg_mean <- 1
    cat(paste0("The background mean was 0 or NA for sample: "), sample_name, " and it was automatically set to 1 \n")
  }
  
  # 4. Normalize distribution by background
  norm_signal <- signal / bg_mean
  
  # 5. Extract peak value from the absolute center window
  # Index 2001 is the exact center (0 bp / TSS). We search a narrow 11bp window around it to catch shifts.
  center_idx <- 2001
  tsse_score <- max(norm_signal[(center_idx - 5):(center_idx + 5)])
  
  # Append to results table
  results <- rbind(results, data.frame(Sample = sample_name, TSSE_Score = round(tsse_score, 3)))
}

# Write the final summary table out to a tab-separated text file
write.table(results, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)

