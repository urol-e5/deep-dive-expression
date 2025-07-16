# run_mRNA_lncRNA_PCC.R
message("Starting PCC correlation analysis in chunks...")

# Load libraries
suppressMessages({
  library(tidyverse)
  library(readr)
})

# === Read in gene counts ===
mRNA_counts <- read.csv("../../output/07-Apul-Hisat/Apul-gene_count_matrix.csv")
rownames(mRNA_counts) <- mRNA_counts$gene_id
mRNA_counts <- mRNA_counts %>%
  select(-gene_id) %>%
  rename("sample140" = `RNA.ACR.140`,
         "sample145" = `RNA.ACR.145`,
         "sample150" = `RNA.ACR.150`,
         "sample173" = `RNA.ACR.173`,
         "sample178" = `RNA.ACR.178`)

# Remove any mRNAs with 0 for all samples
mRNA_counts <- mRNA_counts %>%
  mutate(Total = rowSums(.)) %>%
  filter(Total != 0) %>%
  select(-Total)

# === Read in lncRNA data ===
lncRNA_counts <- read_table(file = "https://raw.githubusercontent.com/urol-e5/deep-dive-expression/refs/heads/main/D-Apul/output/32-Apul-lncRNA-matrix/Apul-lncRNA-counts.txt", skip = 1) %>%
  rename("lncrna_id" = Geneid,
         "sample140" = `../data/32-Apul-lncRNA-matrix/RNA-ACR-140.sorted.bam`,
         "sample145" = `../data/32-Apul-lncRNA-matrix/RNA-ACR-145.sorted.bam`,
         "sample150" = `../data/32-Apul-lncRNA-matrix/RNA-ACR-150.sorted.bam`,
         "sample173" = `../data/32-Apul-lncRNA-matrix/RNA-ACR-173.sorted.bam`,
         "sample178" = `../data/32-Apul-lncRNA-matrix/RNA-ACR-178.sorted.bam`)

lncRNA_counts_df <- as.data.frame(lncRNA_counts) %>%
  select(-Chr, -Start, -End, -Strand, -Length)
rownames(lncRNA_counts_df) <- lncRNA_counts_df$lncrna_id
lncRNA_counts_df <- lncRNA_counts_df %>% select(-lncrna_id)

# Remove any lncRNAs with 0 for all samples
lncRNA_counts_df <- lncRNA_counts_df %>%
  mutate(Total = rowSums(.)) %>%
  filter(Total != 0) %>%
  select(-Total)

# === Normalize counts (RPM) ===
normalize_counts <- function(counts) {
  rpm <- t(t(counts) / colSums(counts)) * 1e6
  return(rpm)
}

mRNA_norm <- normalize_counts(mRNA_counts)
lncRNA_norm <- normalize_counts(lncRNA_counts_df)

# === Set up chunking ===
chunk_size <- 100
lncRNA_ids <- rownames(lncRNA_norm)
num_chunks <- ceiling(length(lncRNA_ids) / chunk_size)

# Output directory
output_dir <- "../../output/21.1-Apul-mRNA-lncRNA-correlation-PCC/"

# Function to calculate PCC
calc_pcc <- function(x, y) {
  result <- cor.test(x, y, method = "pearson")
  return(c(PCC = result$estimate, p_value = result$p.value))
}

# Store chunk filenames
chunk_files <- c()

# === Loop through chunks ===
for (i in seq_len(num_chunks)) {
  message("Processing chunk ", i, " of ", num_chunks)
  
  start_idx <- (i - 1) * chunk_size + 1
  end_idx <- min(i * chunk_size, length(lncRNA_ids))
  lncRNA_chunk_ids <- lncRNA_ids[start_idx:end_idx]
  lncRNA_chunk <- lncRNA_norm[lncRNA_chunk_ids, , drop = FALSE]
  
  pairs <- expand.grid(mRNA = rownames(mRNA_norm), lncRNA = lncRNA_chunk_ids)
  
  chunk_results <- pairs %>%
    rowwise() %>%
    mutate(
      pcc_stats = list(calc_pcc(mRNA_norm[mRNA, ], lncRNA_chunk[lncRNA, ]))
    ) %>%
    unnest_wider(pcc_stats) %>%
    filter(p_value < 0.05)
  
  chunk_file <- file.path(output_dir, paste0("chunk_", i, "_PCC.csv"))
  write.csv(chunk_results, chunk_file, row.names = FALSE)
  chunk_files <- c(chunk_files, chunk_file)
}

# === Combine and adjust p-values ===
message("Combining chunk results...")
all_chunks <- lapply(chunk_files, read.csv)
combined_results <- bind_rows(all_chunks)

message("Adjusting p-values (FDR)...")
combined_results <- combined_results %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "fdr"))

# Save final result
final_file <- file.path(output_dir, "Apul-PCC_mRNA_lncRNA.csv")
write.csv(combined_results, final_file, row.names = FALSE)

message("Finished PCC correlation analysis. Results written to: ", final_file)
