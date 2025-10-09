# run_mRNA_lncRNA_PCC.R
message("Starting PCC correlation analysis in chunks...")

# Load libraries
suppressMessages({
  library(tidyverse)
  library(readr)
})

# === Read in gene counts ===
mRNA_counts <- read.csv("../../output/06.2-Peve-Hisat/Peve-gene_count_matrix.csv")
rownames(mRNA_counts) <- mRNA_counts$gene_id
mRNA_counts <- mRNA_counts %>%
  select(-gene_id) %>%
  rename("sample71" = `RNA.POR.71`,
         "sample73" = `RNA.POR.73`,
         "sample76" = `RNA.POR.76`,
         "sample79" = `RNA.POR.79`,
         "sample82" = `RNA.POR.82`)

# Remove mRNAs with 0 for all samples
mRNA_counts <- mRNA_counts %>%
  mutate(Total = rowSums(.)) %>%
  filter(Total != 0) %>%
  select(-Total)

# === Read in lncRNA data ===
lncRNA_counts<-read_table(file="../../../M-multi-species/output/01.6-lncRNA-pipeline/Peve-lncRNA-counts-filtered.txt") %>%
  rename("lncrna_id"=Geneid, 
         "sample71"=`X.home.shared.8TB_HDD_02.zbengt.github.deep.dive.expression.M.multi.species.data.01.6.Peve.lncRNA.pipeline.RNA.POR.71.sorted.bam`, 
         "sample73"=`X.home.shared.8TB_HDD_02.zbengt.github.deep.dive.expression.M.multi.species.data.01.6.Peve.lncRNA.pipeline.RNA.POR.73.sorted.bam`, 
         "sample76"=`X.home.shared.8TB_HDD_02.zbengt.github.deep.dive.expression.M.multi.species.data.01.6.Peve.lncRNA.pipeline.RNA.POR.76.sorted.bam`, 
         "sample79"=`X.home.shared.8TB_HDD_02.zbengt.github.deep.dive.expression.M.multi.species.data.01.6.Peve.lncRNA.pipeline.RNA.POR.79.sorted.bam`, 
         "sample82"=`X.home.shared.8TB_HDD_02.zbengt.github.deep.dive.expression.M.multi.species.data.01.6.Peve.lncRNA.pipeline.RNA.POR.82.sorted.bam`)

lncRNA_counts_df <- as.data.frame(lncRNA_counts) %>%
  select(-Chr, -Start, -End, -Strand, -Length)
rownames(lncRNA_counts_df) <- lncRNA_counts_df$lncrna_id
lncRNA_counts_df <- lncRNA_counts_df %>% select(-lncrna_id)

# Remove lncRNAs with 0 for all samples
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
output_dir <- "../../output/09.1-Peve-mRNA-lncRNA-correlation-PCC/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Function to calculate PCC
calc_pcc <- function(x, y) {
  result <- cor.test(x, y, method = "pearson")
  return(c(PCC = result$estimate, p_value = result$p.value))
}

# Process each chunk
chunk_files <- c()

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

# Combine and adjust p-values
message("Combining chunk results...")
all_chunks <- lapply(chunk_files, read.csv)
combined_results <- bind_rows(all_chunks)

message("Adjusting p-values (FDR)...")
combined_results <- combined_results %>%
  mutate(adjusted_p_value = p.adjust(p_value, method = "fdr"))

# Save final results
final_file <- file.path(output_dir, "Peve-PCC_mRNA_lncRNA.csv")
write.csv(combined_results, final_file, row.names = FALSE)

message("Finished PCC correlation analysis. Results written to: ", final_file)
