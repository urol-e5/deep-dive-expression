# Debug script to check data structure

# Load libraries
library(dplyr)
library(tidyr)

# Load data
apul_mirna_mrna <- read.csv("D-Apul/output/09-Apul-mRNA-miRNA-interactions/miranda_PCC_miRNA_mRNA.csv")
peve_mirna_mrna <- read.csv("E-Peve/output/10-Peve-mRNA-miRNA-interactions/Peve-miranda_PCC_miRNA_mRNA.csv")
ptuh_mirna_mrna <- read.csv("F-Ptuh/output/11-Ptuh-mRNA-miRNA-interactions/three_prime_interaction/Ptuh-miranda_PCC_miRNA_mRNA.csv")

# Check column names
cat("Apul columns:", colnames(apul_mirna_mrna), "\n")
cat("Peve columns:", colnames(peve_mirna_mrna), "\n")
cat("Ptuh columns:", colnames(ptuh_mirna_mrna), "\n")

# Add species information
apul_mirna_mrna$species <- "A. pulchra"
peve_mirna_mrna$species <- "P. evermanni"
ptuh_mirna_mrna$species <- "P. tuahiniensis"

# Combine data
all_mirna_mrna_interactions <- rbind(apul_mirna_mrna, peve_mirna_mrna, ptuh_mirna_mrna)

# Check combined data
cat("Combined columns:", colnames(all_mirna_mrna_interactions), "\n")
cat("Number of rows:", nrow(all_mirna_mrna_interactions), "\n")

# Check if PCC.cor column exists
if ("PCC.cor" %in% colnames(all_mirna_mrna_interactions)) {
  cat("PCC.cor column found!\n")
  cat("First few values:", head(all_mirna_mrna_interactions$PCC.cor), "\n")
} else {
  cat("PCC.cor column NOT found!\n")
  cat("Available columns:", colnames(all_mirna_mrna_interactions), "\n")
}

