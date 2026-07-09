20-supplementary-files
================
Kathleen Durkin
2026-06-15

- [1 Machinery.fasta reference table](#1-machineryfasta-reference-table)
- [2 Epimachinery category labels for Machinery.fasta
  table](#2-epimachinery-category-labels-for-machineryfasta-table)
- [3 Alignment summary tables](#3-alignment-summary-tables)
  - [3.1 RNA-seq](#31-rna-seq)
    - [3.1.1 1.1 Raw reads](#311-11-raw-reads)
    - [3.1.2 1.2 Trimmed reads](#312-12-trimmed-reads)
    - [3.1.3 1.3 Retained after trimming
      (%)](#313-13-retained-after-trimming-)
    - [3.1.4 1.4 Alignment rate, unique %, multi
      %](#314-14-alignment-rate-unique--multi-)
    - [3.1.5 1.5 Mismatch %](#315-15-mismatch-)
    - [3.1.6 1.6 Number of lncRNAs
      identified](#316-16-number-of-lncrnas-identified)
  - [3.2 sRNA-seq alignment
    statistics](#32-srna-seq-alignment-statistics)
    - [3.2.1 2.1 Species metadata
      (genome)](#321-21-species-metadata-genome)
    - [3.2.2 2.1 Raw reads](#322-21-raw-reads)
    - [3.2.3 2.2 Trimmed reads](#323-22-trimmed-reads)
    - [3.2.4 2.3 Retained after trimming
      (%)](#324-23-retained-after-trimming-)
    - [3.2.5 2.4 ShortStack alignment
      stats](#325-24-shortstack-alignment-stats)
    - [3.2.6 2.5 Mapped miRNA reads](#326-25-mapped-mirna-reads)
    - [3.2.7 2.6 Number of miRNAs
      identified](#327-26-number-of-mirnas-identified)
- [4 Mature miRNA names and sequences for all
  species](#4-mature-mirna-names-and-sequences-for-all-species)
- [5 miRanda binding information + PCC coexpression (miRNA-mRNA and
  miRNA-lncRNA)](#5-miranda-binding-information--pcc-coexpression-mirna-mrna-and-mirna-lncrna)
  - [5.1 miRNA-mRNA](#51-mirna-mrna)
  - [5.2 miRNA-lncRNA](#52-mirna-lncrna)

Code for any compiling/formating of supplementary figures/tables for the
DDE manuscript.

# 1 Machinery.fasta reference table

``` r
protein_names <- c(
  # DNA methylation
  "DNMT1"  = "DNA (cytosine-5-)-methyltransferase 1",
  "DNMT3A" = "DNA (cytosine-5-)-methyltransferase 3 alpha",
  "DNMT3B" = "DNA (cytosine-5-)-methyltransferase 3 beta",
  "DNMT3L" = "DNA (cytosine-5-)-methyltransferase 3 like",
  # TET dioxygenases
  "TET1" = "Methylcytosine dioxygenase TET1",
  "TET2" = "Methylcytosine dioxygenase TET2",
  "TET3" = "Methylcytosine dioxygenase TET3",
  # Methyl-CpG binding domain
  "MBD1" = "Methyl-CpG-binding domain protein 1",
  "MBD2" = "Methyl-CpG-binding domain protein 2",
  "MBD3" = "Methyl-CpG-binding domain protein 3",
  "MBD4" = "Methyl-CpG-binding domain protein 4",
  "MBD5" = "Methyl-CpG-binding domain protein 5",
  "MBD6" = "Methyl-CpG-binding domain protein 6",
  # UHRF
  "UHRF1" = "E3 ubiquitin-protein ligase UHRF1",
  "UHRF2" = "E3 ubiquitin-protein ligase UHRF2",
  # ZBTB readers
  "ZBTB33" = "Zinc finger and BTB domain-containing protein 33 (Kaiso)",
  "ZBTB4"  = "Zinc finger and BTB domain-containing protein 4",
  "ZBTB38" = "Zinc finger and BTB domain-containing protein 38",
  "ZFP57"  = "Zinc finger protein 57",
  "KLF4"   = "Krueppel-like factor 4",
  "EGR1"   = "Early growth response protein 1",
  "WT1"    = "Wilms tumor protein",
  "CTCF"   = "CCCTC-binding factor",
  # HDACs
  "HDAC1"  = "Histone deacetylase 1",
  "HDAC2"  = "Histone deacetylase 2",
  "HDAC3"  = "Histone deacetylase 3",
  "HDAC4"  = "Histone deacetylase 4",
  "HDAC5"  = "Histone deacetylase 5",
  "HDAC6"  = "Histone deacetylase 6",
  "HDAC7"  = "Histone deacetylase 7",
  "HDAC8"  = "Histone deacetylase 8",
  "HDAC9"  = "Histone deacetylase 9",
  "HDAC10" = "Histone deacetylase 10",
  "HDAC11" = "Histone deacetylase 11",
  # Sirtuins
  "SIRT1" = "NAD-dependent protein deacetylase sirtuin-1",
  "SIRT2" = "NAD-dependent protein deacetylase sirtuin-2",
  "SIRT3" = "NAD-dependent protein deacetylase sirtuin-3",
  "SIRT4" = "NAD-dependent protein lysine acylase sirtuin-4",
  "SIRT5" = "NAD-dependent protein lysine desuccinylase sirtuin-5",
  "SIRT6" = "NAD-dependent protein deacylase sirtuin-6",
  "SIRT7" = "NAD-dependent protein desuccinylase sirtuin-7",
  # KATs
  "KAT1-HAT1"    = "Histone acetyltransferase type B catalytic subunit (HAT1/KAT1)",
  "KAT2A"        = "Histone acetyltransferase KAT2A (GCN5)",
  "KAT2B"        = "Histone acetyltransferase KAT2B (PCAF)",
  "KAT3A-CREBBP" = "CREB-binding protein (CREBBP/KAT3A)",
  "KAT3B-EP300"  = "Histone acetyltransferase p300 (EP300/KAT3B)",
  "KAT4-TAF1"    = "Transcription initiation factor TFIID subunit 1 (TAF1/KAT4)",
  "KAT5"         = "Histone acetyltransferase KAT5 (Tip60)",
  "KAT6A"        = "Histone acetyltransferase KAT6A (MOZ)",
  "KAT6B"        = "Histone acetyltransferase KAT6B (MORF)",
  "KAT7"         = "Histone acetyltransferase KAT7 (HBO1)",
  "KAT8"         = "Histone acetyltransferase KAT8 (MOF)",
  "KAT9-ELP3"    = "Elongator complex protein 3 (ELP3/KAT9)",
  "KAT12-GTF3C4" = "Transcription factor IIIC subunit 4 (GTF3C4/KAT12)",
  "KAT13A-NCOA1" = "Nuclear receptor coactivator 1 (NCOA1/KAT13A)",
  "KAT13B-NCOA3" = "Nuclear receptor coactivator 3 (NCOA3/KAT13B)",
  "KAT13C-NCOA2" = "Nuclear receptor coactivator 2 (NCOA2/KAT13C)",
  "KAT13D-CLOCK" = "Circadian locomoter output cycles protein kaput (CLOCK/KAT13D)",
  "KAT14"        = "Histone acetyltransferase KAT14 (CSRP2BP)",
  # KMTs
  "KMT1C-EHMT2"   = "Histone-lysine N-methyltransferase EHMT2 (G9a/KMT1C)",
  "KMT1D-EHMT1"   = "Histone-lysine N-methyltransferase EHMT1 (GLP/KMT1D)",
  "KMT1E-SETDB1"  = "Histone-lysine N-methyltransferase SETDB1 (KMT1E)",
  "KMT1F-SETDB2"  = "Histone-lysine N-methyltransferase SETDB2 (KMT1F)",
  "SUV39H1-KMT1A" = "Histone-lysine N-methyltransferase SUV39H1 (KMT1A)",
  "SUV39H2-KMT1B" = "Histone-lysine N-methyltransferase SUV39H2 (KMT1B)",
  "KMT2A"         = "Histone-lysine N-methyltransferase 2A (MLL1)",
  "KMT2B"         = "Histone-lysine N-methyltransferase 2B (MLL4)",
  "KMT2C"         = "Histone-lysine N-methyltransferase 2C (MLL3)",
  "KMT2D"         = "Histone-lysine N-methyltransferase 2D (MLL2)",
  "KMT2E"         = "Histone-lysine N-methyltransferase 2E (MLL5)",
  "KMT2F-SETD1A"  = "Histone-lysine N-methyltransferase SETD1A (KMT2F)",
  "KMT2G-SETD1B"  = "Histone-lysine N-methyltransferase SETD1B (KMT2G)",
  "KMT2H-ASH1L"   = "Histone-lysine N-methyltransferase ASH1L (KMT2H)",
  "KMT3A-SETD2"   = "Histone-lysine N-methyltransferase SETD2 (KMT3A)",
  "KMT3B-NSD1"    = "Histone-lysine N-methyltransferase NSD1 (KMT3B)",
  "KMT3C-SMYD2"   = "SET and MYND domain-containing protein 2 (SMYD2/KMT3C)",
  "KMT4-DOT1L"    = "Histone-lysine N-methyltransferase DOT1L (KMT4)",
  "KMT5A\u2014SETD8" = "Histone-lysine N-methyltransferase SETD8 (KMT5A)",
  "KMT5B"         = "Histone-lysine N-methyltransferase KMT5B (SUV420H1)",
  "KMT5C"         = "Histone-lysine N-methyltransferase KMT5C (SUV420H2)",
  "KMT6-EZH2"     = "Histone-lysine N-methyltransferase EZH2 (KMT6)",
  "KMT7-SETD7/9"  = "Histone-lysine N-methyltransferase SETD7 (KMT7)",
  "KMT8-PRDM2"    = "PR domain zinc finger protein 2 (PRDM2/KMT8)",
  "ASH2L"         = "Set1/Ash2 histone methyltransferase complex subunit ASH2",
  "WHSC1-NSD2"    = "Histone-lysine N-methyltransferase NSD2 (WHSC1/MMSET)",
  "WHSC1L-NSD3"   = "Histone-lysine N-methyltransferase NSD3 (WHSC1L1)",
  "EZH1"   = "Histone-lysine N-methyltransferase EZH1",
  "SETD3"  = "Histone methyltransferase SETD3",
  "SETD4"  = "SET domain-containing protein 4",
  "SETD5"  = "SET domain-containing protein 5",
  "SETD6"  = "N-lysine methyltransferase SETD6",
  "SETMAR" = "Histone-lysine N-methyltransferase SETMAR (Metnase)",
  "SMYD1"  = "SET and MYND domain-containing protein 1",
  "SMYD3"  = "SET and MYND domain-containing protein 3",
  "SMYD4"  = "SET and MYND domain-containing protein 4",
  "SMYD5"  = "SET and MYND domain-containing protein 5",
  # KDMs
  "KDM1A"         = "Lysine-specific histone demethylase 1A (LSD1)",
  "KDM1B"         = "Lysine-specific histone demethylase 1B (LSD2)",
  "KDM2A"         = "Lysine-specific demethylase 2A",
  "KDM2B"         = "Lysine-specific demethylase 2B",
  "KDM3A"         = "Lysine-specific demethylase 3A",
  "KDM3B"         = "Lysine-specific demethylase 3B",
  "KDM3C-JMJD1C"  = "Lysine-specific demethylase JMJD1C (KDM3C)",
  "KDM4A"         = "Lysine-specific demethylase 4A",
  "KDM4B"         = "Lysine-specific demethylase 4B",
  "KDM4C"         = "Lysine-specific demethylase 4C",
  "KDM4D"         = "Lysine-specific demethylase 4D",
  "KDM5A"         = "Lysine-specific demethylase 5A (JARID1A)",
  "KDM5B"         = "Lysine-specific demethylase 5B (JARID1B)",
  "KDM5C"         = "Lysine-specific demethylase 5C (JARID1C)",
  "KDM5D"         = "Lysine-specific demethylase 5D (JARID1D)",
  "KDM6A"         = "Lysine-specific demethylase 6A (UTX)",
  "KDM6B"         = "Lysine-specific demethylase 6B (JMJD3)",
  "KDM7A"         = "Lysine-specific demethylase 7A (JHDM1D/PHF8-related)",
  "KDM8"          = "Lysine-specific demethylase 8 (JMJD5)",
  "UTY"    = "Inactive histone demethylase UTY",
  "PHF2"   = "PHD finger protein 2 (histone demethylase)",
  "PHF8"   = "Histone lysine demethylase PHF8",
  "JMJD6"  = "Bifunctional arginine demethylase and lysyl-hydroxylase JMJD6",
  "JMJD7"  = "Lysine demethylase JMJD7",
  "JMJD8"  = "JMJD8 protein",
  "JARID2" = "Protein Jumonji (JARID2)",
  "HR"     = "Lysine demethylase and nuclear receptor corepressor HR",
  "RIOX1"  = "Ribosomal oxygenase 1 (NO66)",
  "RIOX2"  = "Ribosomal oxygenase 2 (MINA53)",
  "TYW5"   = "tRNA wybutosine-synthesizing protein 5",
  # PRDM
  "PRDM1"       = "PR domain zinc finger protein 1 (Blimp1)",
  "MECOM-PRDM3" = "EVI1/MDS1 and EVI1 complex locus protein (MECOM/PRDM3)",
  "PRDM4"  = "PR domain zinc finger protein 4",
  "PRDM5"  = "PR domain zinc finger protein 5",
  "PRDM6"  = "PR domain zinc finger protein 6",
  "PRDM7"  = "PR domain zinc finger protein 7",
  "PRDM8"  = "PR domain zinc finger protein 8",
  "PRDM9"  = "PR domain zinc finger protein 9",
  "PRDM10" = "PR domain zinc finger protein 10",
  "PRDM11" = "PR domain zinc finger protein 11",
  "PRDM12" = "PR domain zinc finger protein 12",
  "PRDM13" = "PR domain zinc finger protein 13",
  "PRDM14" = "PR domain zinc finger protein 14",
  "PRDM15" = "PR domain zinc finger protein 15",
  "PRDM16" = "PR domain zinc finger protein 16",
  # PRMTs
  "PRMT2" = "Protein arginine N-methyltransferase 2",
  "PRMT3" = "Protein arginine N-methyltransferase 3",
  "PRMT5" = "Protein arginine N-methyltransferase 5",
  "PRMT6" = "Protein arginine N-methyltransferase 6",
  "PRMT7" = "Protein arginine N-methyltransferase 7",
  "PRMT8" = "Protein arginine N-methyltransferase 8",
  "PRMT9" = "Protein arginine N-methyltransferase 9",
  # PARPs
  "PARP2"          = "Poly [ADP-ribose] polymerase 2",
  "PARP3"          = "Poly [ADP-ribose] polymerase 3",
  "PARP4"          = "Poly [ADP-ribose] polymerase 4",
  "PARP5A-TNKS"    = "Tankyrase-1 (TNKS/PARP5A)",
  "PARP5B-TNKS2"   = "Tankyrase-2 (TNKS2/PARP5B)",
  "PARP6"          = "Poly [ADP-ribose] polymerase 6",
  "PARP7-TIPARP"   = "TCDD-inducible poly [ADP-ribose] polymerase (TIPARP/PARP7)",
  "PARP8"          = "Poly [ADP-ribose] polymerase 8",
  "PARP9"          = "Poly [ADP-ribose] polymerase 9 (BAL1)",
  "PARP10"         = "Poly [ADP-ribose] polymerase 10",
  "PARP11"         = "Poly [ADP-ribose] polymerase 11",
  "PARP12"         = "Poly [ADP-ribose] polymerase 12",
  "PARP13-ZC3HAV1" = "Zinc finger CCCH-type antiviral protein 1 (ZAP/PARP13)",
  "PARP14"         = "Poly [ADP-ribose] polymerase 14 (BAL2)",
  "PARP15"         = "Poly [ADP-ribose] polymerase 15 (BAL3)",
  "PARP16"         = "Poly [ADP-ribose] polymerase 16",
  "PARG"           = "Poly(ADP-ribose) glycohydrolase",
  "ARH1-ADPRH"     = "ADP-ribosylhydrolase ARH1",
  "ARH3-ADPRS"     = "ADP-ribosylhydrolase ARH3",
  "OGT"            = "UDP-N-acetylglucosamine--peptide N-acetylglucosaminyltransferase (OGT)",
  # SUMO
  "SUMO2" = "Small ubiquitin-related modifier 2",
  "SUMO3" = "Small ubiquitin-related modifier 3",
  "SENP1" = "Sentrin-specific protease 1",
  "SENP2" = "Sentrin-specific protease 2",
  "SENP3" = "Sentrin-specific protease 3",
  "SENP5" = "Sentrin-specific protease 5",
  "SENP6" = "Sentrin-specific protease 6",
  "SENP7" = "Sentrin-specific protease 7",
  "SENP8" = "Sentrin-specific protease 8 (NEDD8 deSUMOylation)",
  # Ubiquitin E1/E2/E3
  "UBA1"   = "Ubiquitin-activating enzyme E1",
  "UBA2"   = "SUMO-activating enzyme subunit 2",
  "UBA3"   = "NEDD8-activating enzyme E1 catalytic subunit",
  "UBA5"   = "Ubiquitin-like modifier-activating enzyme 5",
  "UBA6"   = "Ubiquitin-activating enzyme E1-like protein",
  "UBA7"   = "Ubiquitin-activating enzyme E1-like protein 2 (UBE1L)",
  "UBE2A"  = "Ubiquitin-conjugating enzyme E2 A",
  "UBE2D1" = "Ubiquitin-conjugating enzyme E2 D1",
  "UBE2D2" = "Ubiquitin-conjugating enzyme E2 D2",
  "UBE2D3" = "Ubiquitin-conjugating enzyme E2 D3",
  "UBE2V1" = "Ubiquitin-conjugating enzyme E2 variant 1",
  "UBE2Z"  = "Ubiquitin-conjugating enzyme E2 Z",
  "UBE3A"  = "Ubiquitin-protein ligase E3A",
  "UBE4A"  = "Ubiquitin conjugation factor E4 A",
  "UBR2"   = "E3 ubiquitin-protein ligase UBR2",
  "BRCA1"  = "Breast cancer type 1 susceptibility protein (BRCA1)",
  "HUWE1"  = "E3 ubiquitin-protein ligase HUWE1",
  "RNF8"   = "E3 ubiquitin-protein ligase RNF8",
  "RNF20"  = "E3 ubiquitin-protein ligase BRE1A (RNF20)",
  "RNF40"  = "E3 ubiquitin-protein ligase BRE1B (RNF40)",
  "RNF168" = "E3 ubiquitin-protein ligase RNF168",
  "TRIM37" = "E3 ubiquitin-protein ligase TRIM37",
  "DTX3L"  = "Deltex E3 ubiquitin ligase 3L (DTX3L)",
  "CUL4A"  = "Cullin-4A",
  "DDB1"   = "DNA damage-binding protein 1",
  "TOM1"   = "Target of Myb protein 1 (TOM1)",
  # USPs
  "USP1"   = "Ubiquitin carboxyl-terminal hydrolase 1",
  "USP2"   = "Ubiquitin carboxyl-terminal hydrolase 2",
  "USP3"   = "Ubiquitin carboxyl-terminal hydrolase 3",
  "USP4"   = "Ubiquitin carboxyl-terminal hydrolase 4",
  "USP5"   = "Ubiquitin carboxyl-terminal hydrolase 5",
  "USP6"   = "Ubiquitin carboxyl-terminal hydrolase 6",
  "USP7"   = "Ubiquitin carboxyl-terminal hydrolase 7 (HAUSP)",
  "USP8"   = "Ubiquitin carboxyl-terminal hydrolase 8",
  "USP9Y"  = "Ubiquitin carboxyl-terminal hydrolase 9, Y-linked",
  "USP10"  = "Ubiquitin carboxyl-terminal hydrolase 10",
  "USP11"  = "Ubiquitin carboxyl-terminal hydrolase 11",
  "USP12"  = "Ubiquitin carboxyl-terminal hydrolase 12",
  "USP13"  = "Ubiquitin carboxyl-terminal hydrolase 13",
  "USP14"  = "Ubiquitin carboxyl-terminal hydrolase 14",
  "USP15"  = "Ubiquitin carboxyl-terminal hydrolase 15",
  "USP16"  = "Ubiquitin carboxyl-terminal hydrolase 16",
  "USP18"  = "Ubiquitin carboxyl-terminal hydrolase 18",
  "USP19"  = "Ubiquitin carboxyl-terminal hydrolase 19",
  "USP20"  = "Ubiquitin carboxyl-terminal hydrolase 20",
  "USP21"  = "Ubiquitin carboxyl-terminal hydrolase 21",
  "USP22"  = "Ubiquitin carboxyl-terminal hydrolase 22",
  "USP24"  = "Ubiquitin carboxyl-terminal hydrolase 24",
  "USP25"  = "Ubiquitin carboxyl-terminal hydrolase 25",
  "USP26"  = "Ubiquitin carboxyl-terminal hydrolase 26",
  "USP27X" = "Ubiquitin carboxyl-terminal hydrolase 27, X-linked",
  "USP28"  = "Ubiquitin carboxyl-terminal hydrolase 28",
  "USP29"  = "Ubiquitin carboxyl-terminal hydrolase 29",
  "USP30"  = "Ubiquitin carboxyl-terminal hydrolase 30",
  "USP31"  = "Ubiquitin carboxyl-terminal hydrolase 31",
  "USP32"  = "Ubiquitin carboxyl-terminal hydrolase 32",
  "USP33"  = "Ubiquitin carboxyl-terminal hydrolase 33",
  "USP34"  = "Ubiquitin carboxyl-terminal hydrolase 34",
  "USP35"  = "Ubiquitin carboxyl-terminal hydrolase 35",
  "USP36"  = "Ubiquitin carboxyl-terminal hydrolase 36",
  "USP37"  = "Ubiquitin carboxyl-terminal hydrolase 37",
  "USP38"  = "Ubiquitin carboxyl-terminal hydrolase 38",
  "USP39"  = "Ubiquitin carboxyl-terminal hydrolase 39",
  "USP40"  = "Ubiquitin carboxyl-terminal hydrolase 40",
  "USP41"  = "Ubiquitin carboxyl-terminal hydrolase 41",
  "USP42"  = "Ubiquitin carboxyl-terminal hydrolase 42",
  "USP43"  = "Ubiquitin carboxyl-terminal hydrolase 43",
  "USP44"  = "Ubiquitin carboxyl-terminal hydrolase 44",
  "USP45"  = "Ubiquitin carboxyl-terminal hydrolase 45",
  "USP46"  = "Ubiquitin carboxyl-terminal hydrolase 46",
  "USP47"  = "Ubiquitin carboxyl-terminal hydrolase 47",
  "USP48"  = "Ubiquitin carboxyl-terminal hydrolase 48",
  "USP49"  = "Ubiquitin carboxyl-terminal hydrolase 49",
  "USP50"  = "Ubiquitin carboxyl-terminal hydrolase 50",
  "USP51"  = "Ubiquitin carboxyl-terminal hydrolase 51",
  "USP53"  = "Ubiquitin carboxyl-terminal hydrolase 53",
  "USP54"  = "Ubiquitin carboxyl-terminal hydrolase 54",
  "PAN2"   = "Poly(A)-specific ribonuclease subunit PAN2 (USP52)",
  # UCH / OTU / other DUBs
  "UCHL1"   = "Ubiquitin carboxyl-terminal hydrolase isozyme L1",
  "UCHL3"   = "Ubiquitin carboxyl-terminal hydrolase isozyme L3",
  "UCHL4"   = "Ubiquitin carboxyl-terminal hydrolase isozyme L4",
  "UCHL5"   = "Ubiquitin carboxyl-terminal hydrolase isozyme L5 (UCH37)",
  "OTUD1"   = "OTU domain-containing protein 1",
  "OTUD3"   = "OTU domain-containing protein 3",
  "OTUD4"   = "OTU domain-containing protein 4",
  "OTUD5"   = "OTU domain-containing protein 5 (DUBA)",
  "OTUD6A"  = "OTU domain-containing protein 6A",
  "OTUD6B"  = "OTU domain-containing protein 6B",
  "OTUD7A"  = "OTU domain-containing protein 7A (CEZANNE2)",
  "OTUD7B"  = "OTU domain-containing protein 7B (CEZANNE)",
  "BAP1"    = "Ubiquitin carboxyl-terminal hydrolase BAP1",
  "BRCC3"   = "Lys-63-specific deubiquitinase BRCC3",
  "MYSM1"   = "Histone H2A deubiquitinase MYSM1",
  "ATXN3"   = "Ataxin-3 (Josephin domain deubiquitinase)",
  "ATXN3L"  = "Ataxin-3-like protein",
  "JOSD1"   = "Josephin domain-containing protein 1",
  "JOSD2"   = "Josephin domain-containing protein 2",
  "PSMD14"  = "26S proteasome non-ATPase regulatory subunit 14 (POH1/Rpn11)",
  "PSMD7"   = "26S proteasome non-ATPase regulatory subunit 7 (Rpn8)",
  "COPS5"   = "COP9 signalosome complex subunit 5 (JAB1/CSN5)",
  "COPS6"   = "COP9 signalosome complex subunit 6 (CSN6)",
  "STAMBPL1"= "Stam-binding protein-like 1 (AMSH-LP)",
  "ZRANB1"  = "Zinc finger RANBP2-type and RING finger-containing protein 1 (TRABID)",
  "VCPIP1"  = "Deubiquitinating protein VCPIP1",
  "EIF3H"   = "Eukaryotic translation initiation factor 3 subunit H",
  "EIF3F"   = "Eukaryotic translation initiation factor 3 subunit F",
  "MPND"    = "MPN domain-containing protein",
  "PRPF8"   = "Pre-mRNA-processing-splicing factor 8 (PRPF8/SNRNP220)",
  "ARID1A"  = "AT-rich interactive domain-containing protein 1A",
  "ARID1B"  = "AT-rich interactive domain-containing protein 1B",
  # ADAR
  "ADAR"   = "Double-stranded RNA-specific adenosine deaminase (ADAR1)",
  "ADARB1" = "Double-stranded RNA-specific editase B1 (ADAR2)",
  "ADARB2" = "Double-stranded RNA-specific editase B2 (ADAR3)",
  # m6A writers / readers / erasers
  "METTL14"   = "N6-adenosine-methyltransferase non-catalytic subunit (METTL14)",
  "METTL16"   = "RNA N6-methyladenosine methyltransferase METTL16",
  "METTL1"    = "tRNA (guanine(46)-N(7))-methyltransferase METTL1",
  "METTL2A"   = "tRNA (cytosine(34)-N(4))-methyltransferase METTL2A",
  "METTL2B"   = "tRNA (cytosine(34)-N(4))-methyltransferase METTL2B",
  "METTL4"    = "N6-adenosine-methyltransferase METTL4",
  "METTL5"    = "rRNA adenine N(6)-methyltransferase (METTL5)",
  "METTL6"    = "Methyltransferase-like protein 6",
  "METTL8"    = "Methyltransferase-like protein 8",
  "WTAP"      = "pre-mRNA-splicing regulator WTAP (m6A writer complex)",
  "VIRMA"     = "Vir-like m6A methyltransferase associated (KIAA1429)",
  "CBLL1"     = "Hakai (CBLL1); E3 ubiquitin ligase, m6A writer complex",
  "ZC3H13"    = "Zinc finger CCCH domain-containing protein 13 (m6A writer complex)",
  "RBM15"     = "RNA-binding protein 15",
  "RBM15B"    = "RNA-binding protein 15B",
  "PCIF1"     = "mRNA cap-specific m6Am methyltransferase PCIF1",
  "FTO"       = "Alpha-ketoglutarate-dependent dioxygenase FTO (m6A eraser)",
  "ALKBH5"    = "mRNA demethylase ALKBH5 (m6A eraser)",
  "ALKBH1"    = "tRNA N(1)-methyl adenine demethylase ALKBH1",
  "ALKBH3"    = "tRNA demethylase ALKBH3",
  "YTHDC1"    = "YTH domain-containing protein 1 (m6A nuclear reader)",
  "YTHDC2"    = "YTH domain-containing protein 2 (m6A reader/helicase)",
  "YTHDF1"    = "YTH domain-containing family protein 1 (m6A cytoplasmic reader)",
  "YTHDF2"    = "YTH domain-containing family protein 2 (m6A decay reader)",
  "YTHDF3"    = "YTH domain-containing family protein 3 (m6A reader)",
  "HNRNPC"    = "Heterogeneous nuclear ribonucleoprotein C",
  "HNRNPA2B1" = "Heterogeneous nuclear ribonucleoproteins A2/B1",
  "RBMX"      = "RNA-binding motif protein, X-linked (HNRNPG)",
  "IGF2BP1"   = "Insulin-like growth factor 2 mRNA-binding protein 1",
  "IGF2BP2"   = "Insulin-like growth factor 2 mRNA-binding protein 2",
  "IGF2BP3"   = "Insulin-like growth factor 2 mRNA-binding protein 3",
  "FMR1"      = "Fragile X mental retardation protein 1",
  "PRRC2A"    = "Proline-rich coiled-coil-containing protein 2A (m6A reader)",
  "ELF3"      = "ETS-related transcription factor ELF3",
  # RNA methylation (non-m6A)
  "NSUN2"   = "tRNA (cytosine(34)-C(5))-methyltransferase NSUN2",
  "NSUN3"   = "Methyltransferase-like protein NSUN3 (mitochondrial tRNA)",
  "NSUN4"   = "Bifunctional methyltransferase/rRNA maturation factor NSUN4",
  "NSUN5"   = "rRNA methyltransferase NSUN5",
  "NSUN6"   = "tRNA methyltransferase NSUN6",
  "TRDMT1"  = "tRNA (cytosine-5-)-methyltransferase TRDMT1 (DNMT2)",
  "NOP2"    = "Ribosome biogenesis protein NOP2 (NSUN1)",
  "DKC1"    = "H/ACA ribonucleoprotein complex subunit 4 (Dyskerin/DKC1)",
  "TRMT1"   = "tRNA (guanine(26)-N(2))-methyltransferase TRMT1",
  "TRMT61A" = "tRNA (adenine(58)-N(1))-methyltransferase TRMT61A",
  "TRMT61B" = "rRNA (adenine(58)-N(1))-methyltransferase TRMT61B",
  "TRMT10A" = "tRNA (guanine(9)-N(1))-methyltransferase TRMT10A",
  "TRMT10C" = "Mitochondrial tRNA methyltransferase TRMT10C",
  "TRMT112" = "tRNA methyltransferase subunit TRMT112",
  "ALYREF"  = "THO complex subunit ALYREF (Aly/REF export factor)",
  "ZCCHC4"  = "Zinc finger CCHC domain-containing protein 4 (rRNA m6A writer)",
  "NAT10"   = "N-acetyltransferase 10 (RNA cytidine acetyltransferase)",
  "DUS1L"   = "tRNA-dihydrouridine synthase 1-like",
  # Pseudouridine synthases
  "PUSL1"  = "Pseudouridine synthase 1-like",
  "PUS1"   = "tRNA pseudouridine synthase A (PUS1)",
  "PUS3"   = "tRNA pseudouridine synthase C (PUS3)",
  "PUS7"   = "tRNA pseudouridine synthase G (PUS7)",
  "PUS7L"  = "tRNA pseudouridine synthase 7-like",
  "PUS10"  = "tRNA pseudouridine synthase 10 (RPUSD2-related)",
  "RPUSD1" = "RNA pseudouridylate synthase domain-containing 1",
  "RPUSD2" = "RNA pseudouridylate synthase domain-containing 2",
  "RPUSD3" = "RNA pseudouridylate synthase domain-containing 3",
  "RPUSD4" = "RNA pseudouridylate synthase domain-containing 4",
  "TRUB1"  = "Probable tRNA pseudouridine synthase 1 (TRUB1)",
  "TRUB2"  = "Probable tRNA pseudouridine synthase 2 (TRUB2)",
  # NATs / acetyltransferases
  "NAA10"   = "N-alpha-acetyltransferase 10 (NatA catalytic subunit)",
  "NAA50"   = "N-alpha-acetyltransferase 50 (NatE catalytic subunit)",
  "NAA60"   = "N-alpha-acetyltransferase 60 (NatF catalytic subunit)",
  "NAT8"    = "N-acetyltransferase 8 (camello-like 1)",
  "NAT8B"   = "N-acetyltransferase 8B (camello-like 2)",
  "NAT8L"   = "N-acetyltransferase 8-like",
  "ATAT1"   = "Alpha-tubulin acetyltransferase 1",
  "GTF3C2"  = "Transcription factor IIIC 110 kDa subunit (TFIIIC110/GTF3C2)",
  "GTF3C1"  = "Transcription factor IIIC 220 kDa subunit (TFIIIC220/GTF3C1)",
  "DLAT"    = "Dihydrolipoyllysine-residue acetyltransferase (PDC-E2)",
  "ATF2"    = "Cyclic AMP-dependent transcription factor ATF-2",
  "MCM3AP"  = "Germinal-center associated nuclear protein (GANP/MCM3AP)",
  "MAPT"    = "Microtubule-associated protein tau (MAPT)",
  "BRD4"    = "Bromodomain-containing protein 4",
  "ESCO1"   = "N-acetyltransferase ESCO1 (establishment of cohesion)",
  "ESCO2"   = "N-acetyltransferase ESCO2",
  "BLOC1S1" = "Biogenesis of lysosome-related organelles complex 1 subunit 1",
  "ACAT1"   = "Acetyl-CoA acetyltransferase 1 (mitochondrial)",
  "ACAT2"   = "Acetyl-CoA acetyltransferase 2 (cytoplasmic)",
  # Histone variants / linker histones
  "H1-0"  = "Histone H1.0",
  "H1-1"  = "Histone H1.1 (HIST1H1A)",
  "H1-2"  = "Histone H1.2 (HIST1H1C)",
  "H1-3"  = "Histone H1.3 (HIST1H1D)",
  "H1-4"  = "Histone H1.4 (HIST1H1E)",
  "H1-5"  = "Histone H1.5 (HIST1H1B)",
  "H1-6"  = "Histone H1t",
  "H1-7"  = "Histone H1oo (oocyte-specific)",
  "H1-8"  = "Histone H1.8 (H1foo)",
  "H1-10" = "Histone H1.10",
  "H2AX"  = "Histone H2A.x",
  "H2AZ1" = "Histone H2A.z.1",
  "H2AZ2" = "Histone H2A.z.2",
  "H2AW"  = "Histone H2A.W",
  "H2AP"  = "Histone H2A.P",
  "H2AJ"  = "Histone H2A.J",
  "H2AB1" = "Histone H2A-Bbd type 1 (H2A.Bbd)",
  "H2AB2" = "Histone H2A-Bbd type 2",
  "H2AB3" = "Histone H2A-Bbd type 3",
  "MACROH2A1" = "Core histone macro-H2A.1 (macroH2A1/H2AFY)",
  # Chromatin regulators
  "BAZ1B" = "Bromodomain adjacent to zinc finger domain protein 1B (WSTF)",
  "BAZ2B" = "Bromodomain adjacent to zinc finger domain protein 2B",
  "HIRA"  = "Histone cell cycle regulator (HIRA)",
  "WDR5"  = "WD repeat-containing protein 5 (H3K4 methyltransferase complex subunit)",
  # EYA phosphatases
  "EYA1" = "Eyes absent homolog 1 (H2AX phosphatase)",
  "EYA2" = "Eyes absent homolog 2",
  "EYA3" = "Eyes absent homolog 3",
  "EYA4" = "Eyes absent homolog 4",
  # Protein phosphatases
  "PPM1D"    = "Protein phosphatase 1D (WIP1)",
  "PPP6C"    = "Serine/threonine-protein phosphatase 6 catalytic subunit",
  "PPP4C"    = "Serine/threonine-protein phosphatase 4 catalytic subunit",
  "CDCA2"    = "Cell division cycle associated 2 (Repo-Man; PP1 targeting)",
  "PPP1CA"   = "Serine/threonine-protein phosphatase PP1-alpha catalytic subunit",
  "PPP1CB"   = "Serine/threonine-protein phosphatase PP1-beta catalytic subunit",
  "PPP1CC"   = "Serine/threonine-protein phosphatase PP1-gamma catalytic subunit",
  "PPP1R1B"  = "Protein phosphatase 1 regulatory subunit 1B (DARPP-32)",
  "PPP1R1C"  = "Protein phosphatase 1 regulatory subunit 1C (I-3)",
  "PPP1R2"   = "Protein phosphatase inhibitor 2",
  "PPP1R3A"  = "Protein phosphatase 1 regulatory subunit 3A (GM)",
  "PPP1R3B"  = "Protein phosphatase 1 regulatory subunit 3B (GL)",
  "PPP1R3C"  = "Protein phosphatase 1 regulatory subunit 3C (PTG)",
  "PPP1R3D"  = "Protein phosphatase 1 regulatory subunit 3D (R6)",
  "PPP1R3E"  = "Protein phosphatase 1 regulatory subunit 3E",
  "PPP1R3F"  = "Protein phosphatase 1 regulatory subunit 3F",
  "PPP1R3G"  = "Protein phosphatase 1 regulatory subunit 3G",
  "PPP1R7"   = "Protein phosphatase 1 regulatory subunit 7",
  "PPP1R8"   = "Protein phosphatase 1 regulatory subunit 8 (NIPP1)",
  "PPP1R9A"  = "Protein phosphatase 1 regulatory subunit 9A (neurabin-1)",
  "PPP1R9B"  = "Protein phosphatase 1 regulatory subunit 9B (neurabin-2/spinophilin)",
  "PPP1R10"  = "Protein phosphatase 1 regulatory subunit 10 (PNUTS)",
  "PPP1R11"  = "Protein phosphatase 1 regulatory subunit 11 (PNUTS-like)",
  "PPP1R12A" = "Protein phosphatase 1 regulatory subunit 12A (MYPT1)",
  "PPP1R12B" = "Protein phosphatase 1 regulatory subunit 12B (MYPT2)",
  "PPP1R12C" = "Protein phosphatase 1 regulatory subunit 12C (MBS85)",
  "PPP1R13B" = "Protein phosphatase 1 regulatory subunit 13B (ASPP1)",
  "PPP1R14A" = "Protein phosphatase 1 regulatory subunit 14A (CPI-17)",
  "PPP1R14B" = "Protein phosphatase 1 regulatory subunit 14B (PHI-1)",
  "PPP1R14C" = "Protein phosphatase 1 regulatory subunit 14C",
  "PPP1R14D" = "Protein phosphatase 1 regulatory subunit 14D",
  "PPP1R15B" = "Protein phosphatase 1 regulatory subunit 15B (GADD34-related)",
  "PPP1R16A" = "Protein phosphatase 1 regulatory subunit 16A (MYPT3)",
  "PPP1R16B" = "Protein phosphatase 1 regulatory subunit 16B (TA-MYPT)",
  # Kinases
  "AURKB"   = "Aurora kinase B",
  "AURKA"   = "Aurora kinase A",
  "AURKC"   = "Aurora kinase C",
  "JAK2"    = "Tyrosine-protein kinase JAK2",
  "HASPIN"  = "Serine/threonine-protein kinase haspin (GSG2)",
  "BUB1"    = "Mitotic checkpoint serine/threonine-protein kinase BUB1",
  "WEE1"    = "Wee1-like protein kinase",
  "SGO2"    = "Shugoshin-like 2 (SGO2)",
  "SGO1"    = "Shugoshin-like 1 (SGO1)",
  "DLK1"    = "Protein delta homolog 1 (Preadipocyte factor 1/Pref-1)",
  "DAPK3"   = "Death-associated protein kinase 3 (ZIPK)",
  "PKN1"    = "Serine/threonine-protein kinase N1 (PKN1/PRK1)",
  "PKM"     = "Pyruvate kinase PKM",
  "STK4"    = "Serine/threonine-protein kinase 4 (MST1/STK4)",
  "MST1"    = "Hepatocyte growth factor-like protein (MST1)",
  "SLK"     = "STE20-like serine/threonine-protein kinase SLK",
  "PRKAA1"  = "5'-AMP-activated protein kinase catalytic subunit alpha-1 (AMPKa1)",
  "PRKAA2"  = "5'-AMP-activated protein kinase catalytic subunit alpha-2 (AMPKa2)",
  "DBF4"    = "Protein DBF4 homolog A (Cdc7 kinase regulatory subunit)",
  "DYRK1A"  = "Dual specificity tyrosine-phosphorylation-regulated kinase 1A",
  "RPS6KA4" = "Ribosomal protein S6 kinase alpha-4 (MSK2)",
  "RPS6KA5" = "Ribosomal protein S6 kinase alpha-5 (MSK1)",
  "PRKCB"   = "Protein kinase C beta type (PKCb)",
  "JKAMP"   = "JMJD8-associated regulatory complex component (JKAMP)",
  "MAPK8"   = "Mitogen-activated protein kinase 8 (JNK1)",
  "MAPK9"   = "Mitogen-activated protein kinase 9 (JNK2)",
  "MDC1"    = "Mediator of DNA damage checkpoint protein 1",
  "MCPH1"   = "Microcephalin (MCPH1; DNA damage sensor)",
  "PAK2"    = "Serine/threonine-protein kinase PAK2",
  "ATM"     = "Serine-protein kinase ATM (DNA damage checkpoint)",
  "ATR"     = "Serine/threonine-protein kinase ATR",
  "CSNK2B"  = "Casein kinase II subunit beta",
  "PRKCD"   = "Protein kinase C delta type",
  "MAP2K1"  = "Dual specificity mitogen-activated protein kinase kinase 1 (MEK1)",
  "MAP3K12" = "Mitogen-activated protein kinase kinase kinase 12 (DLK/MUK)",
  "MAP3K22" = "Mitogen-activated protein kinase kinase kinase 22",
  "MAPK1"   = "Mitogen-activated protein kinase 1 (ERK2)",
  "MAPK3"   = "Mitogen-activated protein kinase 3 (ERK1)",
  "MAPK8B"  = "Mitogen-activated protein kinase 8b (JNK1b, zebrafish)",
  "MAPK11"  = "Mitogen-activated protein kinase 11 (p38b)",
  "MAPK12"  = "Mitogen-activated protein kinase 12 (p38g/ERK6)",
  "MAPK13"  = "Mitogen-activated protein kinase 13 (p38d/SAPK4)",
  "MAPK14"  = "Mitogen-activated protein kinase 14 (p38a)",
  "PRKDC"   = "DNA-dependent protein kinase catalytic subunit (DNA-PKcs)",
  "PIM1"    = "Proto-oncogene serine/threonine-protein kinase Pim-1",
  "CDK8"    = "Cyclin-dependent kinase 8 (Mediator kinase)",
  "CHUK"    = "Inhibitor of nuclear factor kappa-B kinase subunit alpha (IKKa)",
  "CDK1"    = "Cyclin-dependent kinase 1 (CDC2)",
  "CHEK1"   = "Serine/threonine-protein kinase Chk1",
  "APBB1"   = "Amyloid precursor protein-binding family B member 1 (Fe65)",
  "YWHAQ"   = "14-3-3 protein theta",
  # ZC3H12 ribonucleases
  "ZC3H12A" = "Endoribonuclease ZC3H12A (MCPIP1/Regnase-1)",
  "ZC3H12B" = "Endoribonuclease ZC3H12B (Regnase-2)",
  "ZC3H12C" = "Endoribonuclease ZC3H12C (Regnase-3)",
  "ZC3H12D" = "Endoribonuclease ZC3H12D (Regnase-4)",
  # Other
  "DCP2"   = "mRNA-decapping enzyme subunit DCP2",
  "YBX1"   = "Nuclease-sensitive element-binding protein 1 (YBX1/YB-1)",
  "RAG1"   = "V(D)J recombination-activating protein 1 (RAG1)",
  "ALG13"  = "UDP-N-acetylglucosamine transferase subunit ALG13",
  "HSPBAP1"= "HSPA binding protein 1 (HSPBAP1)"
)

strip_isoform <- function(label) {
  # Remove trailing isoform number, e.g. Dnmt1-201 -> Dnmt1
  sub("-\\d+$", "", label)
}

species_from_id <- function(eid) {
  dplyr::case_when(
    startsWith(eid, "ENSMUSP") ~ "Mus musculus",
    startsWith(eid, "ENSDARP") ~ "Danio rerio",
    TRUE                       ~ "Homo sapiens"
  )
}

lookup_protein_name <- function(gene_name) {
  key <- toupper(gene_name)
  if (!is.na(protein_names[key])) return(unname(protein_names[key]))
  # Strip alias suffixes and try progressively shorter hyphen-joined prefixes
  simplified <- gsub("[\u2014\u2013/].*", "", key)
  parts <- strsplit(simplified, "-")[[1]]
  for (i in seq(length(parts), 1)) {
    candidate <- paste(parts[seq_len(i)], collapse = "-")
    if (!is.na(protein_names[candidate])) return(unname(protein_names[candidate]))
  }
  ""
}

# Parse fasta and build table
lines   <- readLines("../../data/Machinery.fasta")
```

    ## Warning in readLines("../../data/Machinery.fasta"): incomplete final line found
    ## on '../../data/Machinery.fasta'

``` r
headers <- lines[startsWith(lines, ">")]
parsed    <- regmatches(headers, regexec("^>(\\S+)\\s+peptide:\\s+(\\S+)", headers))
label     <- sapply(parsed, `[`, 2)   # e.g. "Dnmt1-201"
eid       <- sapply(parsed, `[`, 3)   # e.g. "ENSMUSP00000004202"
gene_name     <- strip_isoform(label)
protein_name  <- sapply(gene_name, lookup_protein_name, USE.NAMES = FALSE)

machinery_ref <- data.frame(
  accession     = eid,
  protein_name  = protein_name,
  gene_name     = gene_name,
  species       = species_from_id(eid),
  stringsAsFactors = FALSE
)

write.csv(machinery_ref, "../output/20-supplementary-files/Machinery_reference_table.csv", row.names = FALSE)
```

# 2 Epimachinery category labels for Machinery.fasta table

``` r
machinery <- read_csv("../output/20-supplementary-files/Machinery_reference_table.csv",
                      show_col_types = FALSE)

classify_gene <- function(gene) {
  # Uppercase to handle mouse mixed-case gene names
  gene <- toupper(gene)

  cats <- list(
    "ADP-ribosylation" = list(
      prefix = c("PARP", "TNKS", "TARG", "MACROD", "OARD", "ADPRH", "ARH", "TIPARP"),
      exact  = c("PARG", "ARTD1")
    ),
    "DNA methylation & reading" = list(
      prefix = c("DNMT", "TET", "MBD", "UHRF", "MECP",
                 "ZBTB"),
      exact  = c("PRDM14",
                 "ZFP57")
    ),
    "RNA modification" = list(
      prefix = c("METTL", "YTHDF", "YTHDC", "ALKBH", "WTAP", "VIRMA", "ZC3H13",
                 "IGF2BP", "NSUN", "DKC", "TRMT", "PUS", "NAT10", "ADAR",
                 "TRUB", "RPUSD"),
      exact  = c("FTO", "RBM15", "RBM15B", "PCIF1", "BUD23", "WDR4",
                 "RBMX", "TRDMT1", "NOP2", "HNRNPA2B1",
                 "DUS1L",   # tRNA dihydrouridine synthase
                 "TYW5",    # tRNA wybutosine-synthesizing protein (Supp. folder:
                            # "Histone demethylation" — incorrect; TYW5 is a
                            # tRNA modification enzyme, not a histone demethylase)
                 "CBLL1")   # Hakai, E3 ubiquitin ligase component of the m6A
                            # writer complex
    ),
    "ncRNA biogenesis & silencing" = list(
      prefix = c("AGO", "TNRC6", "DICER", "DROSHA", "DGCR8", "PIWI", "EXPORTIN",
                 "TARBP", "PRKRA", "LIN28", "KHSRP", "ZCCHC", "DIS3", "INTS", "EXOSC",
                 "PAN2", "PAN3", "DDX5", "DDX17", "DDX6", "EDC", "DCP", "XRN", "CNOT",
                 "UPF", "SMG", "TUT", "MOV10", "FMR", "PUM", "STAU", "NCBP", "PHAX",
                 "RBM7",
                 "ZC3H12"),  # ZC3H12A-D are endoribonucleases (Regnase-1-4)
                            # that degrade target mRNAs in the cytoplasm —
                            # mRNA turnover/silencing machinery
      exact  = c("RDRP", "MTREX", "SKIV2L2", "MTR4", "ZFC3H1", "SRRT", "ARS2", "XPO5",
                 "XPO1", "CRM1", "RAN", "EWSR1", "FUS", "TARDBP", "HNRNPK", "PARN",
                 "PABPC1", "PAPD5", "DZIP3", "NONO", "SFPQ", "PSPC1",
                 "YBX1", "HSPA8", "HSP90AA1", "HSP90AB1", "PQBP1", "MCM3AP",
                 "PRPF8",  # pre-mRNA splicing factor 8 — core spliceosome
                           # component involved in snRNA processing
                 "EIF3H",  # eIF3 subunit H — translation initiation factor
                           # involved in mRNA silencing/regulation
                 "EIF3F")  # eIF3 subunit F — translation initiation factor;
                           # has MPN deubiquitinase fold but functions in
                           # translation regulation and miRNA-mediated silencing
    ),
    "Histone modification and variants" = list(
      prefix = c("KMT", "KDM", "SETD", "EZH", "SUV", "EHMT", "MLL", "ASH1", "DOT1",
                 "JMJD", "HDAC", "SIRT", "KAT", "PRMT", "CARM", "NSD", "WHSC", "SMYD",
                 "CDYL", "RIOX", "MINA", "NO66",
                 "PRDM",
                 "SUMO", "SENP",
                 "H1-", "H2A", "H2B", "H3", "H4",  # histone variants
                 "NAT"),  # NAT8, NAT8B, NAT8L — N-acetyltransferases (camello
                         # family); Supp. folder groups all NATs under
                         # "Histone acetylation"
      exact  = c("CREBBP", "EP300", "HAT1", "ELP3", "PHF8", "PHF2", "JARID2", "PKN1",
                 "HSPBAP1",
                 "WDR5", "ASH2L", "SETMAR",
                 "OGT", "OGA",
                 "MACROH2A1",  # histone variant macro-H2A.1
                 "ATAT1",      # alpha-tubulin acetyltransferase; Supp. folder:
                               # "Histone acetylation" (grouped with KATs)
                 "GTF3C1", "GTF3C2",  # TFIIIC subunits; Supp. folder:
                                      # "Histone acetylation" (have HAT activity)
                 "NAA10", "NAA50", "NAA60",  # N-alpha-acetyltransferases;
                                             # Supp. folder: "Histone acetylation"
                 "BRD4",       # bromodomain-containing protein 4; reads
                               # acetylated histones; Supp. folder:
                               # "Histone acetylation"
                 "ESCO1", "ESCO2",  # N-acetyltransferases (cohesin establishment);
                                   # Supp. folder: "Histone acetylation"
                 "DLAT",       # dihydrolipoyllysine-residue acetyltransferase;
                               # Supp. folder: "Histone acetylation"
                 "ATF2",       # cyclic AMP-dependent TF; histone acetyltransferase
                               # activity; Supp. folder: "Histone acetylation"
                 "BLOC1S1",    # biogenesis of lysosome-related organelles;
                               # GCN5L1; Supp. folder: "Histone acetylation"
                 "ACAT1", "ACAT2",  # acetyl-CoA acetyltransferases; Supp. folder:
                                   # "Histone acetylation"
                 "UTY",        # inactive histone demethylase (KDM6 family homolog)
                 "HR",         # lysine demethylase and nuclear receptor corepressor
                 "MAPT",       # microtubule-associated protein tau; has histone
                               # H3 binding / chromatin roles
                 "BAZ2B")      # bromodomain adjacent to zinc finger; chromatin
                               # remodeling reader (Supp. folder: not present,
                               # but BAZ1B is in "Histone phosphorylation")
    ),
    "Ubiquitin signaling" = list(
      prefix = c("USP", "UBE", "UBR", "UBA", "UCHL", "UCH", "RNF", "TRIM", "HUWE",
                 "HERC", "MARCH", "RBX", "CUL", "FBX", "SKP", "SOCS", "OTU", "UBQLN",
                 "PSMD", "JOSD", "DDB", "MIB"),
      exact  = c("VCPIP1", "BAP1", "ATXN3",
                 "ATXN3L",    # ataxin-3-like protein; Josephin domain DUB
                 "COPS6", "COPS5", "PARK2", "ZRANB1",
                 "BRCA1", "RAG1", "TOM1", "DTX3L",
                 "BRCC3",     # Lys-63-specific deubiquitinase; Supp. folder:
                              # "Histone Ubiquitination" (Deubiquitination)
                 "MYSM1",     # histone H2A deubiquitinase; Supp. folder:
                              # "Histone Ubiquitination" (Deubiquitination)
                 "STAMBPL1",  # AMSH-LP; deubiquitinase; Supp. folder:
                              # "Histone Ubiquitination" (Deubiquitination)
                 "MPND",      # MPN domain-containing protein; DUB fold;
                              # associated with ubiquitin signaling
                 "ALG13")     # UDP-N-acetylglucosamine transferase subunit;
                              # has deubiquitinase-like fold; Supp. folder:
                              # not present, but classified here per its
                              # role in ubiquitin-related pathways
    ),
    "Chromatin signaling" = list(
      prefix = c("MAP3K", "MAP2K", "MAPK",
                 "EYA",
                 "PPP4", "PPP6"),
      exact  = c("CHUK", "IKBKB", "IKBKG",
                 "PPP1CB", "PPP1CA", "PPP1CC", "PPP1R7",
                 "PPP1R10", "PPP1R8",
                 "PPP2CA", "AURKB", "MSK1", "MSK2", "RPS6KA4", "RPS6KA5", "PARK7",
                 "BUB1", "DLK1", "DAPK3", "PKM", "JAK2", "DBF4", "DYRK1A",
                 "YWHAQ", "PRKCB", "PRKCD", "AURKA", "AURKC", "STK4", "MST1",
                 "MDC1", "MCPH1", "PAK2", "HIRA", "ATM", "ATR", "CSNK2B",
                 "PRKAA1", "PRKAA2", "PRKDC", "PIM1", "CDK8", "CDK1", "CHEK1",
                 "PPM1D", "PPP6C", "PPP4C", "CDCA2",
                 "APBB1", "BAZ1B", "JKAMP", "ARID1A", "ARID1B",
                 "SLK",       # STE20-like kinase; Supp. folder:
                              # "Histone phosphorylation" (slkb)
                 "WEE1",      # Wee1-like kinase; Supp. folder:
                              # "Histone phosphorylation"
                 "HASPIN",    # histone H3T3 kinase (GSG2); chromatin signaling
                 "SGO1", "SGO2")  # Shugoshin proteins; protect centromeric
                                  # cohesion by recruiting PP2A; chromatin-
                                  # associated signaling at kinetochores
    )
  )

  for (cat_name in names(cats)) {
    cat_def <- cats[[cat_name]]
    if (gene %in% cat_def$exact) return(cat_name)
    for (p in cat_def$prefix) {
      if (str_starts(gene, fixed(p))) return(cat_name)
    }
  }
  return("Unclassified")
}

# Apply classification and add category column
machinery_categorized <- machinery %>%
  mutate(epimachinery_category = vapply(gene_name, classify_gene, character(1)))

machinery_categorized
```

    ## # A tibble: 499 × 5
    ##    accession          protein_name       gene_name species epimachinery_category
    ##    <chr>              <chr>              <chr>     <chr>   <chr>                
    ##  1 ENSMUSP00000004202 DNA (cytosine-5-)… Dnmt1     Mus mu… DNA methylation & re…
    ##  2 ENSMUSP00000020991 DNA (cytosine-5-)… Dnmt3a    Mus mu… DNA methylation & re…
    ##  3 ENSMUSP00000051830 DNA (cytosine-5-)… Dnmt3b    Mus mu… DNA methylation & re…
    ##  4 ENSMUSP00000000746 DNA (cytosine-5-)… Dnmt3l    Mus mu… DNA methylation & re…
    ##  5 ENSP00000342812    Ubiquitin carboxy… USP9Y     Homo s… Ubiquitin signaling  
    ##  6 ENSP00000219473    Ubiquitin carboxy… USP10     Homo s… Ubiquitin signaling  
    ##  7 ENSP00000218348    Ubiquitin carboxy… USP11     Homo s… Ubiquitin signaling  
    ##  8 ENSP00000282344    Ubiquitin carboxy… USP12     Homo s… Ubiquitin signaling  
    ##  9 ENSP00000261601    Ubiquitin carboxy… USP14     Homo s… Ubiquitin signaling  
    ## 10 ENSP00000280377    Ubiquitin carboxy… USP15     Homo s… Ubiquitin signaling  
    ## # ℹ 489 more rows

``` r
write_csv(machinery_categorized, "../output/20-supplementary-files/Machinery_categorized.csv")
```

# 3 Alignment summary tables

Compiling per-sample alignment and processing statistics for all three
species

## 3.1 RNA-seq

Set up the metadata

``` r
species_meta <- tibble(
  species = c("D-Apul", "E-Peve", "F-Ptuh"),
  species_full = c("Acropora pulchra", "Porites evermanni",
                   "Pocillopora tuahiniensis"),
  genome = c("Apulcra-genome.fa", "Porites_evermanni_v1.fa",
             "Pocillopora_meandrina_HIv1.assembly.fasta"),
  rna_prefix = c("../../D-Apul", "../../E-Peve", "../../F-Ptuh")
)
species_meta
```

    ## # A tibble: 3 × 4
    ##   species species_full             genome                             rna_prefix
    ##   <chr>   <chr>                    <chr>                              <chr>     
    ## 1 D-Apul  Acropora pulchra         Apulcra-genome.fa                  ../../D-A…
    ## 2 E-Peve  Porites evermanni        Porites_evermanni_v1.fa            ../../E-P…
    ## 3 F-Ptuh  Pocillopora tuahiniensis Pocillopora_meandrina_HIv1.assemb… ../../F-P…

### 3.1.1 1.1 Raw reads

from MultiQC FastQC on raw reads

``` r
raw_reads <- rbind(
  
read_tsv("../../D-Apul/output/01-Apul-RNA-trimming-FastQC/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../E-Peve/output/01-Peve-RNA-trimming-FastQC/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../F-Ptuh/output/01-Ptuh-RNA-trimming-FastQC/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 18 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 20 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 20 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
raw_reads
```

    ## # A tibble: 15 × 2
    ##    Sample                    `Total Sequences`
    ##    <chr>                                 <dbl>
    ##  1 RNA-ACR-140-S1-TP2_R1_001          48199794
    ##  2 RNA-ACR-145-S1-TP2_R1_001          43288869
    ##  3 RNA-ACR-150-S1-TP2_R1_001          44216948
    ##  4 RNA-ACR-173-S1-TP2_R1_001          47976843
    ##  5 RNA-ACR-178-S1-TP2_R1_001          43185440
    ##  6 RNA-POR-71-S1-TP2_R1_001           51248209
    ##  7 RNA-POR-73-S1-TP2_R1_001           51750431
    ##  8 RNA-POR-76-S1-TP2_R1_001           50199897
    ##  9 RNA-POR-79-S1-TP2_R1_001           50299323
    ## 10 RNA-POR-82-S1-TP2_R1_001           49279402
    ## 11 RNA-POC-47-S1-TP2_R1_001           54705696
    ## 12 RNA-POC-48-S1-TP2_R1_001           52036501
    ## 13 RNA-POC-50-S1-TP2_R1_001           55792440
    ## 14 RNA-POC-53-S1-TP2_R1_001           53888583
    ## 15 RNA-POC-57-S1-TP2_R1_001           42934866

### 3.1.2 1.2 Trimmed reads

from MultiQC FastQC on trimmed reads

``` r
trimmed_reads <- rbind(
  
read_tsv("../../D-Apul/output/01-Apul-RNA-trimming-FastQC/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../E-Peve/output/01-Peve-RNA-trimming-FastQC/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("../../F-Ptuh/output/01-Ptuh-RNA-trimming-FastQC/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
trimmed_reads
```

    ## # A tibble: 15 × 2
    ##    Sample                    `Total Sequences`
    ##    <chr>                                 <dbl>
    ##  1 RNA-ACR-140-S1-TP2_R1_001          47710408
    ##  2 RNA-ACR-145-S1-TP2_R1_001          42864294
    ##  3 RNA-ACR-150-S1-TP2_R1_001          43712298
    ##  4 RNA-ACR-173-S1-TP2_R1_001          47501524
    ##  5 RNA-ACR-178-S1-TP2_R1_001          42677752
    ##  6 RNA-POR-71-S1-TP2_R1_001           50831351
    ##  7 RNA-POR-73-S1-TP2_R1_001           51385213
    ##  8 RNA-POR-76-S1-TP2_R1_001           49828147
    ##  9 RNA-POR-79-S1-TP2_R1_001           49976568
    ## 10 RNA-POR-82-S1-TP2_R1_001           48908730
    ## 11 RNA-POC-47-S1-TP2_R1_001           54219233
    ## 12 RNA-POC-48-S1-TP2_R1_001           51579450
    ## 13 RNA-POC-50-S1-TP2_R1_001           55302459
    ## 14 RNA-POC-53-S1-TP2_R1_001           53353608
    ## 15 RNA-POC-57-S1-TP2_R1_001           42513488

### 3.1.3 1.3 Retained after trimming (%)

``` r
read_counts <- raw_reads %>%
  left_join(trimmed_reads, by = c("Sample")) %>%
  mutate(retained_pct = round(trimmed_reads$'Total Sequences' / raw_reads$'Total Sequences' * 100, 2))

read_counts
```

    ## # A tibble: 15 × 4
    ##    Sample                   `Total Sequences.x` `Total Sequences.y` retained_pct
    ##    <chr>                                  <dbl>               <dbl>        <dbl>
    ##  1 RNA-ACR-140-S1-TP2_R1_0…            48199794            47710408         99.0
    ##  2 RNA-ACR-145-S1-TP2_R1_0…            43288869            42864294         99.0
    ##  3 RNA-ACR-150-S1-TP2_R1_0…            44216948            43712298         98.9
    ##  4 RNA-ACR-173-S1-TP2_R1_0…            47976843            47501524         99.0
    ##  5 RNA-ACR-178-S1-TP2_R1_0…            43185440            42677752         98.8
    ##  6 RNA-POR-71-S1-TP2_R1_001            51248209            50831351         99.2
    ##  7 RNA-POR-73-S1-TP2_R1_001            51750431            51385213         99.3
    ##  8 RNA-POR-76-S1-TP2_R1_001            50199897            49828147         99.3
    ##  9 RNA-POR-79-S1-TP2_R1_001            50299323            49976568         99.4
    ## 10 RNA-POR-82-S1-TP2_R1_001            49279402            48908730         99.2
    ## 11 RNA-POC-47-S1-TP2_R1_001            54705696            54219233         99.1
    ## 12 RNA-POC-48-S1-TP2_R1_001            52036501            51579450         99.1
    ## 13 RNA-POC-50-S1-TP2_R1_001            55792440            55302459         99.1
    ## 14 RNA-POC-53-S1-TP2_R1_001            53888583            53353608         99.0
    ## 15 RNA-POC-57-S1-TP2_R1_001            42934866            42513488         99.0

### 3.1.4 1.4 Alignment rate, unique %, multi %

``` r
parse_hisat2_log <- function(lines) {
  text <- paste(lines, collapse = "\n")

  total_pairs <- as.numeric(str_match(text, "^(\\d+) reads")[2])
  total_reads <- total_pairs * 2

  # Overall alignment rate
  overall_rate <- as.numeric(str_match(text, "(\\S+)% overall alignment rate")[2])

  # Concordant pairs 
  conc_0  <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned concordantly 0 times")[2])
  conc_1  <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned concordantly exactly 1 time")[2])
  conc_gt1 <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned concordantly >1 times")[2])

  # Discordant pairs (from concordant-0-times)
  disc_1 <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned discordantly 1 time")[2])
  disc_1 <- if (is.na(disc_1)) 0 else disc_1

  # Individual mates (from pairs that aligned 0 times concordantly AND discordantly)
  mates_total <- as.numeric(str_match(text, "(\\d+) mates make up the pairs")[2])
  mates_total <- if (is.na(mates_total)) 0 else mates_total

  mates_0  <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned 0 times\n")[2])
  mates_0  <- if (is.na(mates_0)) 0 else mates_0

  mates_1  <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned exactly 1 time\n")[2])
  mates_1  <- if (is.na(mates_1)) 0 else mates_1

  mates_gt1 <- as.numeric(str_match(text, "(\\d+) \\(\\S+%\\) aligned >1 times\n")[2])
  mates_gt1 <- if (is.na(mates_gt1)) 0 else mates_gt1

  # Compute READ-level unique/multi (combining all alignment levels)
  # Each concordant pair = 2 reads; each discordant pair = 2 reads
  unique_reads <- (conc_1 * 2) + (disc_1 * 2) + mates_1
  multi_reads  <- (conc_gt1 * 2) + mates_gt1
  unmapped_reads <- mates_0
  mapped_reads <- unique_reads + multi_reads
  mapped_read_pairs <- mapped_reads/2

  # Sanity check: unique + multi + unmapped should ≈ total_reads
  # (may not be exact due to mates from discordant pairs that also aligned individually)

  tibble(
    total_pairs = total_pairs,
    total_reads = total_reads,
    overall_align_rate = overall_rate,
    unique_pct = round(unique_reads / total_reads * 100, 2),
    multi_pct = round(multi_reads / total_reads * 100, 2),
    mapped_read_pairs = mapped_read_pairs
  )
}

# D-Apul: consolidated log (all 5 samples in one file)
parse_hisat2_consolidated <- function(path, samples) {
  lines <- readLines(path)
  block_starts <- which(str_detect(lines, "^\\d+ reads; of these:"))
  results <- list()
  for (i in seq_along(block_starts)) {
    end <- if (i < length(block_starts)) block_starts[i + 1] - 1 else length(lines)
    block <- lines[block_starts[i]:end]
    results[[i]] <- parse_hisat2_log(block) %>%
      mutate(sample = samples[i])
  }
  bind_rows(results)
}

# E-Peve / F-Ptuh: one stderr file per sample
parse_hisat2_per_file <- function(dir_path) {
  files <- list.files(dir_path, pattern = "_hisat\\.stderr$", full.names = TRUE)
  results <- list()
  for (f in files) {
    sname <- basename(f) %>%
      str_remove("_hisat.stderr") %>%
      str_remove("^RNA-")
    results[[sname]] <- parse_hisat2_log(readLines(f)) %>%
      mutate(sample = sname)
  }
  bind_rows(results)
}

hisat2_stats <- bind_rows(
  parse_hisat2_consolidated(
    "../../D-Apul/output/07-Apul-Hisat/hisat.out",
    samples = c("ACR-140", "ACR-145", "ACR-150", "ACR-173", "ACR-178")
  ),
  parse_hisat2_per_file("../../E-Peve/output/06-Peve-Hisat"),
  parse_hisat2_per_file("../../F-Ptuh/output/06-Ptuh-Hisat")
) %>%
  mutate(species = case_when(
    str_detect(sample, "^ACR") ~ "D-Apul",
    str_detect(sample, "^POR") ~ "E-Peve",
    str_detect(sample, "^POC") ~ "F-Ptuh"
  )) %>%
  select(species, sample, everything())

hisat2_stats
```

    ## # A tibble: 15 × 8
    ##    species sample  total_pairs total_reads overall_align_rate unique_pct
    ##    <chr>   <chr>         <dbl>       <dbl>              <dbl>      <dbl>
    ##  1 D-Apul  ACR-140    47710408    95420816               69.3       65.6
    ##  2 D-Apul  ACR-145    42864294    85728588               71.7       66.0
    ##  3 D-Apul  ACR-150    43712298    87424596               54.3       48.4
    ##  4 D-Apul  ACR-173    47501524    95003048               67.4       62.9
    ##  5 D-Apul  ACR-178    42677752    85355504               71.1       64.7
    ##  6 E-Peve  POR-71     50831351   101662702               86.4       67.0
    ##  7 E-Peve  POR-73     51385213   102770426               87.1       63.6
    ##  8 E-Peve  POR-76     49828147    99656294               86.0       65.9
    ##  9 E-Peve  POR-79     49976568    99953136               89.1       63.6
    ## 10 E-Peve  POR-82     48908730    97817460               84.2       66.2
    ## 11 F-Ptuh  POC-47     54219233   108438466               64.2       54.7
    ## 12 F-Ptuh  POC-48     51579450   103158900               60.9       51.0
    ## 13 F-Ptuh  POC-50     55302459   110604918               62.0       56.6
    ## 14 F-Ptuh  POC-53     53353608   106707216               60.6       56.1
    ## 15 F-Ptuh  POC-57     42513488    85026976               50.6       47.6
    ## # ℹ 2 more variables: multi_pct <dbl>, mapped_read_pairs <dbl>

### 3.1.5 1.5 Mismatch %

This requires collecting info from the full \*.sorted.bam files, which
are very large and live on Gannet
(<https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/>).

``` bash

# The files are backed up to Steven's backup directory: 
#     /volume/web/seashell/bu-github/deep-dive-expression

# The BAM paths are:
#   D-Apul: D-Apul/output/07-Apul-Hisat/RNA-ACR-{140,145,150,173,178}.sorted.bam
#   E-Peve: E-Peve/output/06-Peve-Hisat/RNA-POR-{71,73,76,79,82}.sorted.bam
#   F-Ptuh: F-Ptuh/output/06-Ptuh-Hisat/RNA-POC-{47,48,50,53,57}.sorted.bam

mkdir -p ../output/20-supplementary-files/RNAseq_HISAT_samtools_stats

# D-Apul BAMs
wget -r -l1 -nd -nc -A "*.sorted.bam" \
  -P ../output/20-supplementary-files/RNAseq_HISAT_samtools_stats \
  https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/D-Apul/output/07-Apul-Hisat/

# E-Peve BAMs
wget -r -l1 -nd -nc -A "*.sorted.bam" \
  -P ../output/20-supplementary-files/RNAseq_HISAT_samtools_stats \
  https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/E-Peve/output/06-Peve-Hisat/

# F-Ptuh BAMs
wget -r -l1 -nd -nc -A "*.sorted.bam" \
  -P ../output/20-supplementary-files/RNAseq_HISAT_samtools_stats \
  https://gannet.fish.washington.edu/seashell/bu-github/deep-dive-expression/F-Ptuh/output/06-Ptuh-Hisat/

# Then run samtools on each BAM
for bam in ../output/20-supplementary-files/RNAseq_HISAT_samtools_stats/*.sorted.bam; do
  base=$(basename "$bam" .sorted.bam)
  /home/shared/samtools-1.12/samtools stats "$bam" > "../output/20-supplementary-files/RNAseq_HISAT_samtools_stats/${base}.stats.txt"
  /home/shared/samtools-1.12/samtools flagstat "$bam" > "../output/20-supplementary-files/RNAseq_HISAT_samtools_stats/${base}.flagstat.txt"
done
```

``` bash
# Remove BAMs once done (they're quite large)
rm ../output/20-supplementary-files/RNAseq_HISAT_samtools_stats/*.sorted.bam
```

Parse samtools stats output for mismatch rate

``` r
samtools_dir <- "../output/20-supplementary-files/RNAseq_HISAT_samtools_stats"

# Sample configuration
samples_cfg <- list(
  D_Apul = c("ACR-140", "ACR-145", "ACR-150", "ACR-173", "ACR-178"),
  E_Peve = c("POR-71", "POR-73", "POR-76", "POR-79", "POR-82"),
  F_Ptuh = c("POC-47", "POC-48", "POC-50", "POC-53", "POC-57")
)
species_map <- c(D_Apul = "D-Apul", E_Peve = "E-Peve", F_Ptuh = "F-Ptuh")

# Parse a single samtools stats file
parse_samtools_stats <- function(path, sample_name, species_code) {
  if (!file.exists(path)) {
    return(tibble(species = species_code, sample = sample_name,
                  mapped_reads = NA, mismatch_rate_pct = NA))
  }
  lines <- readLines(path)

  # Mapped reads: SN\treads mapped:\t28412345
  mapped <- str_match(lines[grepl("^SN\\treads mapped:", lines)],
                      "^SN\\treads mapped:\\t(\\d+)")[2]
  mapped <- if (!is.na(mapped)) as.numeric(mapped) else NA

  # Mismatch rate (reported as "error rate"):
  # SN\terror rate:\t7.151265e-03\t# mismatches / bases mapped (cigar)
  error_line <- lines[grepl("^SN\\terror rate:", lines)]
  error_rate <- str_match(error_line, "^SN\\terror rate:\\t([\\d.eE+-]+)")[2]
  mismatch_rate_pct <- if (!is.na(error_rate)) as.numeric(error_rate) * 100 else NA

  tibble(species = species_code, sample = sample_name,
         mapped_reads = mapped,
         mismatch_rate_pct = round(mismatch_rate_pct, 4))
}

# Build results from single directory
samtools_results <- bind_rows(lapply(names(samples_cfg), function(sp) {
  species_code <- species_map[sp]
  bind_rows(lapply(samples_cfg[[sp]], function(s) {
    parse_samtools_stats(
      file.path(samtools_dir, paste0("RNA-", s, ".stats.txt")),
      s, species_code
    )
  }))
}))

# Rename for final output
samtools_summary <- samtools_results %>%
  rename(
    Species = species,
    Sample = sample,
    `Mapped reads` = mapped_reads,
    `Mismatch rate (%)` = mismatch_rate_pct
  )

samtools_summary
```

    ## # A tibble: 15 × 4
    ##    Species Sample  `Mapped reads` `Mismatch rate (%)`
    ##    <chr>   <chr>            <dbl>               <dbl>
    ##  1 D-Apul  ACR-140       66166297               0.715
    ##  2 D-Apul  ACR-145       61488934               0.695
    ##  3 D-Apul  ACR-150       47452375               0.668
    ##  4 D-Apul  ACR-173       64013149               0.745
    ##  5 D-Apul  ACR-178       60684881               0.760
    ##  6 E-Peve  POR-71        87835896               0.338
    ##  7 E-Peve  POR-73        89515483               0.268
    ##  8 E-Peve  POR-76        85668912               0.293
    ##  9 E-Peve  POR-79        89024198               0.226
    ## 10 E-Peve  POR-82        82338762               0.351
    ## 11 F-Ptuh  POC-47        69558951               0.954
    ## 12 F-Ptuh  POC-48        62823383               0.928
    ## 13 F-Ptuh  POC-50        68623296               1.08 
    ## 14 F-Ptuh  POC-53        64638718               1.10 
    ## 15 F-Ptuh  POC-57        43011432               1.16

### 3.1.6 1.6 Number of lncRNAs identified

``` r
count_lncRNAs <- function(filepath) {
  
df <- read_tsv(filepath, col_names = TRUE)
sample_cols <- grep(".sorted.bam$", colnames(df), value = TRUE)
sample_names <- str_remove(sample_cols, ".*pipeline\\.RNA\\.") %>% str_remove("\\.sorted\\.bam$")

# Count non-zero entries per sample column
nonzero_counts <- sapply(sample_cols, function(col) {
  sum(df[[col]] > 0, na.rm = TRUE)
})

# Build a summary tibble
summary <- tibble(
  sample = sample_names,
  n_unique_lncRNA = nonzero_counts
)
print(summary)
}

lncRNA_counts <- rbind(
count_lncRNAs("../output/01.6-lncRNA-pipeline/Apul-lncRNA-counts-filtered.txt"),
count_lncRNAs("../output/01.6-lncRNA-pipeline/Peve-lncRNA-counts-filtered.txt"),
count_lncRNAs("../output/01.6-lncRNA-pipeline/Ptuh-lncRNA-counts-filtered.txt")
)
```

    ## Rows: 31490 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Geneid, Chr, Strand
    ## dbl (8): Start, End, Length, X.home.shared.8TB_HDD_02.zbengt.github.deep.div...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 5 × 2
    ##   sample  n_unique_lncRNA
    ##   <chr>             <int>
    ## 1 ACR.140           27515
    ## 2 ACR.145           21807
    ## 3 ACR.150           28248
    ## 4 ACR.173           24585
    ## 5 ACR.178           29713

    ## Rows: 10085 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Geneid, Chr, Strand
    ## dbl (8): Start, End, Length, X.home.shared.8TB_HDD_02.zbengt.github.deep.div...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 5 × 2
    ##   sample n_unique_lncRNA
    ##   <chr>            <int>
    ## 1 POR.71            8577
    ## 2 POR.73            7763
    ## 3 POR.76            9498
    ## 4 POR.79            8007
    ## 5 POR.82            8681

    ## Rows: 16152 Columns: 11
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (3): Geneid, Chr, Strand
    ## dbl (8): Start, End, Length, ...output.01.6.Ptuh.lncRNA.pipeline.RNA.POC.47....
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

    ## # A tibble: 5 × 2
    ##   sample        n_unique_lncRNA
    ##   <chr>                   <int>
    ## 1 POC.47.S1.TP2           12017
    ## 2 POC.48.S1.TP2           13866
    ## 3 POC.50.S1.TP2           14294
    ## 4 POC.53.S1.TP2           14409
    ## 5 POC.57.S1.TP2           12197

``` r
lncRNA_counts
```

    ## # A tibble: 15 × 2
    ##    sample        n_unique_lncRNA
    ##    <chr>                   <int>
    ##  1 ACR.140                 27515
    ##  2 ACR.145                 21807
    ##  3 ACR.150                 28248
    ##  4 ACR.173                 24585
    ##  5 ACR.178                 29713
    ##  6 POR.71                   8577
    ##  7 POR.73                   7763
    ##  8 POR.76                   9498
    ##  9 POR.79                   8007
    ## 10 POR.82                   8681
    ## 11 POC.47.S1.TP2           12017
    ## 12 POC.48.S1.TP2           13866
    ## 13 POC.50.S1.TP2           14294
    ## 14 POC.53.S1.TP2           14409
    ## 15 POC.57.S1.TP2           12197

## 3.2 sRNA-seq alignment statistics

### 3.2.1 2.1 Species metadata (genome)

``` r
# Genome files are referenced in the ShortStack Rmd scripts:
srna_genomes <- tibble(
  species = c("D-Apul", "E-Peve", "F-Ptuh"),
  genome = c("Apulcra-genome.fa", "Porites_evermanni_v1.fa",
             "Pocillopora_meandrina_HIv1.assembly.fasta")
)
srna_genomes
```

    ## # A tibble: 3 × 2
    ##   species genome                                   
    ##   <chr>   <chr>                                    
    ## 1 D-Apul  Apulcra-genome.fa                        
    ## 2 E-Peve  Porites_evermanni_v1.fa                  
    ## 3 F-Ptuh  Pocillopora_meandrina_HIv1.assembly.fasta

### 3.2.2 2.1 Raw reads

``` r
s_raw_reads <- rbind(
  
read_tsv("~/deep-dive/D-Apul/output/08-Apul-sRNAseq-trimming/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("~/deep-dive/E-Peve/output/06-Peve-sRNAseq-trimming/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences'), 

read_tsv("~/deep-dive/F-Pmea/output/08-Pmea-sRNAseq-trimming/raw-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "RNA-")) %>% filter(str_detect(Sample, "_R1_001$")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 6 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 10 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (16): Sample, Filename, File type, Encoding, Total Bases, basic_statisti...
    ## dbl  (7): Total Sequences, Sequences flagged as poor quality, Sequence lengt...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
s_raw_reads$Sample <- str_remove(s_raw_reads$Sample, "_R1_001")

s_raw_reads
```

    ## # A tibble: 13 × 2
    ##    Sample              `Total Sequences`
    ##    <chr>                           <dbl>
    ##  1 sRNA-ACR-140-S1-TP2          18506972
    ##  2 sRNA-ACR-145-S1-TP2          21295979
    ##  3 sRNA-ACR-150-S1-TP2          22586680
    ##  4 sRNA-ACR-173-S1-TP2          19686271
    ##  5 sRNA-ACR-178-S1-TP2          17687362
    ##  6 sRNA-POR-73-S1-TP2           22960890
    ##  7 sRNA-POR-79-S1-TP2           19762559
    ##  8 sRNA-POR-82-S1-TP2           20675792
    ##  9 sRNA-POC-47-S1-TP2           27579504
    ## 10 sRNA-POC-48-S1-TP2           24884997
    ## 11 sRNA-POC-50-S1-TP2           26073937
    ## 12 sRNA-POC-53-S1-TP2           25726329
    ## 13 sRNA-POC-57-S1-TP2           26459749

### 3.2.3 2.2 Trimmed reads

``` r
s_trimmed_reads <- rbind(
  
read_tsv("~/deep-dive/D-Apul/output/08.2-Apul-sRNAseq-trimming-31bp-fastp-merged/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "sRNA-")) %>% select(Sample, 'Total Sequences'), 

read_tsv("~/deep-dive/E-Peve/output/06.2-Peve-sRNAseq-trimming-31bp-fastp-merged/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% select(Sample, 'Total Sequences') %>% mutate(Sample = paste0("sRNA-", Sample)), 

read_tsv("~/deep-dive/F-Pmea/output/08.2-Pmea-sRNAseq-trimming-31bp-fastp-merged/trimmed-fastqc/multiqc_data/multiqc_fastqc.txt", col_names = TRUE) %>% filter(str_detect(Sample, "sRNA-")) %>% select(Sample, 'Total Sequences')
)
```

    ## Rows: 5 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 3 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## Rows: 5 Columns: 23
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr (17): Sample, Filename, File type, Encoding, Total Bases, Sequence lengt...
    ## dbl  (6): Total Sequences, Sequences flagged as poor quality, %GC, total_ded...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
s_trimmed_reads$Sample <- str_remove(s_trimmed_reads$Sample, "-fastp-adapters-polyG-31bp-merged")

s_trimmed_reads
```

    ## # A tibble: 13 × 2
    ##    Sample              `Total Sequences`
    ##    <chr>                           <dbl>
    ##  1 sRNA-ACR-140-S1-TP2          15076001
    ##  2 sRNA-ACR-145-S1-TP2          17265123
    ##  3 sRNA-ACR-150-S1-TP2          18395038
    ##  4 sRNA-ACR-173-S1-TP2          16332044
    ##  5 sRNA-ACR-178-S1-TP2          14295261
    ##  6 sRNA-POR-73-S1-TP2           12238952
    ##  7 sRNA-POR-79-S1-TP2           12259576
    ##  8 sRNA-POR-82-S1-TP2           13933812
    ##  9 sRNA-POC-47-S1-TP2           12066776
    ## 10 sRNA-POC-48-S1-TP2           14045042
    ## 11 sRNA-POC-50-S1-TP2           12141955
    ## 12 sRNA-POC-53-S1-TP2           15732665
    ## 13 sRNA-POC-57-S1-TP2           14136486

### 3.2.4 2.3 Retained after trimming (%)

``` r
s_read_counts <- s_raw_reads %>%
  left_join(s_trimmed_reads, by = c("Sample")) %>%
  mutate(retained_pct = round(s_trimmed_reads$'Total Sequences' / s_raw_reads$'Total Sequences' * 100, 2))

s_read_counts
```

    ## # A tibble: 13 × 4
    ##    Sample              `Total Sequences.x` `Total Sequences.y` retained_pct
    ##    <chr>                             <dbl>               <dbl>        <dbl>
    ##  1 sRNA-ACR-140-S1-TP2            18506972            15076001         81.5
    ##  2 sRNA-ACR-145-S1-TP2            21295979            17265123         81.1
    ##  3 sRNA-ACR-150-S1-TP2            22586680            18395038         81.4
    ##  4 sRNA-ACR-173-S1-TP2            19686271            16332044         83.0
    ##  5 sRNA-ACR-178-S1-TP2            17687362            14295261         80.8
    ##  6 sRNA-POR-73-S1-TP2             22960890            12238952         53.3
    ##  7 sRNA-POR-79-S1-TP2             19762559            12259576         62.0
    ##  8 sRNA-POR-82-S1-TP2             20675792            13933812         67.4
    ##  9 sRNA-POC-47-S1-TP2             27579504            12066776         43.8
    ## 10 sRNA-POC-48-S1-TP2             24884997            14045042         56.4
    ## 11 sRNA-POC-50-S1-TP2             26073937            12141955         46.6
    ## 12 sRNA-POC-53-S1-TP2             25726329            15732665         61.2
    ## 13 sRNA-POC-57-S1-TP2             26459749            14136486         53.4

### 3.2.5 2.4 ShortStack alignment stats

from alignment_details.tsv

NOTE: the Github repo README.md states that the H tag means “H: Very
highly multi-mapped read (\>=50 hits)”. However,the log from running
ShortStack 4.1.0 (e.g., the file
~/deep-dive-expression/D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/shortstack.log)
states that the tag represents “Very highly multi-mapped (\>=20
hits)(H)”. It looks like the \>=50 threshold was used in previous
versions of ShortStack (confirmed by looking at deep-dive ShortStack
logs), but the updated 4.1.0 version switched to a \>=20 threshold. The
Github repo’s README is out of date.

``` r
# alignment_details.tsv columns:
#   readfile  mapping_type  sequence_length  sequence_count  read_count
# mapping_type (from https://github.com/MikeAxtell/ShortStack)
# U: Uniquely mapped (not a multimapper).
# P: Multimapper placed using the method set by option --mmap.
# R: Multimapper placed at random.
# H: Very highly multi-mapped read (>=20 hits).
# N: Unmapped reads.

read_alignment_details <- function(path, species_code) {
  df <- read_tsv(path, show_col_types = FALSE)

  # Extract sample name from readfile path
  df <- df %>%
    mutate(
      sample = str_extract(readfile, "(ACR|POR|POC)-\\d+"),
      species = species_code
    )

  # Aggregate by sample × mapping_type, sum sequence_count
  agg <- df %>%
    group_by(species, sample, mapping_type) %>%
    summarise(seq_count = sum(read_count), .groups = "drop")

  # Pivot to wide format
  agg_wide <- agg %>%
    pivot_wider(names_from = mapping_type, values_from = seq_count,
                values_fill = 0)

  # totals and percentages
  # U=unique, P=guided, R=random, H=very-high, N=not mapped
  agg_wide %>%
    mutate(
      total_reads = U + P + R + H + N,
      mapped_reads = U + P + R + H,
      overall_align_pct = round((mapped_reads) / total_reads * 100, 2),
      unique_pct = round(U / total_reads * 100, 2),
      guided_pct = round(P / total_reads * 100, 2),
      random_pct = round(R / total_reads * 100, 2),
      very_high_pct = round(H / total_reads * 100, 2),
    ) %>%
    select(
      species, sample, total_reads, overall_align_pct, unique_pct, guided_pct, 
      random_pct, very_high_pct, mapped_reads
    ) %>%
    rename(
      `Overall alignment (%)` = overall_align_pct,
      `Uniquely mapped reads (%)` = unique_pct,
      `Multi-mapped reads placed with guidance (%)` = guided_pct,
      `Multi-mapped reads placed randomly (%)` = random_pct,
      `Very high multi-mapping reads (>=50 hits; %)` = very_high_pct,
      `Mapped sRNA reads` = mapped_reads
    )
}

srna_align <- bind_rows(
  read_alignment_details("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/alignment_details.tsv", "D-Apul"),
  read_alignment_details("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/alignment_details.tsv", "E-Peve"),
  read_alignment_details("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/alignment_details.tsv", "F-Ptuh")
)
srna_align
```

    ## # A tibble: 13 × 9
    ##    species sample  total_reads `Overall alignment (%)` Uniquely mapped reads (…¹
    ##    <chr>   <chr>         <dbl>                   <dbl>                     <dbl>
    ##  1 D-Apul  ACR-140    15076180                    80.0                      33.2
    ##  2 D-Apul  ACR-145    17265150                    78.8                      28.9
    ##  3 D-Apul  ACR-150    18395068                    77.8                      35.0
    ##  4 D-Apul  ACR-173    16332641                    78.1                      31.1
    ##  5 D-Apul  ACR-178    14295278                    80.1                      31.9
    ##  6 E-Peve  POR-73     12238972                    73.6                      30.7
    ##  7 E-Peve  POR-79     12259629                    79.3                      40.3
    ##  8 E-Peve  POR-82     13933829                    81.4                      43.7
    ##  9 F-Ptuh  POC-47     12068141                    62.4                      17.8
    ## 10 F-Ptuh  POC-48     14046678                    62.4                      12.9
    ## 11 F-Ptuh  POC-50     12142097                    55.9                      16.2
    ## 12 F-Ptuh  POC-53     15732765                    75.8                      18.3
    ## 13 F-Ptuh  POC-57     14137264                    83.2                      16.8
    ## # ℹ abbreviated name: ¹​`Uniquely mapped reads (%)`
    ## # ℹ 4 more variables: `Multi-mapped reads placed with guidance (%)` <dbl>,
    ## #   `Multi-mapped reads placed randomly (%)` <dbl>,
    ## #   `Very high multi-mapping reads (>=50 hits; %)` <dbl>,
    ## #   `Mapped sRNA reads` <dbl>

### 3.2.6 2.5 Mapped miRNA reads

from ShortStack Counts.txt

``` r
# Counts.txt is a TSV: Coords, Name, MIRNA, <sample1>_condensed, ...
# Filter rows where MIRNA == "Y", then sum each sample column.

read_mirna_counts <- function(path, species_code) {
  df <- read_tsv(path, show_col_types = FALSE)
  sample_cols <- colnames(df)[-(1:3)]  # skip Coords, Name, MIRNA

  # Extract clean sample names from column names
  # e.g. "sRNA-ACR-140-S1-TP2-..._condensed" -> "ACR-140"
  samples <- str_extract(sample_cols, "(ACR|POR|POC)-\\d+")

  mirna_only <- df %>% filter(MIRNA == "Y")

  # Sum each sample column across all miRNA rows
  mapped_mirna <- colSums(mirna_only[, sample_cols])

  tibble(
    species = species_code,
    sample = samples,
    mapped_mirna_reads = mapped_mirna
  )
}

srna_mirna_reads <- bind_rows(
  read_mirna_counts("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Counts.txt", "D-Apul"),
  read_mirna_counts("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt", "E-Peve"),
  read_mirna_counts("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt", "F-Ptuh")
)
srna_mirna_reads
```

    ## # A tibble: 13 × 3
    ##    species sample  mapped_mirna_reads
    ##    <chr>   <chr>                <dbl>
    ##  1 D-Apul  ACR-140             226397
    ##  2 D-Apul  ACR-145             444883
    ##  3 D-Apul  ACR-150             668500
    ##  4 D-Apul  ACR-173             373547
    ##  5 D-Apul  ACR-178             455227
    ##  6 E-Peve  POR-73              373247
    ##  7 E-Peve  POR-79              279129
    ##  8 E-Peve  POR-82              399594
    ##  9 F-Ptuh  POC-47              623726
    ## 10 F-Ptuh  POC-48              368009
    ## 11 F-Ptuh  POC-50              512651
    ## 12 F-Ptuh  POC-53              430756
    ## 13 F-Ptuh  POC-57              178801

### 3.2.7 2.6 Number of miRNAs identified

``` r
count_expressed_mirnas <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE)

  sample_cols <- grep("merged_condensed", colnames(df), value = TRUE)

  as.data.frame(t(df %>%
    filter(MIRNA == "Y") %>%
    summarise(across(all_of(sample_cols), ~ sum(.x > 0, na.rm = TRUE)))))
}

miRNA_counts <- rbind(
count_expressed_mirnas("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Counts.txt"),
count_expressed_mirnas("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt"),
count_expressed_mirnas("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/Counts.txt")
)
colnames(miRNA_counts) <- c("miRNA_count")
miRNA_counts$Sample <- rownames(miRNA_counts)
miRNA_counts <- miRNA_counts %>% select(Sample, miRNA_count) %>% remove_rownames()
miRNA_counts$Sample <- miRNA_counts$Sample %>% str_remove("-S1-TP2-fastp-adapters-polyG-31bp-merged_condensed") %>% str_remove("sRNA-")

miRNA_counts
```

    ##     Sample miRNA_count
    ## 1  ACR-140          38
    ## 2  ACR-145          39
    ## 3  ACR-150          39
    ## 4  ACR-173          39
    ## 5  ACR-178          38
    ## 6   POR-73          44
    ## 7   POR-79          42
    ## 8   POR-82          43
    ## 9   POC-47          37
    ## 10  POC-48          37
    ## 11  POC-50          37
    ## 12  POC-53          37
    ## 13  POC-57          37

# 4 Mature miRNA names and sequences for all species

Already have per-species, just need to combine into single table

``` r
miRNA_details <- rbind(
  read.csv("../../D-Apul/output/11-Apul-sRNA-ShortStack_4.1.0-pulchra_genome/ShortStack_out/Apul_Results_mature_named_miRNAs.csv", header = TRUE) %>% dplyr::select(-X),
  # Including an NA filter here, because it looks like the Peve file has an extra row indicating miRNAs identified in 
  # Ashey et al 2026 which were not also identified here using the new version of ShortStack (4.1.0)
  read.csv("../../E-Peve/output/05-Peve-sRNA-ShortStack_4.1.0/ShortStack_out/Peve_Results_mature_named_miRNAs.csv", header = TRUE) %>% dplyr::select(-X) %>% drop_na(Name),
  read.csv("../../F-Ptuh/output/05-Ptuh-sRNA-ShortStack_4.1.0/ShortStack_out/Ptuh_Results_mature_named_miRNAs.csv", header = TRUE) %>% dplyr::select(-X)
)

miRNA_details
```

    ##                                                        Locus          Name
    ## 1                                 ntLink_6:13351746-13351842  Cluster_1951
    ## 2                                   ntLink_6:4847443-4847535  Cluster_1826
    ## 3                                   ntLink_6:5157537-5157626  Cluster_1832
    ## 4                                   ntLink_6:7263486-7263580  Cluster_1862
    ## 5                               ptg000001l:20063072-20063168  Cluster_2859
    ## 6                                 ptg000001l:5548841-5548934  Cluster_2463
    ## 7                               ptg000002l:14046235-14046328  Cluster_3366
    ## 8                               ptg000002l:14046541-14046634  Cluster_3367
    ## 9                                 ptg000002l:7337508-7337601  Cluster_3250
    ## 10                                ptg000004l:1859862-1859953  Cluster_3437
    ## 11                                ptg000007l:3377283-3377376  Cluster_4254
    ## 12                                  ptg000007l:915905-916002  Cluster_4220
    ## 13                              ptg000008l:10754767-10754862  Cluster_5012
    ## 14                                ptg000009l:4940515-4940609  Cluster_5981
    ## 15                                  ptg000009l:616285-616377  Cluster_5899
    ## 16                                  ptg000009l:616579-616673  Cluster_5900
    ## 17                              ptg000016l:11751358-11751448 Cluster_10093
    ## 18                                ptg000016l:7795508-7795602 Cluster_10051
    ## 19                                ptg000016l:8599833-8599925 Cluster_10057
    ## 20                                ptg000017l:6374963-6375054 Cluster_10207
    ## 21                                ptg000017l:7471118-7471210 Cluster_10228
    ## 22                                ptg000018l:2286780-2286870 Cluster_10419
    ## 23                              ptg000023l:18708903-18708996 Cluster_14402
    ## 24                              ptg000023l:37965243-37965338 Cluster_14768
    ## 25                                ptg000024l:4086202-4086295 Cluster_15316
    ## 26                                ptg000024l:5256454-5256549 Cluster_15340
    ## 27                              ptg000025l:10501030-10501125 Cluster_15851
    ## 28                              ptg000025l:10668901-10668998 Cluster_15854
    ## 29                              ptg000025l:15156415-15156507 Cluster_15928
    ## 30                                ptg000025l:1981692-1981783 Cluster_15671
    ## 31                                ptg000025l:7472531-7472623 Cluster_15775
    ## 32                                ptg000026l:8745540-8745635 Cluster_16409
    ## 33                                ptg000031l:5461305-5461397 Cluster_17776
    ## 34                                ptg000031l:6751935-6752029 Cluster_17791
    ## 35                                ptg000035l:4339550-4339641 Cluster_18711
    ## 36                                ptg000035l:4808338-4808432 Cluster_18723
    ## 37                                ptg000035l:5346032-5346127 Cluster_18728
    ## 38                                ptg000035l:8367726-8367821 Cluster_18772
    ## 39                                    ptg000039l:35764-35856 Cluster_19193
    ## 40              Porites_evermani_scaffold_1128:125526-125618  Cluster_9983
    ## 41                  Porites_evermani_scaffold_1159:6653-6743 Cluster_10060
    ## 42                  Porites_evermani_scaffold_1159:7406-7500 Cluster_10061
    ## 43               Porites_evermani_scaffold_138:127944-128038  Cluster_2787
    ## 44                Porites_evermani_scaffold_1415:72049-72142 Cluster_10934
    ## 45                Porites_evermani_scaffold_1429:47235-47326 Cluster_10965
    ## 46                   Porites_evermani_scaffold_145:5384-5473  Cluster_2854
    ## 47                 Porites_evermani_scaffold_148:94531-94624  Cluster_2882
    ## 48                Porites_evermani_scaffold_1503:47558-47651 Cluster_11134
    ## 49                Porites_evermani_scaffold_1503:47749-47842 Cluster_11135
    ## 50                Porites_evermani_scaffold_16:383386-383478   Cluster_589
    ## 51                Porites_evermani_scaffold_1732:76482-76574 Cluster_11997
    ## 52               Porites_evermani_scaffold_1:1404250-1404342    Cluster_29
    ## 53                Porites_evermani_scaffold_2259:40195-40285 Cluster_13502
    ## 54               Porites_evermani_scaffold_253:200943-201039  Cluster_4079
    ## 55               Porites_evermani_scaffold_253:202235-202327  Cluster_4080
    ## 56               Porites_evermani_scaffold_257:110308-110402  Cluster_4115
    ## 57                Porites_evermani_scaffold_26:382550-382645   Cluster_796
    ## 58                Porites_evermani_scaffold_2738:56836-56926 Cluster_14500
    ## 59                Porites_evermani_scaffold_3072:29428-29521 Cluster_14999
    ## 60                 Porites_evermani_scaffold_316:88415-88506  Cluster_4629
    ## 61               Porites_evermani_scaffold_334:153554-153646  Cluster_4735
    ## 62                Porites_evermani_scaffold_3707:36102-36199 Cluster_15726
    ## 63                  Porites_evermani_scaffold_3893:3797-3890 Cluster_15890
    ## 64               Porites_evermani_scaffold_430:205865-205960  Cluster_5563
    ## 65               Porites_evermani_scaffold_461:215454-215549  Cluster_5882
    ## 66                Porites_evermani_scaffold_47:475972-476066  Cluster_1140
    ## 67                Porites_evermani_scaffold_49:151587-151681  Cluster_1167
    ## 68                Porites_evermani_scaffold_5010:12339-12433 Cluster_16498
    ## 69                 Porites_evermani_scaffold_502:58948-59038  Cluster_6255
    ## 70               Porites_evermani_scaffold_590:199106-199199  Cluster_6904
    ## 71               Porites_evermani_scaffold_590:199370-199465  Cluster_6905
    ## 72               Porites_evermani_scaffold_590:199803-199896  Cluster_6906
    ## 73               Porites_evermani_scaffold_594:158176-158270  Cluster_6914
    ## 74               Porites_evermani_scaffold_613:156453-156546  Cluster_7053
    ## 75                  Porites_evermani_scaffold_6219:6499-6590 Cluster_16738
    ## 76                 Porites_evermani_scaffold_730:81363-81456  Cluster_7657
    ## 77                 Porites_evermani_scaffold_730:82401-82494  Cluster_7658
    ## 78               Porites_evermani_scaffold_768:137973-138066  Cluster_7855
    ## 79                 Porites_evermani_scaffold_866:22803-22895  Cluster_8634
    ## 80               Porites_evermani_scaffold_910:118720-118809  Cluster_8887
    ## 81               Porites_evermani_scaffold_910:139331-139420  Cluster_8888
    ## 82                 Porites_evermani_scaffold_910:99233-99322  Cluster_8884
    ## 83               Porites_evermani_scaffold_942:133648-133739  Cluster_8988
    ## 84                 Porites_evermani_scaffold_984:51832-51924  Cluster_9149
    ## 85  Pocillopora_meandrina_HIv1___Sc0000000:20372416-20372510   Cluster_360
    ## 86    Pocillopora_meandrina_HIv1___Sc0000000:2872019-2872110    Cluster_36
    ## 87      Pocillopora_meandrina_HIv1___Sc0000000:818027-818120    Cluster_21
    ## 88    Pocillopora_meandrina_HIv1___Sc0000001:1459402-1459496   Cluster_390
    ## 89  Pocillopora_meandrina_HIv1___Sc0000001:19145766-19145858   Cluster_757
    ## 90  Pocillopora_meandrina_HIv1___Sc0000002:10845722-10845816  Cluster_1015
    ## 91  Pocillopora_meandrina_HIv1___Sc0000002:15749261-15749351  Cluster_1068
    ## 92  Pocillopora_meandrina_HIv1___Sc0000002:16106559-16106650  Cluster_1080
    ## 93    Pocillopora_meandrina_HIv1___Sc0000002:3841944-3842042   Cluster_925
    ## 94  Pocillopora_meandrina_HIv1___Sc0000003:10129512-10129606  Cluster_1289
    ## 95  Pocillopora_meandrina_HIv1___Sc0000003:10366033-10366130  Cluster_1296
    ## 96      Pocillopora_meandrina_HIv1___Sc0000003:495066-495158  Cluster_1116
    ## 97  Pocillopora_meandrina_HIv1___Sc0000005:10385497-10385597  Cluster_1938
    ## 98  Pocillopora_meandrina_HIv1___Sc0000005:11260734-11260823  Cluster_1952
    ## 99  Pocillopora_meandrina_HIv1___Sc0000005:11261356-11261449  Cluster_1953
    ## 100     Pocillopora_meandrina_HIv1___Sc0000005:601572-601666  Cluster_1793
    ## 101   Pocillopora_meandrina_HIv1___Sc0000008:1783802-1783891  Cluster_2793
    ## 102   Pocillopora_meandrina_HIv1___Sc0000008:3576420-3576514  Cluster_2837
    ## 103   Pocillopora_meandrina_HIv1___Sc0000008:3619341-3619435  Cluster_2839
    ## 104   Pocillopora_meandrina_HIv1___Sc0000008:5387738-5387830  Cluster_2859
    ## 105   Pocillopora_meandrina_HIv1___Sc0000009:3894900-3894992  Cluster_2973
    ## 106   Pocillopora_meandrina_HIv1___Sc0000010:8012916-8013008  Cluster_3392
    ## 107   Pocillopora_meandrina_HIv1___Sc0000012:6455183-6455276  Cluster_3661
    ## 108   Pocillopora_meandrina_HIv1___Sc0000014:2366038-2366132  Cluster_4039
    ## 109   Pocillopora_meandrina_HIv1___Sc0000014:2594743-2594833  Cluster_4040
    ## 110   Pocillopora_meandrina_HIv1___Sc0000014:9365459-9365552  Cluster_4118
    ## 111   Pocillopora_meandrina_HIv1___Sc0000016:7549557-7549654  Cluster_4435
    ## 112   Pocillopora_meandrina_HIv1___Sc0000016:7550593-7550685  Cluster_4437
    ## 113   Pocillopora_meandrina_HIv1___Sc0000016:7553282-7553379  Cluster_4438
    ## 114   Pocillopora_meandrina_HIv1___Sc0000016:7555639-7555733  Cluster_4439
    ## 115   Pocillopora_meandrina_HIv1___Sc0000017:5050908-5051000  Cluster_4572
    ## 116   Pocillopora_meandrina_HIv1___Sc0000018:4564884-4564976  Cluster_4754
    ## 117   Pocillopora_meandrina_HIv1___Sc0000018:6855499-6855592  Cluster_4823
    ## 118   Pocillopora_meandrina_HIv1___Sc0000021:4351817-4351909  Cluster_5253
    ## 119   Pocillopora_meandrina_HIv1___Sc0000024:4808666-4808760  Cluster_5612
    ## 120   Pocillopora_meandrina_HIv1___Sc0000026:1154719-1154813  Cluster_5740
    ## 121   Pocillopora_meandrina_HIv1___Sc0000035:1989820-1989912  Cluster_6382
    ##                                      Chrom    Start      End Length   Reads
    ## 1                                 ntLink_6 13351746 13351842     97   14107
    ## 2                                 ntLink_6  4847443  4847535     93   12648
    ## 3                                 ntLink_6  5157537  5157626     90   11203
    ## 4                                 ntLink_6  7263486  7263580     95    7470
    ## 5                               ptg000001l 20063072 20063168     97    4553
    ## 6                               ptg000001l  5548841  5548934     94   80076
    ## 7                               ptg000002l 14046235 14046328     94    4003
    ## 8                               ptg000002l 14046541 14046634     94    5645
    ## 9                               ptg000002l  7337508  7337601     94     346
    ## 10                              ptg000004l  1859862  1859953     92   10448
    ## 11                              ptg000007l  3377283  3377376     94     332
    ## 12                              ptg000007l   915905   916002     98     414
    ## 13                              ptg000008l 10754767 10754862     96   75971
    ## 14                              ptg000009l  4940515  4940609     95   11540
    ## 15                              ptg000009l   616285   616377     93    2056
    ## 16                              ptg000009l   616579   616673     95   80457
    ## 17                              ptg000016l 11751358 11751448     91   25321
    ## 18                              ptg000016l  7795508  7795602     95   24796
    ## 19                              ptg000016l  8599833  8599925     93   24703
    ## 20                              ptg000017l  6374963  6375054     92   29473
    ## 21                              ptg000017l  7471118  7471210     93    1843
    ## 22                              ptg000018l  2286780  2286870     91    2790
    ## 23                              ptg000023l 18708903 18708996     94     926
    ## 24                              ptg000023l 37965243 37965338     96      11
    ## 25                              ptg000024l  4086202  4086295     94   15328
    ## 26                              ptg000024l  5256454  5256549     96    7544
    ## 27                              ptg000025l 10501030 10501125     96     397
    ## 28                              ptg000025l 10668901 10668998     98     232
    ## 29                              ptg000025l 15156415 15156507     93    1415
    ## 30                              ptg000025l  1981692  1981783     92     837
    ## 31                              ptg000025l  7472531  7472623     93   45209
    ## 32                              ptg000026l  8745540  8745635     96    2854
    ## 33                              ptg000031l  5461305  5461397     93     149
    ## 34                              ptg000031l  6751935  6752029     95  249212
    ## 35                              ptg000035l  4339550  4339641     92   50361
    ## 36                              ptg000035l  4808338  4808432     95    1270
    ## 37                              ptg000035l  5346032  5346127     96 1360205
    ## 38                              ptg000035l  8367726  8367821     96    1830
    ## 39                              ptg000039l    35764    35856     93     579
    ## 40          Porites_evermani_scaffold_1128   125526   125618     93    2705
    ## 41          Porites_evermani_scaffold_1159     6653     6743     91   12683
    ## 42          Porites_evermani_scaffold_1159     7406     7500     95    6558
    ## 43           Porites_evermani_scaffold_138   127944   128038     95    3344
    ## 44          Porites_evermani_scaffold_1415    72049    72142     94     717
    ## 45          Porites_evermani_scaffold_1429    47235    47326     92   13224
    ## 46           Porites_evermani_scaffold_145     5384     5473     90     314
    ## 47           Porites_evermani_scaffold_148    94531    94624     94    7981
    ## 48          Porites_evermani_scaffold_1503    47558    47651     94   54746
    ## 49          Porites_evermani_scaffold_1503    47749    47842     94  208179
    ## 50            Porites_evermani_scaffold_16   383386   383478     93    4060
    ## 51          Porites_evermani_scaffold_1732    76482    76574     93    1108
    ## 52             Porites_evermani_scaffold_1  1404250  1404342     93    9574
    ## 53          Porites_evermani_scaffold_2259    40195    40285     91     458
    ## 54           Porites_evermani_scaffold_253   200943   201039     97   44919
    ## 55           Porites_evermani_scaffold_253   202235   202327     93   62077
    ## 56           Porites_evermani_scaffold_257   110308   110402     95   20997
    ## 57            Porites_evermani_scaffold_26   382550   382645     96   93348
    ## 58          Porites_evermani_scaffold_2738    56836    56926     91     190
    ## 59          Porites_evermani_scaffold_3072    29428    29521     94     609
    ## 60           Porites_evermani_scaffold_316    88415    88506     92    1390
    ## 61           Porites_evermani_scaffold_334   153554   153646     93    1061
    ## 62          Porites_evermani_scaffold_3707    36102    36199     98    1764
    ## 63          Porites_evermani_scaffold_3893     3797     3890     94    1787
    ## 64           Porites_evermani_scaffold_430   205865   205960     96    2539
    ## 65           Porites_evermani_scaffold_461   215454   215549     96    1353
    ## 66            Porites_evermani_scaffold_47   475972   476066     95    3061
    ## 67            Porites_evermani_scaffold_49   151587   151681     95  267836
    ## 68          Porites_evermani_scaffold_5010    12339    12433     95    2367
    ## 69           Porites_evermani_scaffold_502    58948    59038     91     554
    ## 70           Porites_evermani_scaffold_590   199106   199199     94    2976
    ## 71           Porites_evermani_scaffold_590   199370   199465     96    1333
    ## 72           Porites_evermani_scaffold_590   199803   199896     94    2870
    ## 73           Porites_evermani_scaffold_594   158176   158270     95   44619
    ## 74           Porites_evermani_scaffold_613   156453   156546     94    5980
    ## 75          Porites_evermani_scaffold_6219     6499     6590     92     785
    ## 76           Porites_evermani_scaffold_730    81363    81456     94    2935
    ## 77           Porites_evermani_scaffold_730    82401    82494     94    1351
    ## 78           Porites_evermani_scaffold_768   137973   138066     94     353
    ## 79           Porites_evermani_scaffold_866    22803    22895     93     405
    ## 80           Porites_evermani_scaffold_910   118720   118809     90   67563
    ## 81           Porites_evermani_scaffold_910   139331   139420     90   21132
    ## 82           Porites_evermani_scaffold_910    99233    99322     90   35797
    ## 83           Porites_evermani_scaffold_942   133648   133739     92   32177
    ## 84           Porites_evermani_scaffold_984    51832    51924     93     191
    ## 85  Pocillopora_meandrina_HIv1___Sc0000000 20372416 20372510     95  522228
    ## 86  Pocillopora_meandrina_HIv1___Sc0000000  2872019  2872110     92     177
    ## 87  Pocillopora_meandrina_HIv1___Sc0000000   818027   818120     94   12096
    ## 88  Pocillopora_meandrina_HIv1___Sc0000001  1459402  1459496     95    3332
    ## 89  Pocillopora_meandrina_HIv1___Sc0000001 19145766 19145858     93    6322
    ## 90  Pocillopora_meandrina_HIv1___Sc0000002 10845722 10845816     95   12003
    ## 91  Pocillopora_meandrina_HIv1___Sc0000002 15749261 15749351     91     267
    ## 92  Pocillopora_meandrina_HIv1___Sc0000002 16106559 16106650     92    2969
    ## 93  Pocillopora_meandrina_HIv1___Sc0000002  3841944  3842042     99  106356
    ## 94  Pocillopora_meandrina_HIv1___Sc0000003 10129512 10129606     95   11292
    ## 95  Pocillopora_meandrina_HIv1___Sc0000003 10366033 10366130     98  858857
    ## 96  Pocillopora_meandrina_HIv1___Sc0000003   495066   495158     93     990
    ## 97  Pocillopora_meandrina_HIv1___Sc0000005 10385497 10385597    101   58190
    ## 98  Pocillopora_meandrina_HIv1___Sc0000005 11260734 11260823     90   30986
    ## 99  Pocillopora_meandrina_HIv1___Sc0000005 11261356 11261449     94   11335
    ## 100 Pocillopora_meandrina_HIv1___Sc0000005   601572   601666     95  129423
    ## 101 Pocillopora_meandrina_HIv1___Sc0000008  1783802  1783891     90    1051
    ## 102 Pocillopora_meandrina_HIv1___Sc0000008  3576420  3576514     95     301
    ## 103 Pocillopora_meandrina_HIv1___Sc0000008  3619341  3619435     95     366
    ## 104 Pocillopora_meandrina_HIv1___Sc0000008  5387738  5387830     93    1517
    ## 105 Pocillopora_meandrina_HIv1___Sc0000009  3894900  3894992     93     578
    ## 106 Pocillopora_meandrina_HIv1___Sc0000010  8012916  8013008     93   26666
    ## 107 Pocillopora_meandrina_HIv1___Sc0000012  6455183  6455276     94    2417
    ## 108 Pocillopora_meandrina_HIv1___Sc0000014  2366038  2366132     95    1393
    ## 109 Pocillopora_meandrina_HIv1___Sc0000014  2594743  2594833     91     329
    ## 110 Pocillopora_meandrina_HIv1___Sc0000014  9365459  9365552     94   11095
    ## 111 Pocillopora_meandrina_HIv1___Sc0000016  7549557  7549654     98   19764
    ## 112 Pocillopora_meandrina_HIv1___Sc0000016  7550593  7550685     93   88654
    ## 113 Pocillopora_meandrina_HIv1___Sc0000016  7553282  7553379     98    4991
    ## 114 Pocillopora_meandrina_HIv1___Sc0000016  7555639  7555733     95   15915
    ## 115 Pocillopora_meandrina_HIv1___Sc0000017  5050908  5051000     93  132750
    ## 116 Pocillopora_meandrina_HIv1___Sc0000018  4564884  4564976     93     530
    ## 117 Pocillopora_meandrina_HIv1___Sc0000018  6855499  6855592     94    6913
    ## 118 Pocillopora_meandrina_HIv1___Sc0000021  4351817  4351909     93     163
    ## 119 Pocillopora_meandrina_HIv1___Sc0000024  4808666  4808760     95    3154
    ## 120 Pocillopora_meandrina_HIv1___Sc0000026  1154719  1154813     95   23971
    ## 121 Pocillopora_meandrina_HIv1___Sc0000035  1989820  1989912     93    4602
    ##     DistinctSequences FracTop Strand                  MajorRNA MajorRNAReads
    ## 1                 186   0.000      -    UAAAAUGUCGGUUGCUUAAGCU          7452
    ## 2                 343   0.001      -    AUGAUCAUAGCACUUUCUUUGU          3018
    ## 3                 216   0.998      +     CAAUGUUUCGGCUUGUUCCCG          5366
    ## 4                 224   0.000      -  UUUCAAAUUAGGAAGGGAGGUGUU          3004
    ## 5                 186   1.000      +  AAGAACACCCAAAAUAGCUGAGGA           982
    ## 6                 322   0.003      -    UCUCAGAUUACAGUAGUUAAGU         56403
    ## 7                 184   1.000      + UCUGGCAGUAUGUUAUUUUUCCAAU           802
    ## 8                 220   1.000      +   UCUGGCAGUAUGUUAUUUUUCCC          1012
    ## 9                  45   1.000      +    UCAUAACAGUGAGGACCAUUCU           127
    ## 10                227   0.000      -   UAGCAUAACAUUGUAAGAGAUCU          5104
    ## 11                 37   1.000      +    UUUGAUUGCUGUGAUCUGGUUA           143
    ## 12                 39   0.000      -    AAAGGAAGAAAGGAGACCCUGG           139
    ## 13                546   0.000      -     AAAGAAGUACAAGUGGUAGGG         18322
    ## 14                210   0.000      -   CAAGUGAGAGAAGGUUAGUGUGG          5129
    ## 15                 84   0.000      -    AUUUCACUAGAUGAGCGCUAAC          1178
    ## 16                407   0.000      -    UCUGCGUUAUCGGUGAAAUUGU         52366
    ## 17                287   0.000      -    AGUGCACUUUUCUCAGGAUGAA          9236
    ## 18                199   1.000      +    UGGGUGUCAUCUAUUAUGUUUU         18199
    ## 19                398   0.000      -    UAAUGUUCGCAACUGCCUUGUU          5902
    ## 20                383   0.000      -     UAGGCGUAUUUCCGAUUGUCC          9407
    ## 21                 89   1.000      +   UUUGCUAGUUGCUUUUGUCCCGU           912
    ## 22                143   1.000      +    UUCUUUGCACUGUUCAAACCCU           895
    ## 23                 90   1.000      +    UUCGUCUCUUGUAUUUUCCUGG           256
    ## 24                  5   1.000      +     UCGGACACCUGUAAUUGGAUA             5
    ## 25                136   1.000      +   AUUUUUAGCCCGCGGAAGUUGCU          7236
    ## 26                271   0.000      -   UACAUGUGCUGCACUGGAACUCU          2505
    ## 27                 53   0.997      +    AUUGAUGAGCACAACAAAAGGU           143
    ## 28                 25   0.000      -   UAAAGCUUUUGUGAAGAAACACG            94
    ## 29                 87   1.000      +    UUACCACAUUCAGGAUAGACAU           616
    ## 30                 88   0.999      +   UAUUGAUUAGUACUUGAUAGACA           315
    ## 31                444   0.000      -   UACAAAAACAAGAUGAGUGCAGG         21945
    ## 32                 64   0.000      -    UGUGAUUGGAGACUUUUAUCGU          1690
    ## 33                 28   0.000      -    ACGCUAGGAAGGGAUGCCGGGA            84
    ## 34                491   0.000      -  UUAACGAGUAGAUAAAUGAAGAGU        149551
    ## 35                407   0.000      -   AUUGAUUGUAGACAAGCCUCUGA         12975
    ## 36                101   1.000      +    CUCACUGCUUAAAAGAGUCACU           502
    ## 37                589   1.000      +      UCCCGUAGAUCCGAACUUGU        571297
    ## 38                122   0.000      -  UAUAUUGUACGACUCUCAUCGUGU           533
    ## 39                 64   0.000      -    CAAGGAGGAAGCAUGAUACGUA           130
    ## 40                 97   0.001      -     CAUGUGUUAAUUUGGACGAUU          1459
    ## 41                256   0.992      +   UGAGGUAAUUGUUGUUUAAAACU          2958
    ## 42                252   1.000      +   UAAGGUUGUUGUUGCUUAAAACU          1471
    ## 43                133   1.000      +   UCACUGAAGCCGUCUAUCUUGAU           924
    ## 44                 73   1.000      +  ACAGGUGAAUGUUGUGUGAUCGCU           218
    ## 45                205   0.000      -    UGCAUCACUAAAAUUUCUCAAU          6555
    ## 46                 24   0.984      +    UUGGGAAAUCCGGAUUUAGAUU           107
    ## 47                155   0.000      -     UAAUUGUCCUUUAGUGUUUGU          5087
    ## 48                235   1.000      +    UCAAGUCUAGGCUGGUUAGUUU         43003
    ## 49                408   1.000      +    UCAGGUCUAGGCUGGUUAGUUU        151922
    ## 50                194   0.000      -    UUAAGGGCCUGGGAGAGAGAAG          1239
    ## 51                120   0.997      +    UUAUGCAGAAAGUAAUAGUUCA           171
    ## 52                313   0.002      -    UUAUGCUCUGUAUACUGGUAGU          2849
    ## 53                 30   0.856      +    UCUAAAUCCGGAUUUCCCAAUC           324
    ## 54                606   1.000      +  UAGAAAACUCGUGUACGUGACCCU          6501
    ## 55                554   1.000      +      UAAAAAACUCGUCGACGUGA         17855
    ## 56                517   0.000      -     AAGCACAACAUUAAUAGGAGU          8549
    ## 57                537   0.000      -    UCUCAGCUCACCAAUCUCUGCU         55070
    ## 58                 33   1.000      +        UUUCUCUAAUCUCUAUGG            56
    ## 59                 40   1.000      +    AUUUUUAGCCCGCGGAAGUUGC           325
    ## 60                 52   0.000      -   UCCCCGCCAGGAAUUCUUGUAUU           458
    ## 61                119   0.051      -      UGCAGGUACAGUUAUAAGGU            75
    ## 62                 89   0.000      -       UUGAACACUGAAACCCUGU           583
    ## 63                140   0.108      -    UUAAGAAUAUCAGUGAUAGGGA           429
    ## 64                137   0.000      -      UAUAUUGUACGACUCUCAUC           600
    ## 65                 74   1.000      +   CGGUUCUUGUAUGCUGUUCUAAU           605
    ## 66                140   0.000      -   UAAUUAUAUCAGAAGGAAGAAGA          1119
    ## 67                483   0.000      -      UCCCGUAGAUCCGAACUUGU        173986
    ## 68                136   0.995      +    UUAUAUUAUAGACUUUCCUUGU           733
    ## 69                 39   0.000      -         UAGCAUAACAUUGUAAG           345
    ## 70                140   1.000      +   UCACGAAUGGCUAGCGCACGAUU           968
    ## 71                 70   1.000      +   UCACGAAUGGCUAGCGCACGAUU           450
    ## 72                135   1.000      +   UCACGAAUGGCUAGUGCACGAUU          1102
    ## 73                365   1.000      +     AAAGAAGUACAAGUGGUAGGG         12472
    ## 74                 96   1.000      +    UUCCUAAGUUGUUGUAGCUGUU          3276
    ## 75                 94   0.181      -    AUGUGUACAGCUAAGAGCAGUU           232
    ## 76                185   0.000      -    UGUUAGUUUACAGUUAGUGUUU           649
    ## 77                115   0.001      -    UGUUAGUUUACAGUUAGUGUUU           299
    ## 78                 48   1.000      +    UUUUUCGUGUGCUUAUGCUUAU           115
    ## 79                 47   0.000      -    UAUUUCAAACGUCGCGUAACCU           199
    ## 80                167   0.944      +     CAAUGUUUCGGCUUGUUCCCG         59917
    ## 81                222   0.977      +      AAAUGUUUCGGCUUGUUCCC          9474
    ## 82                362   0.976      +     CAAUGUUUCGGCUUGUUCCCA         13239
    ## 83                505   1.000      +   UAGGUAACUCUGUUAGUUUGUUU          8153
    ## 84                 24   0.000      -      UCUGAGCCUGUAGCCCAUUU            45
    ## 85                378   0.000      -    UCAAGUCUAGGCUGGUUAGUUU        361684
    ## 86                 29   1.000      +     AUAGUAUUGUCGACAGGUACU            63
    ## 87                283   0.999      +   UGUUCUCUCUGCAGUAAGCAUGU          3258
    ## 88                118   1.000      +    AAAGAAGUGCAAGUGGUAGAGA          1704
    ## 89                255   1.000      +   CGCUCUUAAAAUCUCCUAGCUAU          1605
    ## 90                340   0.009      -   UUGCCGUUAACCCAAAACCUUCA          4844
    ## 91                 31   1.000      +    ACGCUAGGAAGGGAUGCCGGGA           162
    ## 92                 94   1.000      +   UGUUAUUUUUAUCGUCCAUUAGU          1172
    ## 93                263   0.000      -    UUGUGUAACUCCCUAAGGAAGG         66971
    ## 94                233   0.000      -   CGAUGUACUUCAGAGCAAUUUUU          2869
    ## 95                595   1.000      +      ACCCGUAGAUCCGAACUUGU        302900
    ## 96                 64   1.000      +    UAUAUUGUACGACUCUCAUCGU           438
    ## 97                312   1.000      +   UCAGGGAUUGUGGUGAGUUAGUU         21926
    ## 98                390   1.000      +   UCAGGUUAGGAGCAUGGUAUAGA         11353
    ## 99                284   1.000      +    UAUGCCAGGUCCUAACUUGAAA          3074
    ## 100               404   1.000      +     AAAGAAGUACAAGUGGUAGGG         43823
    ## 101                80   1.000      +      UAUGUUUUGGCUUGUUCCCA           330
    ## 102                45   1.000      +    CGAUGGAACACACAAGUAGAGA            89
    ## 103                65   0.005      -    UGUGCUUUGGUUCAAAAAUGUC           148
    ## 104               131   0.002      -   AACUUUUCUACAAUUCCAGUGGA           274
    ## 105                48   0.002      -         UAGCAUAACAUUGUAAG           330
    ## 106               313   0.991      +    UACUUCAAUGUUAUCGGUUUGU          7987
    ## 107               125   1.000      +   UGAAAUACUCUGACGGAGUCAGC          1099
    ## 108                67   0.000      -    AUUUUUAGCCCGCGGAAGUUGC           641
    ## 109                44   0.997      +    CGUGGUGCGUUAUUGUUUCGGU            54
    ## 110               130   0.000      -   UUCCAGAUUUGUAUCAUUGUUAU          7555
    ## 111               352   1.000      +    UCAGUUCCACCAGCUCACCAAU          5579
    ## 112               512   1.000      +    UCAGUCCCACCAUCUCACCAAU         40428
    ## 113               191   1.000      +    AGGUGAGCUGUGGGAACUCUUA          2482
    ## 114               280   1.000      +    UCAGUCCCACCAGCUCACCAAU          3822
    ## 115               610   0.000      -   AACUGCUGAGAUUCUAUGGAUUU         47680
    ## 116                60   1.000      +   UAACCUUAUUUGUCUUCACAUGU           143
    ## 117               289   1.000      +   GGAUUAAAAUUUGUUGGGUGUGC          1244
    ## 118                39   1.000      +        ACUGAUAUUCACCAAGUG            34
    ## 119               144   1.000      +        AGAACCCAAGAAUCUCGA           894
    ## 120               150   0.000      -   UAAGAUCAAGUUCAUGGUAUUCU         14088
    ## 121               200   1.000      +    UAUUUACAACUCUCAAAACAAC          1726
    ##      Short Long   X21    X22    X23    X24 DicerCall MIRNA
    ## 1      653    7   468  12598    318     63        22     Y
    ## 2     1048    9   758   5669   4675    489        22     Y
    ## 3     3224   34  7253    334    344     14         N     Y
    ## 4      124 2136   131    115    826   4138         N     Y
    ## 5      329  144  1035    476   1435   1134        23     Y
    ## 6     1600   61  2455  70964   4578    418        22     Y
    ## 7      229 1294   834    736    795    115         N     Y
    ## 8      291 1076  1032    939   1507    800         N     Y
    ## 9       37    0    20    231     49      9        22     Y
    ## 10     475  296    99   1109   7865    604        23     Y
    ## 11      22    0    80    167     59      4        22     Y
    ## 12      71    0    20    293     30      0        22     Y
    ## 13   24349   10 33761  15714   1522    615         N     Y
    ## 14     209   73  2546    285   7999    428        23     Y
    ## 15     187    5    33   1452    377      2        22     Y
    ## 16    2468 1367  2466  71781   2281     94        22     Y
    ## 17    4881    7  5372  10326   4494    241        22     Y
    ## 18     462    7  1053  21591   1646     37        22     Y
    ## 19     607  211  1680  11437   8547   2221        22     Y
    ## 20    3824  189 22149   2297    717    297        21     Y
    ## 21      24   28    33     58   1485    215        23     Y
    ## 22     147   63   251   1457    279    593        22     Y
    ## 23      69    2   105    549    185     16        22     Y
    ## 24       3    0     7      0      1      0         N     Y
    ## 25     176   86   447   3022  11118    479        23     Y
    ## 26     289  381  1410    796   4070    598        23     Y
    ## 27      57    0    54    243     42      1        22     Y
    ## 28      10    1     5      7    177     32        23     Y
    ## 29     185   36   436    668     89      1        22     Y
    ## 30      98   10    31    189    482     27        23     Y
    ## 31    5021  267  4534   2938  31330   1119        23     Y
    ## 32     133    3    43   2337    336      2        22     Y
    ## 33       7    0    13    110     18      1        22     Y
    ## 34    3846 1813 14249  11217  37270 180817        24     Y
    ## 35    3270 3059  4025   3782  27285   8940        23     Y
    ## 36     134    7    46    733    321     29        22     Y
    ## 37  758793  417 31482 556861  12194    458         N     Y
    ## 38     215  104    25    160    571    755        24     Y
    ## 39      55    0   179    205    136      4        22     Y
    ## 40     798    0  1727    157     20      3         N     Y
    ## 41    2218   44  1357   1788   7247     29        23     Y
    ## 42     492   82  1349    346   4053    236        23     Y
    ## 43     878  108   125    887   1116    230         N     Y
    ## 44     144   24    70     48    172    259         N     Y
    ## 45    4484    5   227   7951    495     62         N     Y
    ## 46      61    9    71    155      9      9         N     Y
    ## 47     481    8  6125    241    968    158        21     Y
    ## 48    1714   37  2952  48770   1092    181        22     Y
    ## 49    8215  251 19527 172650   7006    530        22     Y
    ## 50    1297   46   158   1819    699     41         N     Y
    ## 51     262  101   279    335     96     35         N     Y
    ## 52    3870   80   468   4427    685     44         N     Y
    ## 53      53    7    61    333      4      0        22     Y
    ## 54   15384 2459  1381   8809   6453  10433         N     Y
    ## 55   48007  477  3649   7018   1002   1924         N     Y
    ## 56    2899  133 11315   2395   3478    777        21     Y
    ## 57    3610 1891  2578  71538  12756    975        22     Y
    ## 58      87    1    27     51     11     13         N     Y
    ## 59      22    0    29    514     44      0        22     Y
    ## 60     155   12   102    531    578     12        23     Y
    ## 61     500    3   196    266     90      6         N     Y
    ## 62     974    9   102    590     38     51         N     Y
    ## 63     375    0   404    574    358     76         N     Y
    ## 64    1259   46    46    107    480    601         N     Y
    ## 65     397   44    24    174    681     33         N     Y
    ## 66     334  114   218    268   1384    743        23     Y
    ## 67  216104   92 14714  35566   1169    191         N     Y
    ## 68     599    1   589    995    138     45         N     Y
    ## 69     391    2    16     42    103      0         N     Y
    ## 70     715   18   131    963   1133     16         N     Y
    ## 71     305    0    68    272    629     59         N     Y
    ## 72     491    1   215    823   1309     31        23     Y
    ## 73   24239    0 16335   3821    216      8         N     Y
    ## 74     662   19   447   4484    368      0        22     Y
    ## 75     225    0   132    373     40     15         N     Y
    ## 76     414  160   691    756    577    337        22     Y
    ## 77     180   32   353    352    303    131        21     Y
    ## 78      80    0    14    214     37      8         N     Y
    ## 79      79    6    27    251     37      5         N     Y
    ## 80    3615   66 62951    413    503     15        21     Y
    ## 81   19234    4  1717    108     57     12         N     Y
    ## 82   18134  120 16648    693    173     29         N     Y
    ## 83    7669  302  3095   5335  10466   5310         N     Y
    ## 84      59   28     7     56     35      6         N     Y
    ## 85    1294   19 20216 498308   2255    136        22     Y
    ## 86      16    0   116     42      3      0        21     Y
    ## 87     393   78   202   3329   5487   2607        23     Y
    ## 88     629    3   340   2293     56     11        22     Y
    ## 89    1316   82   585   1212   2501    626         N     Y
    ## 90     286  129   442   5105   5860    181        23     Y
    ## 91       9    4    32    201     19      2        22     Y
    ## 92    1055    9    61    209   1607     28         N     Y
    ## 93    4181    3  2786  98832    535     19        22     Y
    ## 94     397  119  3756    240   5015   1765        23     Y
    ## 95  453621   35 21036 380606   3365    194         N     Y
    ## 96      71    6    12    681    173     47        22     Y
    ## 97   10279  126  1422  10174  32234   3955        23     Y
    ## 98   11069  138  1185   4477  13443    674         N     Y
    ## 99    2510   37  1887   4105   2298    498         N     Y
    ## 100  39962    4 68302  20915    235      5         N     Y
    ## 101    914   11   105     14      5      2         N     Y
    ## 102     31    9    14    185     60      2        22     Y
    ## 103     65    6    15    212     56     12        22     Y
    ## 104    246   37   133    391    660     50        23     Y
    ## 105    351    7     3    170     45      2         N     Y
    ## 106    397  129   142  12953  11855   1190        22     Y
    ## 107    309  151    42    196   1432    287        23     Y
    ## 108     54    1    97    930    301     10        22     Y
    ## 109     60    0    17    167     49     36        22     Y
    ## 110     84   19   144    240  10206    402        23     Y
    ## 111   2156   16  1034   9200   7318     40        22     Y
    ## 112  12735  185  7610  60346   7494    284        22     Y
    ## 113    362   35   782   3195    529     88        22     Y
    ## 114   3453   42  3693   5606   3059     62         N     Y
    ## 115   4764  285  1096  22198 101160   3247        23     Y
    ## 116     35    0    42    119    278     56        23     Y
    ## 117    738  316   522    428   3428   1481        23     Y
    ## 118     72    4    45     26     13      3         N     Y
    ## 119   2822    1   301     26      3      1         N     Y
    ## 120     61   53    21    570  22819    447        23     Y
    ## 121    614   49   387   2009   1506     37        22     Y
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            known_miRNAs
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  <NA>
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  <NA>
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ami-spi-miR-L-temp-14-5p_Acropora_millepora_Praher_et_al._2021_NA
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Adi-Mir-Novel-2_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2026-5p
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 nve-miR-9425;Adi-Mir-9425_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-9425;_nve-miR-9425;ami-nve-F-miR-9425-5p_Acropora_millepora_Praher_et_al._2021_NA;miR-9425_Nematostella_vectensis_Moran_et_al._2014_NA
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  <NA>
    ## 7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       adi-miR-P-novel-2-3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       adi-miR-P-novel-3-3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       adi-miR-P-novel-7-3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 10                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       Adi-Mir-2030_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2030-5p;_nve-miR-2030-5p;_spi-miR-temp-40;ami-nve-F-miR-2030-5p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-2030_Acropora_digitifera__Praher_et_al._2021_NA;spi-mir-temp-40_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2030;eca-nve-F-miR-2030_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;miR-2030_Nematostella_vectensis_Moran_et_al._2014_NA
    ## 11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              apa-mir-2050_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-2050;_Nve;_Spis;_Adi;Adi-Mir-2050_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2050-3p;_nve-miR-2050-3p;ami-nve-F-miR-2050-3p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-2050_Acropora_digitifera__Praher_et_al._2021_NA
    ## 12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 13                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  nve-miR-2023-3p;spi-mir-temp-4_Stylophora_pistillata_Liew_et_al._2014_Exact_match_of_nve-miR-2023.;apa-mir-2023_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-2023;_Nve;_Spis;_Adi;mse-nve-F-miR-2023_Metridium_senile_Praher_et_al._2021_Transcriptome-level;sca-nve-miR-2023-3p_Scolanthus_callimorphus_Praher_et_al._2021_Transcriptome-level;Adi-Mir-2023_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2023-3p;_nve-miR-2023-3p;_spi-miR-temp-4;ami-nve-F-miR-2023-3p_Acropora_millepora_Praher_et_al._2021_NA;epa-nve-F-miR-2023_Exaiptasia_pallida_Praher_et_al._2021_NA;spi-nve-F-miR-2023_Stylophora_pistillata_Praher_et_al._2021_NA;miR-2023_Nematostella_vectensis_Moran_et_al._2014_NA;avi-miR-temp-2023_Anemonia_viridis_Urbarova_et_al._2018_NA;eca-nve-F-miR-2023_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level
    ## 14                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Adi-Mir-Novel-5_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_NA;ami-Adi-MiR-G-Novel-5-3p_Acropora_millepora_Praher_et_al._2021_NA;Adi-MiR-G-Novel-5_3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 15                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 17                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 18                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 19                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Adi-Mir-Novel-3_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2028-5p
    ## 20                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 21                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         spi-mir-temp-25_Stylophora_pistillata_Liew_et_al._2014_Considered_bona_fideas_close_match_to_nve-_and_hma-miR-2022.;Adi-Mir-2022_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2022-3p;_nve-miR-2022-3p;_spi-miR-temp-25;ami-nve-F-miR-2022-3p_Acropora_millepora_Praher_et_al._2021_NA;eca-nve-F-miR-2022_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level
    ## 22                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 23                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 24                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ami-miR-P-novel-10-3p_Acropora_millepora_Praher_et_al._2021_NA;ami-miR-P-novel-11-3p_Acropora_millepora_Praher_et_al._2021_NA
    ## 25                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Adi-Mir-2025_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2025-3p;_nve-miR-2025-3p;ami-nve-F-miR-2025-3p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-2025_Acropora_digitifera__Praher_et_al._2021_NA
    ## 26                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 27                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      adi-miR-P-novel-4-3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 28                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ami-miR-P-novel-17-3p_Acropora_millepora_Praher_et_al._2021_NA
    ## 29                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 30                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 31                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ami-miR-P-novel-3-5p_Acropora_millepora_Praher_et_al._2021_NA
    ## 32                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       spi-mir-temp-20_Stylophora_pistillata_Liew_et_al._2014_NA;apa-mir-2037_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-2037;_Nve;_Spis;mse-nve-F-miR-2037-3p_Metridium_senile_Praher_et_al._2021_Transcriptome-level;sca-nve-F-miR-2037-3p_Scolanthus_callimorphus_Praher_et_al._2021_Transcriptome-level;ami-nve-F-miR-2037-3p_Acropora_millepora_Praher_et_al._2021_NA;spi-nve-F-miR-2037_Stylophora_pistillata_Praher_et_al._2021_NA;avi-miR-temp-2037_Anemonia_viridis_Urbarova_et_al._2018_NA;eca-nve-F-miR-2037_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level
    ## 33                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 34                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Adi-Mir-Novel-1_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-9437;ami-miR-P-novel-1-3p_Acropora_millepora_Praher_et_al._2021_NA;adi-miR-P-novel-1-3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 35                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 36                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 37                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Adi-Mir-100_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_mir-100;_spi-miR-temp-1;ami-nve-F-miR-100-5p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-100_Acropora_digitifera__Praher_et_al._2021_NA
    ## 38                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          spi-mir-temp-30_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2036.;spi-nve-F-miR-2036_Stylophora_pistillata_Praher_et_al._2021_NA;Adi-Mir-2036_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2036-3p;_nve-miR-2036-3p;_spi-miR-temp-30;ami-nve-F-miR-2036-3p_Acropora_millepora_Praher_et_al._2021_NA;eca-nve-F-miR-2036_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level
    ## 39                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ami-miR-P-novel-13-3p_Acropora_millepora_Praher_et_al._2021_NA
    ## 40                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 41                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 42                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 43                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 44                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 45                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 46                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 47                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 48                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              spi-mir-temp-15_Stylophora_pistillata_Liew_et_al._2014_NA;spi-miR-L-temp-15_Stylophora_pistillata_Praher_et_al._2021_NA
    ## 49                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 50                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 51                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 52                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 53                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 54                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 55                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 56                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 57                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 58                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 59                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              ami-nve-F-miR-2025-3p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-2025_Acropora_digitifera__Praher_et_al._2021_NA;Adi-Mir-2025_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2025-3p;_nve-miR-2025-3p
    ## 60                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 61                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 62                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 63                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 64                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          spi-mir-temp-30_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2036.;spi-nve-F-miR-2036_Stylophora_pistillata_Praher_et_al._2021_NA;ami-nve-F-miR-2036-3p_Acropora_millepora_Praher_et_al._2021_NA;Adi-Mir-2036_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2036-3p;_nve-miR-2036-3p;_spi-miR-temp-30;eca-nve-F-miR-2036_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level
    ## 65                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 66                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 67                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          ami-nve-F-miR-100-5p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-100_Acropora_digitifera__Praher_et_al._2021_NA;Adi-Mir-100_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_mir-100;_spi-miR-temp-1
    ## 68                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 69                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       ami-nve-F-miR-2030-5p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-2030_Acropora_digitifera__Praher_et_al._2021_NA;Adi-Mir-2030_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2030-5p;_nve-miR-2030-5p;_spi-miR-temp-40;spi-mir-temp-40_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2030;eca-nve-F-miR-2030_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;miR-2030_Nematostella_vectensis_Moran_et_al._2014_NA
    ## 70                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 71                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 72                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 73                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  nve-miR-2023-3p;spi-mir-temp-4_Stylophora_pistillata_Liew_et_al._2014_Exact_match_of_nve-miR-2023.;apa-mir-2023_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-2023;_Nve;_Spis;_Adi;eca-nve-F-miR-2023_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;mse-nve-F-miR-2023_Metridium_senile_Praher_et_al._2021_Transcriptome-level;sca-nve-miR-2023-3p_Scolanthus_callimorphus_Praher_et_al._2021_Transcriptome-level;ami-nve-F-miR-2023-3p_Acropora_millepora_Praher_et_al._2021_NA;epa-nve-F-miR-2023_Exaiptasia_pallida_Praher_et_al._2021_NA;spi-nve-F-miR-2023_Stylophora_pistillata_Praher_et_al._2021_NA;miR-2023_Nematostella_vectensis_Moran_et_al._2014_NA;avi-miR-temp-2023_Anemonia_viridis_Urbarova_et_al._2018_NA;Adi-Mir-2023_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2023-3p;_nve-miR-2023-3p;_spi-miR-temp-4
    ## 74                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 75                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 76                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 77                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 78                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 79                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 80                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          adi-miR-P-novel-6-3p_Acropora_digitifera__Praher_et_al._2021_NA;spi-mir-temp-14_Stylophora_pistillata_Liew_et_al._2014_NA;ami-spi-miR-L-temp-14-5p_Acropora_millepora_Praher_et_al._2021_NA
    ## 81                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ami-spi-miR-L-temp-14-5p_Acropora_millepora_Praher_et_al._2021_NA
    ## 82                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      adi-miR-P-novel-6-3p_Acropora_digitifera__Praher_et_al._2021_NA
    ## 83                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 84                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 85                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              spi-mir-temp-15_Stylophora_pistillata_Liew_et_al._2014_NA;spi-miR-L-temp-15_Stylophora_pistillata_Praher_et_al._2021_NA
    ## 86                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 87                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 88                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 89                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 90                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 91                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 92                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 93                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 94                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            spi-mir-temp-23_Stylophora_pistillata_Liew_et_al._2014_NA
    ## 95  hsa-miR-100-5p;mmu-miR-100-5p;rno-miR-100-5p;gga-miR-100-5p;aga-miR-100;dre-miR-100-5p;lla-miR-100;ggo-miR-100;age-miR-100;ppa-miR-100;ppy-miR-100;ptr-miR-100;mml-miR-100-5p;sla-miR-100;tni-miR-100;fru-miR-100;xtr-miR-100;mdo-miR-100-5p;ame-miR-100-5p;oan-miR-100-5p;tca-miR-100-5p;bta-miR-100;lgi-miR-100;sko-miR-100;bmo-miR-100;eca-miR-100;ssc-miR-100;bma-miR-100b;aae-miR-100;cqu-miR-100-5p;tgu-miR-100-5p;nvi-miR-100;pma-miR-100a-5p;aca-miR-100;sha-miR-100;ola-miR-100;cgr-miR-100-5p;mse-miR-100;ccr-miR-100;ipu-miR-100;pmi-miR-100-5p;bbe-miR-100-5p;ssa-miR-100a-5p;cpi-miR-100-5p;ami-miR-100-5p;cli-miR-100-5p;pbv-miR-100-5p;chi-miR-100-5p;tch-miR-100-5p;pal-miR-100-5p;tcf-miR-100;abu-miR-100;mze-miR-100;nbr-miR-100;oni-miR-100;pny-miR-100;gmo-miR-100a-5p;pte-miR-100b-5p;xla-miR-100-5p;cpo-miR-100-5p;dno-miR-100-5p;ocu-miR-100-5p;mmr-miR-100;dma-miR-100;cja-miR-100;nle-miR-100;sbo-miR-100;pha-miR-100;oga-miR-100;dqu-miR-100-5p;pca-miR-100-5p;spi-mir-temp-1_Stylophora_pistillata_Liew_et_al._2014_Matches_miR-100_family.;eca-nve-F-miR-100_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;mse-nve-F-miR-100_Metridium_senile_Praher_et_al._2021_Transcriptome-level;sca-nve-F-miR-100_Scolanthus_callimorphus_Praher_et_al._2021_Transcriptome-level;apa-mir-100_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-100;_Nve;_Spis;_Adi;epa-nve-F-miR-100_Exaiptasia_pallida_Praher_et_al._2021_NA;spi-nve-F-miR-100_Stylophora_pistillata_Praher_et_al._2021_NA;avi-miR-temp-100_Anemonia_viridis_Urbarova_et_al._2018_NA;miR-100_Nematostella_vectensis_Moran_et_al._2014_NA
    ## 96                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          spi-mir-temp-30_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2036.;eca-nve-F-miR-2036_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;ami-nve-F-miR-2036-3p_Acropora_millepora_Praher_et_al._2021_NA;spi-nve-F-miR-2036_Stylophora_pistillata_Praher_et_al._2021_NA;Adi-Mir-2036_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2036-3p;_nve-miR-2036-3p;_spi-miR-temp-30
    ## 97                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             spi-mir-temp-3_Stylophora_pistillata_Liew_et_al._2014_NA
    ## 98                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 99                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 <NA>
    ## 100                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 nve-miR-2023-3p;spi-mir-temp-4_Stylophora_pistillata_Liew_et_al._2014_Exact_match_of_nve-miR-2023.;eca-nve-F-miR-2023_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;mse-nve-F-miR-2023_Metridium_senile_Praher_et_al._2021_Transcriptome-level;sca-nve-miR-2023-3p_Scolanthus_callimorphus_Praher_et_al._2021_Transcriptome-level;ami-nve-F-miR-2023-3p_Acropora_millepora_Praher_et_al._2021_NA;apa-mir-2023_Exaiptasia_pallida_Baumgarten_et_al._2017_miR-2023;_Nve;_Spis;_Adi;epa-nve-F-miR-2023_Exaiptasia_pallida_Praher_et_al._2021_NA;spi-nve-F-miR-2023_Stylophora_pistillata_Praher_et_al._2021_NA;miR-2023_Nematostella_vectensis_Moran_et_al._2014_NA;Adi-Mir-2023_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2023-3p;_nve-miR-2023-3p;_spi-miR-temp-4;avi-miR-temp-2023_Anemonia_viridis_Urbarova_et_al._2018_NA
    ## 101                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 102                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 103                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 104                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 105                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      adi-nve-F-miR-2030_Acropora_digitifera__Praher_et_al._2021_NA;ami-nve-F-miR-2030-5p_Acropora_millepora_Praher_et_al._2021_NA;Adi-Mir-2030_5p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2030-5p;_nve-miR-2030-5p;_spi-miR-temp-40;eca-nve-F-miR-2030_Edwardsiella_carnea_Praher_et_al._2021_Transcriptome-level;spi-mir-temp-40_Stylophora_pistillata_Liew_et_al._2014_Close_match_of_nve-miR-2030;miR-2030_Nematostella_vectensis_Moran_et_al._2014_NA
    ## 106                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 107                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 108                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             ami-nve-F-miR-2025-3p_Acropora_millepora_Praher_et_al._2021_NA;adi-nve-F-miR-2025_Acropora_digitifera__Praher_et_al._2021_NA;Adi-Mir-2025_3p_Acropora_digitifera_Gajigan_&_Conaco_2017_nve-miR-2025-3p;_nve-miR-2025-3p
    ## 109                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 110                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 111                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 112                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 113                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 114                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 115                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 116                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 117                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 118                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           spi-mir-temp-34_Stylophora_pistillata_Liew_et_al._2014_NA
    ## 119                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 120                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ## 121                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                <NA>
    ##                           given_miRNA_name
    ## 1                         apul-mir-novel-3
    ## 2                         apul-mir-novel-6
    ## 3                        apul-mir-novel-27
    ## 4                         apul-mir-novel-5
    ## 5                            apul-mir-9425
    ## 6                        apul-mir-novel-13
    ## 7                        apul-mir-novel-28
    ## 8                        apul-mir-novel-29
    ## 9                        apul-mir-novel-21
    ## 10                           apul-mir-2030
    ## 11                           apul-mir-2050
    ## 12                       apul-mir-novel-10
    ## 13                           apul-mir-2023
    ## 14                        apul-mir-novel-9
    ## 15                       apul-mir-novel-30
    ## 16    apul-mir-novel-8a; apul-mir-novel-8b
    ## 17                       apul-mir-novel-20
    ## 18                       apul-mir-novel-15
    ## 19                           apul-mir-2028
    ## 20                        apul-mir-novel-4
    ## 21                           apul-mir-2022
    ## 22                       apul-mir-novel-14
    ## 23                       apul-mir-novel-22
    ## 24                       apul-mir-novel-26
    ## 25                           apul-mir-2025
    ## 26                       apul-mir-novel-19
    ## 27                       apul-mir-novel-11
    ## 28                       apul-mir-novel-12
    ## 29                       apul-mir-novel-31
    ## 30                       apul-mir-novel-32
    ## 31                       apul-mir-novel-33
    ## 32                       apul-mir-novel-34
    ## 33                        apul-mir-novel-7
    ## 34                        apul-mir-novel-2
    ## 35                       apul-mir-novel-16
    ## 36                       apul-mir-novel-17
    ## 37                            apul-mir-100
    ## 38                           apul-mir-2036
    ## 39                       apul-mir-novel-23
    ## 40                       peve-mir-novel-25
    ## 41                       peve-mir-novel-26
    ## 42                       peve-mir-novel-27
    ## 43                       peve-mir-novel-39
    ## 44                       peve-mir-novel-28
    ## 45                       peve-mir-novel-29
    ## 46                       peve-mir-novel-40
    ## 47                        peve-mir-novel-6
    ## 48                       peve-mir-novel-30
    ## 49                       peve-mir-novel-31
    ## 50                        peve-mir-novel-2
    ## 51                       peve-mir-novel-32
    ## 52                        peve-mir-novel-1
    ## 53                       peve-mir-novel-41
    ## 54                        peve-mir-novel-7
    ## 55                        peve-mir-novel-8
    ## 56                        peve-mir-novel-9
    ## 57                        peve-mir-novel-3
    ## 58                       peve-mir-novel-34
    ## 59                           peve-mir-2025
    ## 60                       peve-mir-novel-10
    ## 61                       peve-mir-novel-11
    ## 62                       peve-mir-novel-35
    ## 63                       peve-mir-novel-36
    ## 64                           peve-mir-2036
    ## 65                       peve-mir-novel-13
    ## 66                        peve-mir-novel-4
    ## 67                            peve-mir-100
    ## 68                       peve-mir-novel-38
    ## 69                           peve-mir-2030
    ## 70  peve-mir-novel-14a; peve-mir-novel-14b
    ## 71  peve-mir-novel-14a; peve-mir-novel-14b
    ## 72                       peve-mir-novel-15
    ## 73                           peve-mir-2023
    ## 74                       peve-mir-novel-16
    ## 75                       peve-mir-novel-42
    ## 76  peve-mir-novel-17a; peve-mir-novel-17b
    ## 77  peve-mir-novel-17a; peve-mir-novel-17b
    ## 78                       peve-mir-novel-18
    ## 79                       peve-mir-novel-19
    ## 80                       peve-mir-novel-21
    ## 81                       peve-mir-novel-22
    ## 82                       peve-mir-novel-20
    ## 83                       peve-mir-novel-23
    ## 84                       peve-mir-novel-24
    ## 85                       ptuh-mir-novel-32
    ## 86                       ptuh-mir-novel-27
    ## 87                       ptuh-mir-novel-31
    ## 88                       ptuh-mir-novel-33
    ## 89                        ptuh-mir-novel-1
    ## 90                        ptuh-mir-novel-5
    ## 91                        ptuh-mir-novel-7
    ## 92                        ptuh-mir-novel-2
    ## 93                       ptuh-mir-novel-24
    ## 94                        ptuh-mir-novel-8
    ## 95                            ptuh-mir-100
    ## 96                           ptuh-mir-2036
    ## 97                        ptuh-mir-novel-6
    ## 98                        ptuh-mir-novel-4
    ## 99                        ptuh-mir-novel-3
    ## 100                          ptuh-mir-2023
    ## 101                      ptuh-mir-novel-20
    ## 102                      ptuh-mir-novel-26
    ## 103                      ptuh-mir-novel-25
    ## 104                      ptuh-mir-novel-18
    ## 105                          ptuh-mir-2030
    ## 106                      ptuh-mir-novel-10
    ## 107                      ptuh-mir-novel-17
    ## 108                          ptuh-mir-2025
    ## 109                      ptuh-mir-novel-28
    ## 110                       ptuh-mir-novel-9
    ## 111   ptuh-mir-novel-11; ptuh-mir-novel-15
    ## 112                      ptuh-mir-novel-14
    ## 113                      ptuh-mir-novel-13
    ## 114                      ptuh-mir-novel-12
    ## 115                      ptuh-mir-novel-19
    ## 116                      ptuh-mir-novel-22
    ## 117                      ptuh-mir-novel-16
    ## 118                      ptuh-mir-novel-23
    ## 119                      ptuh-mir-novel-21
    ## 120                      ptuh-mir-novel-30
    ## 121                      ptuh-mir-novel-29

``` r
write.csv(miRNA_details, "../output/20-supplementary-files/AllSpecies_Results_mature_named_miRNAs.csv", col.names = TRUE, row.names = FALSE)
```

    ## Warning in write.csv(miRNA_details,
    ## "../output/20-supplementary-files/AllSpecies_Results_mature_named_miRNAs.csv",
    ## : attempt to set 'col.names' ignored

# 5 miRanda binding information + PCC coexpression (miRNA-mRNA and miRNA-lncRNA)

## 5.1 miRNA-mRNA

``` r
consolidate_regions <- function(path_3UTR, path_CDS, path_5UTR){
# Load for each region
data_3UTR <- read.csv(path_3UTR) %>% dplyr::select(-X, -X.1)
data_CDS <- read.csv(path_CDS) %>% dplyr::select(-X)
data_5UTR <- read.csv(path_5UTR) %>% dplyr::select(-X)

# Format and combine into pooled targets
data_3UTR$region <- "3UTR"

colnames(data_CDS) <- c("miRNA", "mRNA_coord", "score", "energy", "query_start", "query_end", "subject_start", "subject_end", "total_bp_shared", "query_similar", "subject_similar", "mRNA", "PCC.cor", "p_value", "adjusted_p_value")
data_CDS$query_start_end <- paste0(data_CDS$query_start, " ", data_CDS$query_end)
data_CDS$subject_start_end <- paste0(data_CDS$subject_start, " ", data_CDS$subject_end)
data_CDS$region <- "CDS"
data_CDS <- data_CDS %>% dplyr::select(colnames(data_3UTR))

colnames(data_5UTR) <- c("miRNA", "mRNA_coord", "score", "energy", "query_start", "query_end", "subject_start", "subject_end", "total_bp_shared", "query_similar", "subject_similar", "mRNA", "PCC.cor", "p_value", "adjusted_p_value")
data_5UTR$query_start_end <- paste0(data_5UTR$query_start, " ", data_5UTR$query_end)
data_5UTR$subject_start_end <- paste0(data_5UTR$subject_start, " ", data_5UTR$subject_end)
data_5UTR$region <- "5UTR"
data_5UTR <- data_5UTR %>% dplyr::select(colnames(data_3UTR))

# Confirm all 3 region tables have same column names and orders
ifelse(identical(colnames(data_3UTR), colnames(data_CDS)) == TRUE,
       ifelse(identical(colnames(data_3UTR), colnames(data_5UTR)) == TRUE,
              data <- rbind(data_3UTR, data_CDS, data_5UTR),
              data <- NULL)
)
data
}

apul_miRanda_PCC_mRNA <- consolidate_regions(
  "../../D-Apul/output/09-Apul-mRNA-miRNA-interactions/miranda_PCC_miRNA_mRNA.csv",
  "../../D-Apul/output/09.01-Apul-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_CDS.csv",
  "../../D-Apul/output/09.01-Apul-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_5UTR.csv")

peve_miRanda_PCC_mRNA <- consolidate_regions(
  "../../E-Peve/output/10-Peve-mRNA-miRNA-interactions/Peve-miranda_PCC_miRNA_mRNA.csv",
  "../../E-Peve/output/10.01-Peve-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_CDS.csv",
  "../../E-Peve/output/10.01-Peve-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_5UTR.csv"
)

ptuh_miRanda_PCC_mRNA <- consolidate_regions(
  "../../F-Ptuh/output/11-Ptuh-mRNA-miRNA-interactions/three_prime_interaction/Ptuh-miranda_PCC_miRNA_mRNA.csv",
  "../../F-Ptuh/output/11.01-Ptuh-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_CDS.csv",
  "../../F-Ptuh/output/11.01-Ptuh-mRNA-miRNA-interactions-CDS_5UTR/miRanda_PCC_miRNA_5UTR.csv"
)

# Combine
miRanda_PCC_mRNA <- rbind(apul_miRanda_PCC_mRNA, peve_miRanda_PCC_mRNA, ptuh_miRanda_PCC_mRNA)

# Filter to only significant PCC (this consitutes our putative interaction set, which is both 
# predicted to bind via miRanda and significantly coexpressed)
miRanda_PCC_mRNA_sig <- miRanda_PCC_mRNA %>% filter(p_value < 0.05)

# Save
write.csv(miRanda_PCC_mRNA_sig, "../output/20-supplementary-files/miRanda_PCC_mRNA_sig.csv", row.names = FALSE)
write.csv(miRanda_PCC_mRNA, "../output/20-supplementary-files/miRanda_PCC_mRNA.csv", row.names = FALSE)
```

## 5.2 miRNA-lncRNA

``` r
apul_miRanda_PCC_lncRNA <- read.csv("../../D-Apul/output/28-Apul-miRNA-lncRNA-interactions/miranda_PCC_miRNA_lncRNA.csv") %>% dplyr::select(-X, -X.1)
peve_miRanda_PCC_lncRNA <- read.csv("../../E-Peve/output/15-Peve-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv") %>% dplyr::select(-X, -X.1)
ptuh_miRanda_PCC_lncRNA <- read.csv("../../F-Ptuh/output/15-Ptuh-miRNA-lncRNA-PCC/miranda_PCC_miRNA_lncRNA.csv") %>% dplyr::select(-X, -X.1)

# Ensure all three tables have matching column names and order
identical(colnames(apul_miRanda_PCC_lncRNA), colnames(peve_miRanda_PCC_lncRNA))
```

    ## [1] TRUE

``` r
identical(colnames(apul_miRanda_PCC_lncRNA), colnames(ptuh_miRanda_PCC_lncRNA))
```

    ## [1] TRUE

``` r
# Combine
miRanda_PCC_lncRNA <- rbind(apul_miRanda_PCC_lncRNA, peve_miRanda_PCC_lncRNA, ptuh_miRanda_PCC_lncRNA)

# Reorder columns to have same order as miRNA-mRNA table (above)
miRanda_PCC_lncRNA <- miRanda_PCC_lncRNA %>% dplyr::select(miRNA, lncRNA, PCC.cor, p_value, adjusted_p_value,
                                                           score, energy, query_start_end, subject_start_end, 
                                                           total_bp_shared, query_similar, subject_similar)

# Filter to only significant PCC (this constitutes our putative interaction set, which is both 
# predicted to bind via miRanda and significantly coexpressed)
miRanda_PCC_lncRNA_sig <- miRanda_PCC_lncRNA %>% filter(p_value < 0.05)

# Save
write.csv(miRanda_PCC_lncRNA_sig, "../output/20-supplementary-files/miRanda_PCC_lncRNA_sig.csv", row.names = FALSE)
write.csv(miRanda_PCC_lncRNA, "../output/20-supplementary-files/miRanda_PCC_lncRNA.csv", row.names = FALSE)
```
