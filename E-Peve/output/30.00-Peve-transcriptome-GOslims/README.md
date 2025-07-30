`deep-dive-expression/E-Peve/output/30.00-Peve-transcriptome-GOslims`

All files in this directory are generated from [`30.00-Peve-transcriptome-GOslims.Rmd`](../../code/30.00-Peve-transcriptome-GOslims.Rmd), annotating expressed genes with Gene Ontology (GO) slim terms.

---

- [`Pevermanni-expressed-genes.blastx.outfmt6`](./Pevermanni-expressed-genes.blastx.outfmt6): DIAMOND BLASTx output of expressed genes against the UniProt database.

- [`Pevermanni-genes.bed`](./Pevermanni-genes.bed): BED file containing the genomic coordinates of _P.evermanni_ genes.

- [`Pevermanni-genes.fasta.fai`](./Pevermanni-genes.fasta.fai): Index file for the _P.evermanni_ gene FASTA file.

- [`gene-SPIDs-GOIDs.tsv`](./gene-SPIDs-GOIDs.tsv): Tab-separated file mapping genes to SPIDs and GO IDs.

- [`gene-SPIDs.txt`](./gene-SPIDs.txt): List of genes with their associated SPIDs.

- [`GOslim-counts.tsv`](./GOslim-counts.tsv): Counts of GO terms assigned to  GOslims associated with expressed _P.evermanni_ genes.

- [`SPIDs.txt`](./SPIDs.txt): List of unique SPIDs (Sequence Protein Identifiers) associated with expressed genes. Used as input to [`../../../M-multi-species/code/uniprot-retrieval.py`](../../../M-multi-species/code/uniprot-retrieval.py).