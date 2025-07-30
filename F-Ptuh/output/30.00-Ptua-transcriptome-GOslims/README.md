`deep-dive-expression/F-Ptuh/output/30.00-Ptua-transcriptome-GOslims`

All files in this directory are generated from [`30.00-Ptua-transcriptome-GOslims.Rmd`](../../code/30.00-Ptua-transcriptome-GOslims.Rmd), annotating expressed genes with Gene Ontology (GO) slim terms.

---

- [`Ptuahiniensis-expressed-genes.blastx.outfmt6`](./Ptuahiniensis-expressed-genes.blastx.outfmt6): DIAMOND BLASTx output of expressed genes against the UniProt database.

- [`Ptuahiniensis-genes.bed`](./Ptuahiniensis-genes.bed): BED file containing the genomic coordinates of _P.tuahiniensis_ genes.

- [`Ptuahiniensis-genes.fasta.fai`](./Ptuahiniensis-genes.fasta.fai): Index file for the _P.tuahiniensis_ gene FASTA file.

- [`gene-SPIDs-GOIDs.tsv`](./gene-SPIDs-GOIDs.tsv): Tab-separated file mapping genes to SPIDs and GO IDs.

- [`gene-SPIDs.txt`](./gene-SPIDs.txt): List of genes with their associated SPIDs.

- [`GOslim-counts.tsv`](./GOslim-counts.tsv): Counts of GO terms assigned to  GOslims associated with expressed _P.tuahiniensis_ genes.

- [`SPIDs.txt`](./SPIDs.txt): List of unique SPIDs (Sequence Protein Identifiers) associated with expressed genes. Used as input to [`../../../M-multi-species/code/uniprot-retrieval.py`](../../../M-multi-species/code/uniprot-retrieval.py).