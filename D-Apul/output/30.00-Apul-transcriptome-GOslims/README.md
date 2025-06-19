`deep-dive-expression/D-Apul/output/30.00-Apul-transcriptome-GOslims`

---

- [`Apulchra-expressed-genes.blastx.outfmt6`](./Apulchra-expressed-genes.blastx.outfmt6): DIAMOND BLASTx output of expressed genes against the UniProt database.

- [`Apulchra-genes.bed`](./Apulchra-genes.bed): BED file containing the genomic coordinates of _A.pulchra_ genes.

- [`Apulchra-genes.fasta.fai`](./Apulchra-genes.fasta.fai): Index file for the _A.pulchra_ gene FASTA file.

- [`gene-SPIDs-GOIDs.tsv`](./gene-SPIDs-GOIDs.tsv): Tab-separated file mapping genes to SPIDs and GO IDs.

- [`gene-SPIDs.txt`](./gene-SPIDs.txt): List of genes with their associated SPIDs.

- [`GOslim-counts.tsv`](./GOslim-counts.tsv): Counts of GO terms assigned to  GOslims associated with expressed _A.pulchra_ genes.

- [`SPIDs.txt`](./SPIDs.txt): List of unique SPIDs (Sequence Protein Identifiers) associated with expressed genes. Used as input to [`../../../M-multi-species/code/uniprot-retrieval.py`](../../../M-multi-species/code/uniprot-retrieval.py).

- [`uniprot-retrieval.tsv`](./uniprot-retrieval.tsv): Tab-separated file containing UniProt information for the SPIDs associated with expressed _A.pulchra_ genes, including gene names, descriptions, and GO terms.