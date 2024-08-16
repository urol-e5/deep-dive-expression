# expression

Analyzing gene expression and gene/ncRNA coexpression in three species of stony coral. This work builds on a completed summary of the ncRNA landscape of these species ([deep-dive](https://github.com/urol-e5/deep-dive)).

## Specific sub-efforts

As in [deep-dive](https://github.com/urol-e5/deep-dive), this expression/co-expression work will focus on the following three species of stony corals:

D)  Acropora pulchra

E)  Porites evermanni

F)  Pocillopora tuahiniensis

## How to work in this repo
### (file structure)

Top level directories are associated with each sub-effort categorized by species. For instance:

```         
A-Pver
B-Mcap
C-Pacu
D-Apul
E-Peve
F-Ptua
```

Within each top level directory there should be 3 directories:

```         
data
code
output
```

Any document in `code` should start with a 2 number prefix (e.g., 01-methylation-explore.Rmd). All output from that code should be in a sub-directory of `output` named the same as the code. For example the output of **01-methylation-explore.Rmd** would be in **A-pver/output/01-methylation-explore/**.

Please use **Relative Paths**. Commit and Push often.

Links to other data types (e.g. FastQs, BAMs) can be found in the [project wiki](https://github.com/urol-e5/deep-dive/wiki).

## More
### Genomes of interest

All genomes of interest can be found in our [species descriptions and genomic resources wiki page](https://github.com/urol-e5/deep-dive/wiki/Species-Characteristics-and-Genomic-Resources).
