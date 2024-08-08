# epaclades
Clade-level summary of results of the Evolutionary Placement Algorithm.

**epaclades** contains a set of R functions (R Core Team 2024) post-processing outputs from a phylogenetic placement analysis, performed by *Evolutionary Placement Algorithm* (EPA, Berger et al. 2011) or *pplacer* (Matsen et al. 2010).
They summarize probability of a query sequence placement on separate branches of the tree, so the support for its placement within a whole clade (or lineage) is obtained.

The functions perform successively three basic tasks:
1. read_jplace imports the output file in .jplace format (Matsen et al. 2012)
2. classify_jplace classifies branches as belonging to clades, which are defined in a suitable data frame
3. classify_sequences calculates probabilities of query sequences belonging to pre-specified clades

In addition, root_jplace roots the tree (and modifies branch labels accordingly), using specified outgroup sequences

The functions depend on R package ape (Paradis & Schliep 2019). 'plot.jplace' can save the plot as a .pdf file, the output of 'classify_sequences' can be saved as a tab-delimited file.

### **Package installation**
It can be installed using devtools:

```
library(devtools)
devtools::install_github("onmikula/epaclades", dependencies = TRUE)
```

**References:**
- Berger SA, Krompass D, Stamatakis A (2011) Performance, accuracy, and web server for evolutionary placement of short sequence reads under maximum likelihood. Systematic Biology, 60: 291–302. https://doi.org/10.1093/sysbio/syr010
- Matsen FA, Hoffman NG, Gallagher A, Stamatakis A (2012) A format for phylogenetic placements. PLoS ONE 7: e31009. https://doi.org/10.1371/journal.pone.0031009
- Matsen FA, Kodner R, Armbrust E (2010) pplacer: linear time maximum-likelihood and Bayesian phylogenetic placement of sequences onto a fixed reference tree. BMC Bioinformatics 11: 538. https://doi.org/10.1186/1471-2105-11-538
- Paradis E, Schliep K (2019) ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics 35, 526–528. https://doi.org/10.1093/ bioinformatics/bty633
- R Core Team (2024) R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org
