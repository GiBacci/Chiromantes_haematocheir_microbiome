# *Chiromantes haematocheir* microbiome <br/> [![Generic badge](https://img.shields.io/badge/Made_with-R_Markdown-blue.svg)](https://shields.io/) [![Generic badge](https://img.shields.io/github/license/gibacci/Chiromantes_haematocheir_microbiome)](https://shields.io/)

<img src="ACR_01_Color_Imeco_Black.png" width="106" height="60" align="right">
<img src="COMBOmod_final.png" width="115" height="60" align="right">

This repository contains all the codes used to generate files and table reported in the work:

Bacci, Fratini, Meriggi et al. (2021) [Organ-specific microbiota enhances the terrestrial lifestyle of a brachyuran crab](https://www.biorxiv.org/content/10.1101/2021.03.30.437674v2) *bioRxiv*

If you use code from this repo, please cite our paper as follows (BibTex):

```BibTeX
@article{Bacci2021.03.30.437674,
  author = {Bacci, Giovanni and Fratini, Sara and Meriggi, Niccolo and Cheng, Christine L.Y. and Ng, Ka Hei and Pindo, Massimo and Iannucci, Alessio and Mengoni, Alessio and Cavalieri, Duccio and Cannicci, Stefano},
  title = {Organ-specific microbiota enhances the terrestrial lifestyle of a brachyuran crab},
  year = {2021},
  doi = {10.1101/2021.03.30.437674},
  publisher = {Cold Spring Harbor Laboratory},	
  URL = {https://www.biorxiv.org/content/early/2021/04/01/2021.03.30.437674},
  eprint = {https://www.biorxiv.org/content/early/2021/04/01/2021.03.30.437674.full.pdf},
  journal = {bioRxiv}}
```

The main objective of this project is to unveil the evolutionary role of associated and/or symbiotic bacteria and fungi in the successful conquest of land by crabs and, ultimately, to share light on the importance of such association/symbioses for the evolution of terrestrial habits by animals in general. Additional information about all people and institutions involved in this work can be found here:

1. [Stefano Cannicci's lab](https://www.imeco-lab.com/)
2. [Duccio Cavalieri's home page](https://www.unifi.it/p-doc2-2015-0-A-2b333d2e342d-0.html)
3. [Florence Computational Biology Group](https://github.com/combogenomics)

Sequences were deposited in the European Nucleotide Archive (ENA) under the accession [PRJEB43930](https://www.ebi.ac.uk/ena/browser/view/PRJEB43930).

## Project abstract

The transition to terrestrial environments by formerly aquatic species has occurred repeatedly in all animal groups and lead to the vast biodiversity of terrestrial forms that we are observing nowadays. The differences between aquatic and terrestrial habitats are enormous and nearly every facet of life was altered to cope with the water-to-land transition. In morphological and physiological terms, animals belonging to rather different and unrelated taxa are known to share very similar adaptations to land, but nothing is known about the critical role of commensal and symbiotic microbiota in such transition. The holobiont theory of evolution indicates that evolution acts on both the genes of the host and of the commensal microbiota, termed microbiome. Terrestrial crabs are a perfect model to study the evolutionary pathways and the ecological role of animal-microbiome symbioses (the so called symbiomes), since the evolution of crabs towards the conquest of terrestrial environments is happening right now through a number of non-related groups.

## File description

Microbial counts and taxonomic assignments are stored in `16s` (Bacteria/Archaea) and `its` (Fungi) folders:

  - `track_back.tsv` : number of reads retained after each step of analysis.
  - `seqTabs/all_seqtab_nochim.rds` : raw microbial counts.
  - `seqTabs/tax_assignments.rds` : taxonomic classification of amplicon sequence variants detected.

All figures and tables included in the paper are saved in the `Figures` and `Tables` folder, with the following modifications:

  - Figure 1: this map was manually created and it can not be reproduced within the R environment. A copy of the image is provided in the `Figures` folder.
  - Figure 5: stripes were manually added and are not automatically generated.
  - Figure 6: venn diagrams were manually aligned to the panels.
  - Table S1: the table was manually compiled. A copy of the file is provided in the `Tables` folder.

Files generated during the analysis are saved in the `outputs` folder. To repeat one or more steps simply delete output files and launch the pipeline again.

Scripts used are reported in the `scripts` folder and launched in the following order:

  1. `formatDataAndCorrelations.R` : raw data format and correlation between fungal replicates.
  2. `alphaBetaDiversity.R` : alpha and beta diversity analyses.
  3. `deseq2Clustering.R` : amplicon sequence variant clustering and likelihood-ratio test.
      - `hypergeometricTest.R` : hypergeometric test used for variant enrichment.
      - `sunburstPlot.R` : sunburst plot generation.
  4. `geneEnrichment.R` : gene enrichment in variant clusters.

Data coming from other resources are available in the `data` folder:

  - `ec2go.txt` : map file for Enzyme Commission numbers (EC numbers).
  - `gene_asv_matrix.rds` : abundance of GO terms in all sequence variants detected.
  - `gene_count.rds` : number of copies for each GO term detected in each cluster.

## Instructions

The file `report.Rmd` can be rendered into `pdf` by using the `knitr` plugin integrated into [RStudio](https://rstudio.com/?_ga=2.50552553.1339302526.1611745574-1183453795.1578408315). Alternatively, one can directly render documents from terminal by launching this command:

```shell
R -e "rmarkdown::render('report.Rmd', output_file='report.pdf')"
```

It is important to note that `pdf` documents are generated using `pandoc` package integrated in RStudio. In Debian/Linux systems the library is usually in `/usr/lib/rstudio/bin/pandoc/` and could differ from the version installed at system level. If during the rendering you get a message like:

```
Error: pandoc version 1.XY.Z or higher is required and was not found (see the help page ?rmarkdown::pandoc_available).
```

then you probably need to tell `R` where the `pandoc` library integrated with RStudio is located in your system. You can do this with the command: 

```shell
R -e "Sys.setenv(RSTUDIO_PANDOC='/usr/lib/rstudio/bin/pandoc/'); rmarkdown::render('report.Rmd', output_file='report.pdf')"
```
