# *Chiromantes haematocheir* microbiome <br/> [![Generic badge](https://img.shields.io/badge/Made_with-R_Markdown-blue.svg)](https://shields.io/) [![Generic badge](https://img.shields.io/github/license/gibacci/Chiromantes_haematocheir_microbiome)](https://shields.io/)

<img src="ACR_01_Color_Imeco_Black.png" width="106" height="60" align="right">
<img src="COMBOmod_final.png" width="115" height="60" align="right">

This repository contains all the codes used to generate files and table reported in the work:

Bacci, Fratini, Meriggi et al. (2021) [Organ-specific microbiota enhances the terrestrial lifestyle of a brachyuran crab](https://www.biorxiv.org/content/10.1101/2021.03.30.437674v1) *bioRxiv*

If you use code from this repo, please cite our paper as follows (BibTex):

```BibTeX
@article{Bacci2021.03.30.437674,
  author = {Bacci, Giovanni and Fratini, Sara and Meriggi, Niccolo and Cheng, Christine L.Y. and Ng, Ka Hei and Iannucci, Alessio and Mengoni, Alessio and Cavalieri, Duccio and Cannicci, Stefano},
  title = {Organ-specific microbiota enhances the terrestrial lifestyle of a brachyuran crab},
  year = {2021},
  doi = {10.1101/2021.03.30.437674},
  publisher = {Cold Spring Harbor Laboratory},	
  URL = {https://www.biorxiv.org/content/early/2021/03/31/2021.03.30.437674},
  eprint = {https://www.biorxiv.org/content/early/2021/03/31/2021.03.30.437674.full.pdf},
  journal = {bioRxiv}}
```

The main objective of this project is to unveil the evolutionary role of associated and/or symbiotic bacteria and fungi in the successful conquest of land by crabs and, ultimately, to share light on the importance of such association/symbioses for the evolution of terrestrial habits by animals in general. Additional information about all people and institutions involved in this work can be found here:

1. [Stefano Cannicci's lab](https://www.imeco-lab.com/)
2. [Duccio Cavalieri's home page](https://www.unifi.it/p-doc2-2015-0-A-2b333d2e342d-0.html)
3. [Florence Computational Biology Group](https://github.com/combogenomics)

## Project abstract

The transition to terrestrial environments by formerly aquatic species has occurred repeatedly in all animal groups and lead to the vast biodiversity of terrestrial forms that we are observing nowadays. The differences between aquatic and terrestrial habitats are enormous and nearly every facet of life was altered to cope with the water-to-land transition. In morphological and physiological terms, animals belonging to rather different and unrelated taxa are known to share very similar adaptations to land, but nothing is known about the critical role of commensal and symbiotic microbiota in such transition. The holobiont theory of evolution indicates that evolution acts on both the genes of the host and of the commensal microbiota, termed microbiome. Terrestrial crabs are a perfect model to study the evolutionary pathways and the ecological role of animal-microbiome symbioses (the so called symbiomes), since the evolution of crabs towards the conquest of terrestrial environments is happening right now through a number of non-related groups.

## File description

Microbial counts and taxonomic assignments are stored in `16s` (Bacteria/Archae) and `its` (Fungi) folders:

  - `track_back.tsv` : number of reads retained after each step of analysis
  - `seqTabs/all_seqtab_nochim.rds` : raw microbial counts
  - `seqTabs/tax_assignments.rds` : taxonomic classification of amplicon sequence variants detected

All figures and tables included in the paper are saved in the `Figures` and `Tables` folder, with the following modifications:

  - Figure 1: this map was manually created and it can not be reproduced within the R environment. A copy of the image is provided in the `Figures` folder
  - Figure 5: stripes were manually added and are not automatically generated
  - Figure 6: venn diagrams were manually aligned to the panels
  - Table S1: the table was manually compiled. A copy of the file is provided in the `Tables` foder

Files generated during the analysis are saved in the `outputs` folder. To repeat one or more steps simply delete output files and launch the pipeline again.

Scripts used in the pipeline are reported in the `scripts` folder and launched in the following order:

  1. `formatDataAndCorrelations.R` : raw data format and correlation between fungal replictaes
  2. `alphaBetaDiversity.R` : alpha and beta diversity analyses
  3. `deseq2Clustering.R` : amplicon sequence variant clustering and likelihood-ratio test
  4. `geneEnrichment.R` : gene enrichemnt in variant clusters


