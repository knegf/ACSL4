1. ACSL4-programmed sphingomyelin metabolism rejuvenates tumor-infiltrating T-cells to enhance efficacy of anti-PD-1 immunotherapy

2. [GitHub all releases](https://github.com/NiEntropy/CYM)

3. This repository contains the source code for the scRNA-seq analysis of CD8+ and CD4+ TILs from glioblastoma (GSE154795) and squamous cell carcinoma (GSE123813) patients pre- and post-anti-PD-1 therapy.

4. For questions, please contact the corresponding author at \<1620184503@cpu.edu.cn \>.

5. Dependencies

* R-4.1.1
* R package 
  * 'Seurat', 
  * 'ggplot2',
  * 'viridisLite', 
  * 'viridis',
  * 'dplyr', 
  * 'openxlsx', 
  * 'tibble'

6. To Run

`scRNA-seq analysis_code.R`

the steps contained in the script are listed as belows:

* a. Data preparing
* b. Data integrating
* c. Dimensional reduction
* d. clustering
* e. cluster DE genes
* f. Annotation
* g. list CD4 T cells treatment DE genes
* h. list CD8 T cells treatment DE genes

