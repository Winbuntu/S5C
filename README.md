# S5C
**Single-Cell Cluster based Cell-Cycle Correction**

Cell-cycle difference between single cells can be a source of heterogenetiy, 
which masked the interested variation in single cell RNA-seq data. Here we developed Single-Cell Cluster based Cell-Cycle Correction (S5C), 
a package optimized for non-droplet based scRNA-seq data, which seqeunce a fewer cells 
 (10~1000) but generate ~1M reads and cover ~10K genes per cell, such as SMART-seq, SMART-seq2, CEL-seq, C1 fluidigm, etc. 

Compared to previous packages, such as ccRemover and Seurat that can also remove unwanted heterogeneity from cell cycle from 
single cell RNA-seq data, S5C has following features makes it more suitable for non-droplet based scRNA-seq protocols:

- Can take not only read counts, but also any other properly normalized metrics as input, including RPM, TPM and RPKM, 
which is generated from software commenly used for non-droplet based scRNA-seq data such as Kallisto and Salmon 
(Up to today Seurat only accepts read count as input). 

- Automatically detect outlier cells and exclude them when computing the coefficient, 
which improve the robustness of "regressing-out" cell cycle bias, especially when the cell number is small 
(ccRemover do not remove outliers; Seurat rely on large cell numbers to metigate the influence of outliers). 

- Parameters were optimized for non-droplet based scRNA-seq data, 
which generate \~1M reads and cover \~10K genes from 10~1000 cells 
(Seurat were designed for droplet-based methods, therefore more suitable for large cell number (more than 1k), and fewer reads per sample (10k-50k) ).
