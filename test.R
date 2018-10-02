library(Seurat)

# Read in the expression matrix The first row is a header row, the first
# column is rownames
exp.mat <- read.table(file = "./cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt", 
                      header = TRUE, as.is = TRUE, row.names = 1)


exp.mat.norm = read_count(exp.mat)


# Also read in a list of cell cycle markers, from Tirosh et al, 2015
#cc.genes <- readLines(con = "~/Downloads/seurat_resources/regev_lab_cell_cycle_genes.txt")

# We can segregate this list into markers of G2/M phase and markers of S
# phase
#s.genes <- cc.genes[1:43]
#g2m.genes <- cc.genes[44:97]

# Create our Seurat object and complete the initalization steps
marrow <- CreateSeuratObject(raw.data = exp.mat)
marrow <- NormalizeData(object = marrow)
marrow <- FindVariableGenes(object = marrow, do.plot = FALSE, display.progress = FALSE)
marrow <- ScaleData(object = marrow, display.progress = FALSE)

marrow <- RunPCA(object = marrow, pc.genes = marrow@var.genes, pcs.print = 1:4, 
                 genes.print = 10)


marrow <- RunTSNE(object = marrow, dims.use = 1:15, do.fast = TRUE)


TSNEPlot(marrow)

marrow <- FindClusters(marrow)


TSNEPlot(marrow)


FeaturePlot(marrow, features.plot = c("GATA1","CD48","MPL"), reduction.use = "tsne")

VlnPlot(marrow, features.plot = c("GATA1","CD48","MPL"))

marrow <- CellCycleScoring(object = marrow, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes, 
                           set.ident = F)


VlnPlot(object = marrow, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
        nCol = 2)

marrow.S = SubsetData(marrow, cells.use =  rownames(marrow@meta.data)[marrow@meta.data$Phase=="S" ], do.scale = T, do.center = T)


marrow <- RunPCA(object = marrow, pc.genes = c(cc.genes$s.genes, cc.genes$g2m.genes), do.print = FALSE)
PCAPlot(object = marrow)

plot(marrow@dr$pca@cell.embeddings[,1], marrow@dr$pca@cell.embeddings[,2], col = factor(assignments) )

#################################################

marrow.2 = Add_info(marrow, genes.list = list(cc.genes$s.genes))

cor(
  genes.scores.use[,1],

marrow@meta.data$S.Score,method = "spearman"
)


