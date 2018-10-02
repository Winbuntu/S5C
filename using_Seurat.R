library(data.table)



RC = as.data.frame(fread("counts.txt"))
rownames(RC) = RC$V1
RC=RC[,-1]


colnames(RC)

day.data = sapply(strsplit(  colnames(RC), split='.',fixed=TRUE),function(x) x[1])

names(day.data) = colnames(RC)
#######################
library(Seurat)

mr <- CreateSeuratObject(raw.data = RC, min.cells = 10,min.genes = 5000, is.expr = 0, project = "Pr")

mr

mr <- AddMetaData(object = mr, metadata = day.data  , 
                  col.name = c("Day")  )


VlnPlot(object = mr, features.plot = c("nGene", "nUMI"), 
        nCol = 3)

#################

mr <- FilterCells(object = mr, subset.names = c("nUMI"), 
                  low.thresholds = c(500000), high.thresholds = c(Inf))

mr


#############################

mr <- NormalizeData(object = mr, normalization.method = "LogNormalize", 
                    scale.factor = 10000)

mr <- FindVariableGenes(object = mr, mean.function = ExpMean,
                        dispersion.function = LogVMR, 
                        x.low.cutoff = 0.125, x.high.cutoff = 3, y.cutoff = 1,do.plot = T)

length(mr@var.genes)


#############

mr <- ScaleData(object = mr, vars.to.regress = c("nUMI"))


lineages = read.csv("mmc3-2.csv",header = T, stringsAsFactors = F)$ACE

mr <- RunPCA(object = mr, do.print = TRUE, pcs.print = 1:5, 
             genes.print = 5,pcs.compute = 50)


PCAPlot(object = mr, dim.1 = 1, dim.2 = 2, group = "Day")

PCAPlot(object = mr, dim.1 = 1, dim.2 = 2)




mr <- JackStraw(object = mr, num.replicate = 100, display.progress = FALSE)

JackStrawPlot(object = mr, PCs = 1:12)

mr <- RunTSNE(object = mr, dims.use = 1:10, genes.use = lineages,do.fast = TRUE,perplexity=10)

TSNEPlot(object = mr, group = "Day")

mr <- CellCycleScoring(object = mr, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes, 
                       set.ident = TRUE)


PCAPlot(object = mr, dim.1 = 1, dim.2 = 2)

FeaturePlot(mr, features.plot = c("CDX2"), reduction.use = "pca")

RidgePlot(object = mr, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
          nCol = 2)

VlnPlot(object = mr, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
          nCol = 2, group.by = "Day")


VlnPlot(object = mr, features.plot = c("GATA2","GATA3","DPPA5"), 
        nCol = 2, group.by = "Phase")
#######################################

mr.G2m = SubsetData(mr, cells.use =  rownames(mr@meta.data)[mr@meta.data$Phase=="G2M" ], do.scale = T, do.center = T)

VlnPlot(object = mr.G2m, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
        nCol = 2, group.by = "Day")

mr.S = SubsetData(mr, cells.use =  rownames(mr@meta.data)[mr@meta.data$Phase=="S" ], do.scale = T, do.center = T)
VlnPlot(object = mr.S, features.plot = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
        nCol = 2, group.by = "Day")



##########################################

mr.E67 = SubsetData(mr, cells.use =  rownames(mr@meta.data)[mr@meta.data$Day=="E6" | mr@meta.data$Day=="E7"], do.scale = T, do.center = T)


mr.E67 <- CellCycleScoring(object = mr.E67, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes, 
                       set.ident = TRUE)

PCAPlot(object = mr.E67, dim.1 = 1, dim.2 = 2)




FeaturePlot(mr, features.plot = c("TOP2A","MCM6","MKI67"), reduction.use = 'pca', pt.size = 2)

##########################

mr.E34 = SubsetData(mr, cells.use =  rownames(mr@meta.data)[mr@meta.data$Day=="E3" | mr@meta.data$Day=="E4"], do.scale = T, do.center = T)


mr.E34 <- CellCycleScoring(object = mr.E34, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes, 
                           set.ident = TRUE)

PCAPlot(object = mr.E34, dim.1 = 1, dim.2 = 2)


