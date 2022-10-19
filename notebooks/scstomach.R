---
title: "Paper_Figures_plotting"
author: "Peng Zhang"
date: "Oct 24th,2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## R Markdown

## preprocessing data and annotating the cell clusters

```{r preprocessing and summary of data, eval= T, echo=T,warning=F}

library(dplyr)
library(reshape2)
library(Seurat)

## Not Run
data_temp <- Read10X("/heartdata8t_A/zhangpeng/Final_data_all/raw_gene_bc_matrices_mex/hg19")
dt <- CreateSeuratObject(raw.data = data_temp, min.cells = 3, min.genes = 300,project = "10X_IM")

# # annotating the batch effect
dt@meta.data$batch <- "NAG1"
dt@meta.data$batch[grepl('-10',dt@cell.names)] <- 'IMS3'
dt@meta.data$batch[grepl('-12',dt@cell.names)] <- 'CAN1'
dt@meta.data$batch[grepl('-13',dt@cell.names)] <- 'EGC'
dt@meta.data$batch[grepl('-11',dt@cell.names)] <- 'IMS4'
dt@meta.data$batch[grepl('-9',dt@cell.names)] <- 'IMS2'
dt@meta.data$batch[grepl('-6',dt@cell.names)] <- 'IMW1'
dt@meta.data$batch[grepl('-7',dt@cell.names)] <- 'IMW2'
dt@meta.data$batch[grepl('-4',dt@cell.names)] <- 'CAG3'
dt@meta.data$batch[grepl('-5',dt@cell.names)] <- 'NAG2'
dt@meta.data$batch[grepl('-3',dt@cell.names)] <- 'CAG2'
dt@meta.data$batch[grepl('-2',dt@cell.names)] <- 'CAG1'
dt@meta.data$batch[grepl('-8',dt@cell.names)] <- 'IMS1'

valid.cells <- dt@cell.names[-which(dt@meta.data$batch == "CAN1")]
dt <- SubsetData(dt,cells.use = valid.cells,do.clean = T)


## Ribosomal & Mitochondrial genes
ribo.genes <- grep(pattern = "^RPL", x = rownames(x = dt@data), value = TRUE)
mito.genes <- grep(pattern = "^MT-", x = rownames(x = dt@data), value = TRUE)
percent.ribo <- Matrix::colSums(dt@raw.data[ribo.genes, ])/Matrix::colSums(dt@raw.data)
percent.mito <- Matrix::colSums(dt@raw.data[mito.genes, ])/Matrix::colSums(dt@raw.data)
dt <- AddMetaData(object = dt, metadata = percent.ribo, col.name = "percent.ribo")
dt <- AddMetaData(object = dt, metadata = percent.mito, col.name = "percent.mito")

# # Selecting the valid genes and UMIs
dt <- FilterCells(object = dt, subset.names = c("nGene", "percent.mito","percent.ribo"),
                  low.thresholds = c(400, -Inf, -Inf), high.thresholds = c(7000, 0.2, 0.2))
dt <- NormalizeData(object = dt, normalization.method = "LogNormalize",
                    scale.factor = 10000)
dt <- FindVariableGenes(object = dt, mean.function = ExpMean, dispersion.function = LogVMR,
                        x.low.cutoff = 0.05, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
dt <- ScaleData(object = dt, do.scale = T,do.center = T,vars.to.regress = c("nUMI", "percent.mito","percent.ribo","batch"))
dt <- RunPCA(object = dt, pc.genes = dt@var.genes, do.print = TRUE, pcs.print = 1:5,
             genes.print = 5,pcs.compute = 120)
dt <- FindClusters(object = dt, reduction.type = "pca", dims.use = 1:50,
                   resolution = 2, print.output = 0, save.SNN = TRUE)
dt <- RunTSNE(object = dt, dims.use = 1:50, do.fast = TRUE)
# }


save(dt,file = "data.temp.res.2.Rdata")

load("data.temp.res.2.Rdata")


pdf(file = "SupplementaryFigure1-1.pdf", width = 20, height = 20)
FeaturePlot(object = dt, features.plot = c("EPCAM","VIM","MUC6","MUC5AC","PGA4","CHGA","MKI67","CEACAM6","OLFM4","FABP1","MUC2","TFF3"), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 0.5)
dev.off()


pdf(file = "SupplementaryFigure1-2.pdf",width = 20, height = 20)
FeaturePlot(object = dt, features.plot = c("CD79A","CEACAM5","CSF1R","TPSAB1","CD68","DCN","ACTA2","VWF","CD4","CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 0.5)
dev.off()
pdf(file = "SupplementaryFigure1-6.pdf", width = 20, height = 20)
FeaturePlot(object = dt, features.plot = c("HES6","MUC2","SPINK4","SOX4","KLK10","ITLN1"), cols.use = c("grey", "blue"), reduction.use = "tsne",pt.size = 0.5)
dev.off()

all.cluster.markers <- FindAllMarkers(dt,logfc.threshold = log2(1.5),only.pos = T)
# save(all.cluster.markers,file = "all.cluster.marker.Rdata")
write.csv(all.cluster.markers,file = "all.markers.csv")

# load("data.temp.res.2.Rdata")

{
current.cluster.ids <- 0:38
new.cluster.ids <- c("PMC","PMC","T cell","PMC","Stem-like cell","Enterocyte","PMC","GMC","Enteroendocrine","PC","PMC","B cell","PMC","B cell","Neck cell","Fibroblast","Enterocyte","PMC","Cancer cell","Enterocyte","Enteroendocrine","EC","Goblet cell","PMC","Enterocyte","Fibroblast","PMC","Enteroendocrine","Macrophage","PMC","PMC","Enteroendocrine","SM cell","EC","PMC","Mast cell","GMC","Chief cell","NA")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}

pdf(file = "SupplementaryFigure1-4.pdf", width = 20, height = 20)
TSNEPlot(object = dt,group.by = "batch")
dev.off()


{
current.cluster.ids <- c(2,5,8,60,68,33,34,75,66,44,65,50,79,7,9,64)
new.cluster.ids <- c("B cell","T cell","Enteroendocrine","Fibroblast","SMC","Endothelial cell","Goblet cell","Chief cell","Mast cell","Macrophage","GMC","GMC","GMC","GMC","GMC","GMC")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}


{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}


save(dt, file = "data.temp.res.2.annotation.Rdata")


## Supplementary figures for batch effects across multiple samples
NAG1 <- dt@cell.names[which(dt@meta.data$batch == "NAG1")]
NAG1 <- SubsetData(dt,cells.use = NAG1)
pdf("S5-1.pdf",width = 5,height = 5)
TSNEPlot(object = NAG1,pt.size = 1,do.label = F) 
dev.off()


NAG2 <- dt@cell.names[which(dt@meta.data$batch == "NAG2")]
NAG2 <- SubsetData(dt,cells.use = NAG2)
pdf("S5-2.pdf",width = 5,height = 5)
TSNEPlot(object = NAG2,pt.size = 1,do.label = F) 
dev.off()

CAG1 <- dt@cell.names[which(dt@meta.data$batch == "CAG1")]
CAG1 <- SubsetData(dt,cells.use = CAG1)
pdf("S5-3.pdf",width = 5,height = 5)
TSNEPlot(object = CAG1,pt.size = 1,do.label = F) 
dev.off()


CAG2 <- dt@cell.names[which(dt@meta.data$batch == "CAG2")]
CAG2 <- SubsetData(dt,cells.use = CAG2)
pdf("S5-4.pdf",width = 5,height = 5)
TSNEPlot(object = CAG2,pt.size = 1,do.label = F) 
dev.off()


CAG3 <- dt@cell.names[which(dt@meta.data$batch == "CAG3")]
CAG3 <- SubsetData(dt,cells.use = CAG3)
pdf("S5-5.pdf",width = 5,height = 5)
TSNEPlot(object = CAG3,pt.size = 1,do.label = F) 
dev.off()

IMW1 <- dt@cell.names[which(dt@meta.data$batch == "IMW1")]
IMW1 <- SubsetData(dt,cells.use = IMW1)
pdf("S5-6.pdf",width = 5,height = 5)
TSNEPlot(object = IMW1,pt.size = 1,do.label = F) 
dev.off()

IMW2 <- dt@cell.names[which(dt@meta.data$batch == "IMW2")]
IMW2 <- SubsetData(dt,cells.use = IMW2)
pdf("S5-7.pdf",width = 5,height = 5)
TSNEPlot(object = IMW2,pt.size = 1,do.label = F) 
dev.off()


IMS1 <- dt@cell.names[which(dt@meta.data$batch == "IMS1")]
IMS1 <- SubsetData(dt,cells.use = IMS1)
pdf("S5-8.pdf",width = 5,height = 5)
TSNEPlot(object = IMS1,pt.size = 1,do.label = F) 
dev.off()

IMS2 <- dt@cell.names[which(dt@meta.data$batch == "IMS2")]
IMS2 <- SubsetData(dt,cells.use = IMS2)
pdf("S5-9.pdf",width = 5,height = 5)
TSNEPlot(object = IMS2,pt.size = 1,do.label = F) 
dev.off()

IMS3 <- dt@cell.names[which(dt@meta.data$batch == "IMS3")]
IMS3 <- SubsetData(dt,cells.use = IMS3)
pdf("S5-10.pdf",width = 5,height = 5)
TSNEPlot(object = IMS3,pt.size = 1,do.label = F) 
dev.off()

IMS4 <- dt@cell.names[which(dt@meta.data$batch == "IMS4")]
IMS4 <- SubsetData(dt,cells.use = IMS4)
pdf("S5-11.pdf",width = 5,height = 5)
TSNEPlot(object = IMS4,pt.size = 1,do.label = F) 
dev.off()

CAN2 <- dt@cell.names[which(dt@meta.data$batch == "EGC")]
CAN2 <- SubsetData(dt,cells.use = CAN2)
pdf("S5-12.pdf",width = 5,height = 5)
TSNEPlot(object = CAN2,pt.size = 1,do.label = F) 
dev.off()

```



####### Preparing the data ########

#### isolating the raw data for each cell type
```{r isolating the raw data for each cell type}

load("data.temp.res.2.annotation.Rdata")

## Epithelial cells all (Including intestinal cell lineages)
Epithelial.cells.all <- SubsetData(dt,ident.use = c("PMC","GMC","Stem-like cell","Enterocyte","Enteroendocrine","PC","Neck cell","Chief cell","Cancer cell","Goblet cell"),do.center = T,do.scale = T)
  

## Epithelial cells (Excluding intestinal cell lineages)
Epithelial.cells <- SubsetData(dt,ident.use = c("PMC","GMC","Stem-like cell","PC","Neck cell","Chief cell","Cancer cell","Enterocyte","Enteroendocrine","Goblet cell"),do.center = T,do.scale = T)
Epithelial.cells@meta.data$orig.ident <- Epithelial.cells@ident
Epithelial.cells <- FindVariableGenes(object = Epithelial.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = Epithelial.cells@var.genes)
Epithelial.cells <- RunPCA(object = Epithelial.cells,pc.genes = Epithelial.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
Epithelial.cells <- FindClusters(object = Epithelial.cells, reduction.type = "pca", dims.use = 1:20,
    resolution = c(2), print.output = 0, save.SNN = TRUE)
Epithelial.cells <- RunTSNE(object = Epithelial.cells, dims.use = 1:20, do.fast = TRUE)




## IIC epithelial cells vs cancer cells
Epithelial.cells.IIC <- SubsetData(dt,ident.use = c("PMC","GMC"),do.center = T,do.scale = T)
Epithelial.cells.IIC@meta.data$orig.ident <- Epithelial.cells.IIC@ident
Epithelial.cells.IIC <- FindVariableGenes(object = Epithelial.cells.IIC, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = Epithelial.cells.IIC@var.genes)
Epithelial.cells.IIC <- RunPCA(object = Epithelial.cells.IIC,pc.genes = Epithelial.cells.IIC@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
Epithelial.cells.IIC <- FindClusters(object = Epithelial.cells.IIC, reduction.type = "pca", dims.use = 1:20,
    resolution = c(2), print.output = 0, save.SNN = TRUE)
Epithelial.cells.IIC <- RunTSNE(object = Epithelial.cells.IIC, dims.use = 1:20, do.fast = TRUE)


GMC.PMC.markers <- FindMarkers(dt,ident.1 = "PMC",ident.2 = "GMC",logfc.threshold = log2(1.5))
write.csv(GMC.PMC.markers,file = "GMC.PMC.csv",sep = "\t")

## PMC cells
PMC.cells  <- SubsetData(dt,ident.use = c("PMC"),do.center = T,do.scale = T)
PMC.cells  <- FindVariableGenes(object = PMC.cells , mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = PMC.cells @var.genes)
PMC.cells  <- RunPCA(object = PMC.cells ,pc.genes = PMC.cells @var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
PMC.cells  <- FindClusters(object = PMC.cells , reduction.type = "pca", dims.use = 1:20,
    resolution = c(1), print.output = 0, save.SNN = TRUE)
PMC.cells  <- SubsetData(PMC.cells,ident.use = c(0:9,10,11,13),do.center = T,do.scale = T)
PMC.cells  <- FindVariableGenes(object = PMC.cells , mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = PMC.cells @var.genes)
PMC.cells  <- RunPCA(object = PMC.cells ,pc.genes = PMC.cells @var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
PMC.cells  <- FindClusters(object = PMC.cells , reduction.type = "pca", dims.use = 1:15,
    resolution = c(1), print.output = 0, save.SNN = TRUE)
PMC.cells  <- RunTSNE(object = PMC.cells , dims.use = 1:15, do.fast = TRUE)

{
current.cluster.ids <- c(0,1,2,3,11,14,5,4,6,7,10,8,9,12,13)
new.cluster.ids <- c(rep("CAG",7),rep("NAG",4),rep("IM",3),"Can")
PMC.cells@ident <- plyr::mapvalues(x = PMC.cells@ident, from = current.cluster.ids, to = new.cluster.ids)
}



## GMC cells
GMC.cells.1 <- Epithelial.cells@cell.names[which(Epithelial.cells@ident == 22)]
GMC.cells.2 <- dt@cell.names[which(dt@ident == "GMC")]
GMC.cells <- union(GMC.cells.1,GMC.cells.2)
GMC.cells  <- SubsetData(Epithelial.cells,cells.use = GMC.cells,do.center = T,do.scale = T)
GMC.cells  <- FindVariableGenes(object = GMC.cells , mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.5, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = GMC.cells @var.genes)
GMC.cells  <- RunPCA(object = GMC.cells ,pc.genes = GMC.cells @var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
GMC.cells  <- FindClusters(object = GMC.cells , reduction.type = "pca", dims.use = 1:20,
    resolution = c(1), print.output = 0, save.SNN = TRUE)
GMC.cells  <- RunTSNE(object = GMC.cells , dims.use = 1:20, do.fast = TRUE)





## cells
PC.cells  <- SubsetData(dt,ident.use = "PC",do.center = T,do.scale = T)
PC.cells  <- FindVariableGenes(object = PC.cells , mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = PC.cells @var.genes)
PC.cells  <- RunPCA(object = PC.cells ,pc.genes = PC.cells @var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
PC.cells  <- FindClusters(object = PC.cells , reduction.type = "pca", dims.use = 1:20,
    resolution = c(1), print.output = 0, save.SNN = TRUE)
PC.cells  <- RunTSNE(object = PC.cells , dims.use = 1:20, do.fast = TRUE)



## H.pylori.infections
H.pylori.cells <- dt@cell.names[which(dt@meta.data$batch %in% c("IMW1","IMW2","IMS1","IMS2","IMS3","IMS4"))] 
  H.pylori.cells <- SubsetData(dt,cells.use = H.pylori.cells,do.center = T,do.scale = T)
H.pylori.cells <- SubsetData(H.pylori.cells,ident.use = c("PMC","GMC"))
H.pylori.cells@meta.data$orig.ident <- "Hp+"
H.pylori.cells@meta.data$orig.ident[which(H.pylori.cells@meta.data$batch %in% c("IMW1","IMS1","IMS2"))] <- "Hp-"
H.pylori.cells  <- FindVariableGenes(object = H.pylori.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = H.pylori.cells@var.genes)
H.pylori.cells  <- RunPCA(object = H.pylori.cells ,PC.genes = SLC.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
H.pylori.cells  <- FindClusters(object = H.pylori.cells, reduction.type = "pca", dims.use = 1:20,
    resolution = c(3), print.output = 0, save.SNN = TRUE)
H.pylori.cells <- SubsetData(H.pylori.cells,ident.remove = 5,do.center = T,do.scale = T)
H.pylori.cells  <- FindVariableGenes(object = H.pylori.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = H.pylori.cells@var.genes)
H.pylori.cells  <- RunPCA(object = H.pylori.cells ,PC.genes = SLC.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
H.pylori.cells  <- FindClusters(object = H.pylori.cells, reduction.type = "pca", dims.use = 1:20,
    resolution = c(3), print.output = 0, save.SNN = TRUE)
H.pylori.cells  <- RunTSNE(object = H.pylori.cells, dims.use = 1:20, do.fast = TRUE)


## Goblet cells
goblet.cell <- SubsetData(dt,ident.use = "Goblet cell",do.scale  = T,do.center = T)
valid.cells <- goblet.cell@cell.names[-which(goblet.cell@meta.data$batch == "EGC")]  
goblet.cell <- SubsetData(goblet.cell,cells.use = valid.cells,do.center = T,do.scale = T) 
goblet.cell <- FindVariableGenes(object = goblet.cell, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 1,do.plot = T)
length(x = goblet.cell@var.genes)
goblet.cell <- RunPCA(object = goblet.cell,pc.genes = goblet.cell@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)

goblet.cell <- FindClusters(object = goblet.cell, reduction.type = "pca", dims.use = 1:10,
    resolution = c(0.6), print.output = 0, save.SNN = TRUE)
goblet.cell <- RunTSNE(object = goblet.cell, dims.use = 1:10, do.fast = TRUE)
 

## Cancer cells
Cancer.cells <- SubsetData(dt,ident.use = c("Cancer cell"),do.center = T,do.scale = T)
Cancer.cells <- FindVariableGenes(object = Cancer.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 1,do.plot = T)
length(x = Cancer.cells@var.genes)
Cancer.cells <- RunPCA(object = Cancer.cells,pc.genes = Cancer.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
Cancer.cells <- FindClusters(object = Cancer.cells, reduction.type = "pca", dims.use = 1:20,resolution = 1, print.output = 0, save.SNN = TRUE)
Cancer.cells <- RunTSNE(object = Cancer.cells, dims.use = 1:20, do.fast = TRUE)


## Enteroendocrine (EEC) cells 
EEC.cells <- SubsetData(dt,ident.use = c("Enteroendocrine"),do.center = T,do.scale = T)
EEC.cells <- FindVariableGenes(object = EEC.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 3, y.cutoff = 1,do.plot = T)
length(x = EEC.cells@var.genes)
EEC.cells <- RunPCA(object = EEC.cells,pc.genes = EEC.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
EEC.cells <- FindClusters(object = EEC.cells, reduction.type = "pca", dims.use = 1:20,resolution = 1, print.output = 0, save.SNN = TRUE)
EEC.cells <- RunTSNE(object = EEC.cells, dims.use = 1:20, do.fast = TRUE)

## Microenvironment cell lineages
# B cells
B.cells <- SubsetData(dt,ident.use = c("B cell"),do.center = T,do.scale = T)
B.cells <- FindVariableGenes(object = B.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.01, x.high.cutoff = 3, y.cutoff = 1,do.plot = T)
length(x = B.cells@var.genes)
B.cells <- RunPCA(object = B.cells,pc.genes = B.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
B.cells <- FindClusters(object = B.cells, reduction.type = "pca", dims.use = 1:20,resolution = 2, print.output = 0, save.SNN = TRUE)
B.cells <- RunTSNE(object = B.cells, dims.use = 1:20, do.fast = TRUE)


# Fibroblasts
Fibroblast <- SubsetData(dt,ident.use = c("Fibroblast"),do.center = T,do.scale = T)
Fibroblast <- FindVariableGenes(object = Fibroblast, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.01, x.high.cutoff = 3, y.cutoff = 1,do.plot = T)
length(x = Fibroblast@var.genes)
Fibroblast <- RunPCA(object = Fibroblast,pc.genes = Fibroblast@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
Fibroblast <- FindClusters(object = Fibroblast, reduction.type = "pca", dims.use = 1:20,resolution = 1, print.output = 0, save.SNN = TRUE)
Fibroblast <- RunTSNE(object = Fibroblast, dims.use = 1:20, do.fast = TRUE)

# T cells
T.cells <- SubsetData(dt,ident.use = c("T cell"),do.center = T,do.scale = T)
T.cells <- FindVariableGenes(object = T.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.01, x.high.cutoff = 3, y.cutoff = 1,do.plot = T)
length(x = T.cells@var.genes)
T.cells <- RunPCA(object = T.cells,pc.genes = T.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
T.cells <- FindClusters(object = T.cells, reduction.type = "pca", dims.use = 1:20,resolution = 1, print.output = 0, save.SNN = TRUE)
T.cells <- RunTSNE(object = T.cells, dims.use = 1:20, do.fast = TRUE)

```


##### Figure plotting #####
# color code: NAG: #7EA0C2; CAG: #F4B185; IM:#FDB3AC; CAN:#FE4A47

## Figure 1 ## 
## Overall brid's-eye view
```{r}

library(dplyr)
library(reshape2)
library(Seurat)

load("data.temp.res.2.annotation.Rdata")

{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}
# figure 1a

# Imported from predetermined figures

# figure 1b
dt1 <- DBClustDimension(dt,G.use = 0.7)
invalid.cells <- dt1@cell.names[which(dt1@ident == 1)]
valid.cells  <-  setdiff(dt@cell.names,invalid.cells)
dt2 <- SubsetData(dt,cells.use = valid.cells)

pdf(file = "Figure1b.pdf",width = 10,height = 11)
plot.order <- c("Enteroendocrine","GMC","Stem-like cell","Fibroblast","Mast cell","EC","Cancer cell","Chief cell","Macrophage","PC","T cell","Goblet cell","SM cell","Neck cell","Enterocyte","B cell","PMC")
gg <- TSNEPlot(object = dt2,do.label = F,pt.size = 1,do.return = T,plot.order = plot.order) 
gg <- gg + theme(panel.border = element_blank(),
                 axis.line = element_line(),
                 legend.text = element_text(size = 20),
                 legend.position = 'bottom',
                 legend.key.height = unit(0.8,"cm"),
                 legend.key.width = unit(0.6,"cm"),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 20),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20))
gg
dev.off()


# figure 1c
library(dplyr)
all.markers <- FindAllMarkers(dt,logfc.threshold = 1,only.pos = T)
save(all.markers,file = "markers.celltype.Rdata")
load("markers.celltype.Rdata")
all.markers <- read.csv("markers.celltype.csv",stringsAsFactors = F)
top3 <- all.markers %>% group_by(cluster) %>% top_n(3,avg_logFC)
dt.sample <- SubsetData(dt,cells.use = sample(dt@cell.names,6000))
pdf(file = "Figure1c.pdf",height = 16,width = 16)
DoHeatmap(dt.sample,genes.use = top3$gene,slim.col.label = TRUE,remove.key = TRUE,cex.row = 18,rotate.key = T)
dev.off()
rm(dt.sample)

# figure 1d
{
load("markers.celltype.Rdata")
marker.genes <- all.markers$gene

load("pathway.gene.kegg.RData")
NFKB.markers <- valid.final.pathway.element.symbol$hsa04210
NFKB.markers <- intersect(NFKB.markers,marker.genes)
NFKB.markers <- c("NFKB1","LY96","IL6","PTGS2","TNFSF13B","TNFSF10")

cytokine.chemokine.genes <- valid.final.pathway.element.symbol$hsa04060
Chemokine.genes <- valid.final.pathway.element.symbol$hsa04062
Chemokine.genes <- intersect(Chemokine.genes,marker.genes)

cytokine.genes <- read.table("cytokine.txt",sep="\t",stringsAsFactors = F)
cytokine.genes <- cytokine.genes$V1
cytokine.genes <- intersect(cytokine.genes,marker.genes)
cytokine.genes <- c("IL1B","IFNG","TNF","BMP4","TGFBR2")

# dt.temp <- SubsetData(dt,ident.remove = c("PMC","Stem-like cell","Enterocyte","GMC","Enteroendocrine","PC","Neck cell","Cancer cell","Goblet cell","NA"),do.scale = T,do.center = T)

# {
# current.cluster.ids <- c("NA")
# new.cluster.ids <- c("Fibroblast")
# dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
# }

dt@meta.data$groups <- " "

pdf("figure1d-1.pdf",height = 8,width = 5,onefile = F)
SplitDotPlotGG(dt, grouping.var = "groups",NFKB.markers,cols.use = "red", x.lab.rot = T, plot.legend = F, dot.scale = 10,do.return = T) + theme(legend.text = element_text(size = 10),
                 legend.position = 'bottom',
                 legend.key.height = unit(1,"cm"),
                 legend.key.width = unit(1,"cm"),
                 legend.title =element_text(size = 15),
                 axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15),
                 axis.title.x = element_text(size = 0),
                 axis.title.y = element_text(size = 0))
dev.off()

pdf("figure1d-2.pdf",height = 8,width = 4,onefile = F)
SplitDotPlotGG(dt, grouping.var = "groups",cytokine.genes,cols.use = "red", x.lab.rot = T, plot.legend = F, dot.scale = 10,do.return = T) + theme(legend.text = element_text(size = 10),
                 legend.position = 'bottom',
                 legend.key.height = unit(1,"cm"),
                 legend.key.width = unit(1,"cm"),
                 legend.title =element_text(size = 15),
                 axis.text.x = element_text(size = 15),
                 axis.text.y = element_text(size = 15),
                 axis.title.x = element_text(size = 0),
                 axis.title.y = element_text(size = 0))
dev.off()
}

## Focusing on the fatty acid metabolism
library(org.Hs.eg.db)
library(annotate)

metabolism.fatty.acid <- read.table("lipid2genes.txt",sep = "\t",stringsAsFactors = F)
metabolism.fatty.acid.genes <- metabolism.fatty.acid[9,]
metabolism.fatty.acid.genes <- as.character(metabolism.fatty.acid.genes[-which(is.na(metabolism.fatty.acid.genes))][-1])
metabolism.fatty.acid.genes.symbol <- metabolism.fatty.acid.genes
for(i in 1:length(metabolism.fatty.acid.genes.symbol)){
  metabolism.fatty.acid.genes.symbol[i] <- lookUp(metabolism.fatty.acid.genes[i], 'org.Hs.eg', 'SYMBOL')[[1]]
}
metabolism.fatty.acid.genes.symbol <- unique(metabolism.fatty.acid.genes.symbol)

## plot the TSNE-plot for these fatty acid metabolism-related genes

```



## Figure 2 ##
## Focusing the epithelial cells
```{r}
library(dplyr)
library(reshape2)
library(Seurat)

load("data.temp.res.2.annotation.Rdata")
{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}




### 

### Figure 2a
#color code: NAG: #7EA0C2; CAG: #F4B185; IM:#FDB3AC; CAN:#FE4A47

Epithelial.cells <- SubsetData(dt,ident.use = c("PMC","GMC","Stem-like cell","PC","Neck cell","Chief cell","Cancer cell","Enterocyte","Goblet cell","Enteroendocrine"),do.center = T,do.scale = T)

TSNEPlot(object = Epithelial.cells,do.label = T,label.size = 10) 

dt1 <- DBClustDimension(Epithelial.cells,G.use = 0.8)
invalid.cells <- dt1@cell.names[which(dt1@ident == 1)]
valid.cells  <-  setdiff(Epithelial.cells@cell.names,invalid.cells)
dt2 <- SubsetData(Epithelial.cells,cells.use = valid.cells)

{
dt2@meta.data$batch[dt2@meta.data$batch %in% c("NAG1","NAG2")] <- "NAG"
dt2@meta.data$batch[dt2@meta.data$batch %in% c("CAG1","CAG2","CAG3")] <- "CAG"
dt2@meta.data$batch[dt2@meta.data$batch %in% c("IMS1","IMS2","IMS3","IMS4","IMW1","IMW2")] <- "IM"   
dt2@meta.data$batch[dt2@meta.data$batch %in% c("CAN2")] <- "EGC"
  
pdf(file = "Figure2a-2.pdf",width = 8,height = 8)
plot.order <- c("EGC","IM","CAG","NAG")
cols.use <- c("#FE4A47","#FDB3AC","#F4B185","#7EA0C2")
cols.use <- rev(cols.use)
gg <- TSNEPlot(object = dt2,group.by = "batch",do.label = F,pt.size = 1.5,do.return = T,plot.order = plot.order,colors.use = cols.use) 
gg <- gg + theme(panel.border = element_blank(),
                 axis.line = element_line(),
                 legend.text = element_text(size = 25),
                 legend.position = 'bottom',
                 legend.key.height = unit(1.4,"cm"),
                 legend.key.width = unit(1.4,"cm"),
                 axis.text.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 25))
gg
dev.off() 

}  

{
pdf(file = "Figure2a-3.pdf",width = 8,height = 8)
plot.order <- c("Enteroendocrine","GMC","Stem-like cell","Cancer cell","Chief cell","PC","Goblet cell","Neck cell","Enterocyte","PMC")
cols.use <- c("#FF8AB1","#FF74CA","#F166E8","#71C7F5","#00BCD6","#27C597","#45B500","#C6AF61","#D09402","#F8766D")
cols.use <- rev(cols.use)
gg <- TSNEPlot(object = dt2,do.label = F,pt.size = 1.5,do.return = T,plot.order = plot.order,colors.use = cols.use) 
gg <- gg + theme(
                 panel.border = element_blank(),
                 axis.line = element_line(),
                 legend.text = element_text(size = 15),
                 legend.position = 'bottom',
                 legend.key.height = unit(0.5,"cm"),
                 legend.key.width = unit(0.5,"cm"),
                 axis.text.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 25))
gg
dev.off()

}



  
 
### Figure 2b --supplied by YiDing Zhang
{

rename<-function(name,from,to){
  name[which(name %in% from)]<-to
  return(name)
}
lesion<-dt@meta.data$batch
lesion<-rename(lesion,"NAG1","NAG")
lesion<-rename(lesion,c("CAG1","CAG2","CAG3"),"CAG")
lesion<-rename(lesion,c("IMW1","IMW2","IMW3","IMS1","IMS2","IMS3","IMS4"),"IMS")
lesion<-rename(lesion,c("CAN2"),"EGC")
dt@meta.data$lesion<-factor(lesion,levels = c("NAG","CAG","IMS","EGC"))

#####extract epi cell
epi<- SubsetData(dt,ident.use =c("PMC","GMC","Enteroendocrine","Cancer cell","Enterocyte","PC","Stem-like cell" ,"Goblet cell" ,"Chief cell" ,"Neck cell"),do.scale = T,do.center = T)

# TSNEPlot(object = epi,do.label = T,label.size = 5 )
# 
# epi <- FindVariableGenes(object = epi, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 1,do.plot = T)
# 
# length(x = epi@var.genes)
# # cancer.cell <- ScaleData(object = cancer.cell, do.scale = T,do.center = T,vars.to.regress = c("nUMI", "percent.mito"))
# epi<- RunPCA(object =epi, pc.genes = epi@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5,pcs.compute = 120)
# 
# epi <- FindClusters(object = epi, reduction.type = "pca", dims.use = 1:20,
#     resolution = 0.6, print.output = 0, save.SNN = TRUE)
# 
# epi<- RunTSNE(object =epi, dims.use = 1:20, do.fast = TRUE)

# library(RColorBrewer)
# color<-c(brewer.pal(8,"Paired"),brewer.pal(11,"RdGy"))
# barplot(rep(1,19),col=color)
# TSNEPlot(object = epi,group.by="lesion",do.label = T,label.size = 5,colors.use = color[c(3,4,1,2,11)])
# saveRDS(epi,file = "./epi.rds")

##

t<-factor(epi@meta.data$lesion,levels = c("NAG","CAG","IMS","EGC"))
names(t)<-epi@cell.names
epi@ident<-t
epi.nag<-SubsetData(epi,ident.use =c("NAG"),do.scale = T,do.center = T)
epi.cag<-SubsetData(epi,ident.use =c("CAG"),do.scale = T,do.center = T)
epi.ims<-SubsetData(epi,ident.use =c("IMS"),do.scale = T,do.center = T)
epi.can<-SubsetData(epi,ident.use =c("EGC"),do.scale = T,do.center = T)

nag.num <- data.frame(table(epi.nag@meta.data$cell))
cag.num <- data.frame(table(epi.cag@meta.data$cell))
ims.num <- data.frame(table(epi.ims@meta.data$cell))
can.num <- data.frame(table(epi.can@meta.data$cell))

epi.nag.cells <- epi.nag$ident

cellnum<-data.frame(CellType=nag.num$Var1,NAG=nag.num$Freq,CAG=cag.num$Freq,IMS=ims.num$Freq,EGC=can.num$Freq)

for (i in 1:nrow(cellNum)){
  SUM<-sum(cellNum[i,-1])
  for (j in 2:ncol(cellNum)){
    cellNum[i,j]<-cellNum[i,j]/SUM
  }
}
colnames(cellNum)[1]<-"Lesions"
library(reshape2)
cellNum<- melt(cellNum, id.var="Lesions")
cellNum$Lesions<-factor(cellNum$Lesions,levels=c("NAG","CAG","IM","EGC"))
colnames(cellNum)<-c("Lesions","Cell.type","Proportion")
library(RColorBrewer)
#levels(cellNum$Cell.type)<-c("PMC","Neck.cell","GMC","PC","Stem.like.cell","Cancer.cell","Enteroendocrine","Chief.cell" ,"Enterocyte" ,"Goblet.cell")
cellNum$Cell.type<-factor(cellNum$Cell.type,levels =c("PMC","Neck.cell","GMC","PC","Enteroendocrine","Chief.cell" ,"Enterocyte" ,"Goblet.cell","Stem.like.cell","Cancer.cell") )
#c("#F8766D","#C6AF61","#FF74CA","#27C597","#F166E8","#71C7F5","#FF8AB1","#00BCD6","#D09402","#45B500")
colors=c("#F8766D","#C6AF61","#FF74CA","#27C597","#FF8AB1","#00BCD6","#D09402","#45B500","#F166E8","#71C7F5")
library(ggplot2)
pdf("figure2b.pdf",width = 8,height = 8,onefile = F)
ggplot(cellNum, aes(x =Lesions, y = Proportion, fill = Cell.type)) + 
  geom_bar(stat = "identity")+scale_fill_manual(values=colors )+theme(axis.text=element_text(size=25), axis.title=element_text(size=25),legend.text = element_text(size=25),legend.title=element_text(size=25),legend.key.size =  unit(1, "in"))
dev.off()

}

### Supplementary_Figure 2b
{
library(pheatmap)
load("data.temp.res.2.annotation.Rdata")
{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}

Epithelial.cells <- SubsetData(dt,ident.use = c("PMC","GMC","Stem-like cell","PC","Neck cell","Chief cell","Cancer cell","Enterocyte","Goblet cell","Enteroendocrine"),do.center = T,do.scale = T) 
Epithelial.cells@meta.data$batch[Epithelial.cells@meta.data$batch %in% c("NAG1","NAG2")] <- "NAG"
Epithelial.cells@meta.data$batch[Epithelial.cells@meta.data$batch %in% c("CAG1","CAG2","CAG3")] <- "CAG"
Epithelial.cells@meta.data$batch[Epithelial.cells@meta.data$batch %in% c("IMS1","IMS2","IMS3","IMS4","IMW1","IMW2")] <- "IM"   
Epithelial.cells@meta.data$batch[Epithelial.cells@meta.data$batch %in% c("EGC")] <- "EGC"
Epithelial.cells@meta.data$orig.ident <- Epithelial.cells@ident
Epithelial.cells@meta.data$orig.ident <- factor(paste(Epithelial.cells@meta.data$batch,Epithelial.cells@ident))
Epithelial.cells <- SetAllIdent(Epithelial.cells,"orig.ident")


all.cells <- as.numeric(table(Epithelial.cells@meta.data$batch))
Epithelial.cells.<- SubsetData(Epithelial.cells,cells.use = Epithelial.cells@cell.names[grep("PC",Epithelial.cells@ident)]) 
PC.cells <- as.numeric(table(Epithelial.cells.PC@meta.data$batch))
rm(Epithelial.cells.PC)

Epithelial.cells.Stem <- SubsetData(Epithelial.cells,cells.use = Epithelial.cells@cell.names[grep("Stem",Epithelial.cells@ident)]) 
Stem.cells <- as.numeric(table(Epithelial.cells.Stem@meta.data$batch))
rm(Epithelial.cells.Stem)  

pdf("Supplementary_Figure2b-Stem.pdf",width = 6,height = 6)
propor.cells <- as.vector(Stem.cells/all.cells)
names(propor.cells) <- names(table(Epithelial.cells@meta.data$batch))
names(propor.cells)[2] <- "EGC"
propor.cells <- propor.cells[c("NAG","CAG","IM","EGC")]
barplot(propor.cells,ylab = "The proportion of stem-like cells")
dev.off() 
}


# Epithelial.cells.markers <- FindAllMarkers(Epithelial.cells,only.pos = T,min.cells.group = 80)
epithelial.cells.markers <- read.csv("Epithelial.cells.marker.csv",stringsAsFactors = F)

### Figure 2c
{
valid.cells <- Epithelial.cells@cell.names[which(Epithelial.cells@ident %in% names(table(Epithelial.cells@ident))[which(as.numeric(table(Epithelial.cells@ident))>80)])]
# sampling data
valid.cells <- sample(valid.cells, round(length(valid.cells)/5))

epithelial.heatmap <- SubsetData(Epithelial.cells,cells.use = valid.cells,do.center = T,do.scale = T)


# predetermine the groups 
pheatmap.data <- as.matrix(epithelial.heatmap@scale.data)
cells.order <- order(ordered(epithelial.heatmap@ident, levels = c("NAG PMC","NAG PC","NAG Neck cell","NAG Enteroendocrine","CAG PMC","CAG Neck cell","CAG GMC","CAG PC","CAG Enteroendocrine","IM GMC","IM PMC","IM PC","IM Stem-like cell","IM Enteroendocrine","IM Goblet cell","IM Enterocyte","EGC PMC","EGC PC","EGC Goblet cell","EGC Enteroendocrine","EGC Stem-like cell","EGC Cancer cell")))
cells.names.order <- epithelial.heatmap@cell.names[cells.order]

top30 <- epithelial.cells.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
gene.order <- order(ordered(top30$cluster, levels = c("NAG PMC","NAG PC","NAG Neck cell","NAG Enteroendocrine","CAG PMC","CAG Neck cell","CAG GMC","CAG PC","CAG Enteroendocrine","IM GMC","IM PMC","IM PC","IM Stem-like cell","IM Enteroendocrine","IM Goblet cell","IM Enterocyte","CAN PMC","CAN PC","CAN Goblet cell","CAN Enteroendocrine","CAN Stem-like cell","CAN Cancer cell")))

top30 <- top30[gene.order,]
genes.names.order <- as.character(top30$gene)

pheatmap.data <- as.matrix(pheatmap.data)[genes.names.order,cells.names.order]
pheatmap.data <- t(as.data.frame(pheatmap.data))

annotation_col = data.frame(CellType = factor(sapply(strsplit(as.character(epithelial.heatmap@ident[cells.order]),split = " "),function(x){paste(x[-1],collapse =  " ")})), Stage = factor(epithelial.heatmap@meta.data$batch[cells.order]))
rownames(annotation_col) <- rownames(pheatmap.data)

col_temp <- data.frame(apply(annotation_col,1,function(x){x = x[1];return(paste(strsplit(x, split = "-")[[1]],collapse = "."))}))
col_temp <- data.frame(apply(col_temp,1,function(x){x = x[1];return(paste(strsplit(x, split = " ")[[1]],collapse = "."))}))
annotation_col$CellType <- col_temp[,1]

# Setting gaps
gap_col <- as.numeric(table(annotation_col$Stage)[order(ordered(names(table(annotation_col$Stage)), levels = c("NAG","CAG","IM","EGC")))])
gap_col[1] <- gap_col[1]
gap_col[2] <- sum(gap_col[c(1,2)])
gap_col[3] <- sum(gap_col[c(2,3)])
gap_col[4] <- sum(gap_col[c(3,4)])
gap_row <- c(80,180,320,440)

# Setting colors
#color code: NAG: #7EA0C2; CAG: #F4B185; IM:#FDB3AC; CAN:#FE4A47
# plot.order <-  <- c("Enteroendocrine","GMC","Stem-like cell","Cancer cell","Chief cell","PC","Goblet cell","Neck cell","Enterocyte","PMC")
# ")
# cols.use <-e <- c("#FF8AB1","#FF74CA","#F166E8","#71C7F5","#00BCD6","#27C597","#45B500","#C6AF61","#D09402","#F8766D")

ann_colors = list(Stage = c(NAG = "#7EA0C2", CAG = "#F4B185",IM = "#FDB3AC",CAN = "#FE4A47"),
                  CellType = c(PMC = "#F8766D",GMC = "#FF74CA",Enterocyte = "#D09402",= "#27C597",Neck.cell = "#C6AF61",Enteroendocrine = "#FF8AB1",Goblet.cell = "#45B500",Stem.like.cell = "#F166E8",Cancer.cell = "#71C7F5"))

# plot the heatmap
# pheatmap(pheatmap.data, annotation_col = annotation_col,annotation_colors = ann_colors,cluster_rows = F,cluster_col = F,gaps_col = gap_col)
# pdf("figure2d-1.pdf",width = 20,height = 20)
# pheatmap(pheatmap.data, annotation_col = annotation_col,cluster_rows = F,cluster_col = F)
# dev.off()
breaks = c(seq(-3.5, 0, length.out = 101),seq(0.01, 4, length.out = 101))
breaks2 <- breaks

#setting colors
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 201)

pdf("figure2c.pdf",width = 20,height = 20,onefile = F)
pheatmap(t(pheatmap.data), cluster_rows = F,cluster_col = F,show_rownames = F,show_colnames = F,gaps_col = gap_col,gaps_row = c(),annotation_col = annotation_col,breaks = breaks2, color = my_palette,fontsize = 25,annotation_colors = ann_colors)
dev.off()

tiff("figure2d.tiff",width = 200,height = 200,res = 50,units = "mm")
pheatmap(t(pheatmap.data), cluster_rows = F,cluster_col = F,show_rownames = F,show_colnames = F,gaps_col = gap_col,gaps_row = c(),annotation_col = annotation_col,breaks = breaks2, color = my_palette,fontsize = 25,annotation_colors = ann_colors)
dev.off()

pdf("figure2c-xx.pdf",width = 20,height = 20,onefile = F)
pheatmap(t(pheatmap.data)[,data.temp], cluster_rows = F,cluster_col = F,show_rownames = F,show_colnames = F,annotation_col = annotation_col[data.temp,],breaks = breaks2, color = my_palette,fontsize = 25,annotation_colors = ann_colors)
dev.off()
}



##Figure 2d
{
{
library(GEOquery)
library(annotate)
library(org.Hs.eg.db) 
  
GSE2669.rawdata <- getGEO("GSE2669")
GSE2669.rawdata <- GSE2669.rawdata$GSE2669_series_matrix.txt.gz

GSE2669.symbol <- data.frame(matrix(nrow =  as.numeric(nrow(GSE2669.rawdata@featureData)),ncol = 1))

for(i in 1:nrow(GSE2669.symbol)){
  print(i)
  gene_i <- GSE2669.rawdata@featureData@data$GENE[i]
  GSE2669.symbol[i,1] <- lookUp(gene_i, 'org.Hs.eg', 'SYMBOL')[[1]]
}
colnames(GSE2669.symbol)  <- "Gene"

GSE2669.expressiondata <- GSE2669.rawdata@assayData$exprs

GSE2669.phenodata <-   GSE2669.rawdata@phenoData@data
GSE2669.phenodata.data <- rownames(GSE2669.phenodata)[which(GSE2669.phenodata$description %in% c("Biopsy; Laurens::INTESTINAL","Biopsy; Chronic Gastritis","Biopsy; Intestinal Metaplasia","Biopsy; Normal Stomach"))]
GSE2669.phenodata <- data.frame(GSE2669.phenodata[GSE2669.phenodata.data,"description"])
rownames(GSE2669.phenodata) <- as.character(GSE2669.phenodata.data)

invalid.index1 <- which(is.na(GSE2669.symbol$Gene))
invalid.index2 <- which(duplicated(GSE2669.symbol$Gene))
invalid.index <- union(invalid.index1,invalid.index2)
GSE2669.symbol <- GSE2669.symbol[-invalid.index,]

GSE2669.expressiondata <- GSE2669.expressiondata[-invalid.index,rownames(GSE2669.phenodata)]
rownames(GSE2669.expressiondata) <- GSE2669.symbol

result.all <- list(expression = GSE2669.expressiondata,phenotype = GSE2669.phenodata)
}

save(result.all,file = "GSE2669.Rdata")


# annova test
load("C:/Users/Peng-Zhang/Desktop/胃炎癌-多部位研究/Results/Rdata/GSE2669.Rdata")

GSE2669.expressiondata <- result.all$expression
GSE2669.phenodata <- result.all$phenotype


GSE2669.expressiondata <- data.frame(t(GSE2669.expressiondata))

GSE2669.expressiondata$group <- as.vector(as.character(GSE2669.phenodata$GSE2669.phenodata.GSE2669.phenodata.data...description..))


chg.index <- grep("Gastritis",GSE2669.expressiondata$group);non.chg.index <- grep("Normal Stomach",GSE2669.expressiondata$group)
chg.p <- data.frame(apply(GSE2669.expressiondata[,1:(ncol(GSE2669.expressiondata)-1)],2,function(x){t.test(x[chg.index],x[-chg.index])$p.val}))
chg.q <- data.frame(p.adjust(chg.p$apply.GSE2669.expressiondata...1..ncol.GSE2669.expressiondata....));rownames(chg.q) <- rownames(chg.p)
chg.fold.change <- data.frame(apply(GSE2669.expressiondata[,1:(ncol(GSE2669.expressiondata)-1)],2,function(x){mean(x[chg.index]/mean(x[non.chg.index]))}))
chg.result <- cbind(chg.q,chg.fold.change);colnames(chg.result) <- c("p.val","fold.change")
chg.deg <- rownames(chg.result)[which((chg.result$p.val<0.01)|(chg.result$fold.change>1.5))]


IM.index <- grep("Intestinal Metaplasia",GSE2669.expressiondata$group);non.IM.index <- setdiff(1:nrow(GSE2669.phenodata),IM.index)
IM.p <- data.frame(apply(GSE2669.expressiondata[,1:(ncol(GSE2669.expressiondata)-1)],2,function(x){t.test(x[IM.index],x[non.IM.index])$p.val}))
IM.q <- data.frame(p.adjust(IM.p$apply.GSE2669.expressiondata...1..ncol.GSE2669.expressiondata....));rownames(IM.q) <- rownames(IM.p)
IM.fold.change <- data.frame(apply(GSE2669.expressiondata[,1:(ncol(GSE2669.expressiondata)-1)],2,function(x){mean(x[IM.index]/mean(x[non.IM.index]))}))
IM.result <- cbind(IM.q,IM.fold.change);colnames(IM.result) <- c("p.val","fold.change")
IM.deg <- rownames(IM.result)[which((IM.result$p.val<0.01)|(IM.result$fold.change>1.5))]

CAN.index <- grep("Laurens::INTESTINAL",GSE2669.expressiondata$group);non.CAN.index <-setdiff(1:nrow(GSE2669.phenodata),CAN.index)
CAN.p <- data.frame(apply(GSE2669.expressiondata[,1:(ncol(GSE2669.expressiondata)-1)],2,function(x){t.test(x[CAN.index],x[non.CAN.index])$p.val}))
CAN.q <- data.frame(p.adjust(CAN.p$apply.GSE2669.expressiondata...1..ncol.GSE2669.expressiondata....));rownames(CAN.q) <- rownames(CAN.p)
CAN.fold.change <- data.frame(apply(GSE2669.expressiondata[,1:(ncol(GSE2669.expressiondata)-1)],2,function(x){mean(x[CAN.index]/mean(x[non.CAN.index]))}))
CAN.result <- cbind(CAN.q,CAN.fold.change);colnames(CAN.result) <- c("p.val","fold.change")
CAN.deg <- rownames(CAN.result)[which((CAN.result$p.val<0.01)|(CAN.result$fold.change>1.5))]




# comparison with single-cell dataset
epithelial.cells.markers <- read.csv("C:/Users/Peng-Zhang/Desktop/胃炎癌-多部位研究/Cell_Reports_revision_v0320/返修最终提交版/Figures/Article/Figure2/Epithelial.cells.marker.csv",stringsAsFactors = F)
epithelial.cells.markers <- epithelial.cells.markers[which(epithelial.cells.markers$avg_logFC > log2(1.5)),]
epithelial.cells.markers.chg <- epithelial.cells.markers[grep("CAG",epithelial.cells.markers$cluster),]
sc.chg.deg <- epithelial.cells.markers.chg$gene
# sc.chg.deg <- intersect(epithelial.cells.markers.chg$gene,rownames(chg.p))

epithelial.cells.markers.im <- epithelial.cells.markers[grep("IM",epithelial.cells.markers$cluster),]
sc.im.deg <- epithelial.cells.markers.im$gene
# sc.im.deg <- intersect(epithelial.cells.markers.im$gene,rownames(chg.p))

epithelial.cells.markers.can <- epithelial.cells.markers[grep("CAN Cancer cell",epithelial.cells.markers$cluster),]
sc.can.deg <- epithelial.cells.markers.can$gene
# sc.can.deg <- intersect(epithelial.cells.markers.can$gene,rownames(chg.p))

# cipher_prediction
gc.prediction <- read.delim("C:/Users/Peng-Zhang/Desktop/胃炎癌-多部位研究/Cell_Reports_revision_v0320/返修最终提交版/Figures/Article/Figure2/Figure2e/cipher_gastric_cancer.txt",sep = "\t",stringsAsFactors = F,header = T)
gc.prediction <- unique(gc.prediction$Gene) #[-which(is.na(gc.prediction$Score))])
# gc.prediction <- intersect(gc.prediction,rownames(chg.p))

gs.prediction <- read.delim("C:/Users/Peng-Zhang/Desktop/胃炎癌-多部位研究/Cell_Reports_revision_v0320/返修最终提交版/Figures/Article/Figure2/Figure2e/cipher_gastritis.txt",sep = "\t",stringsAsFactors = F,header = F)
gs.prediction <- unique(gs.prediction$V1) #[-which(is.na(gs.prediction$Score))])
# gs.prediction <- intersect(gs.prediction,rownames(chg.p))
cipher.prediction <- union(gc.prediction,gs.prediction)
epithelial.cells.markers.genes <- as.character(epithelial.cells.markers$gene)

## Boxplot for the distributon of CIPHER-inferred high-risk genes in Bulk RNA-seq Datasets and Single-cell datasets

load("E:/SC_rawdata_20190313/data.temp.res.2.annotation.Rdata")
epithelial.cells <- c("PMC","GMC","Enteroendocrine","PC","Neck cell","Chief cell","Cancer cell")

library(dplyr)
library(reshape2)
library(Seurat)

candidate.genes <- intersect(intersect(gc.prediction,sc.can.deg),CAN.deg)
## Bulk datasets
data.temp <- GSE2669.expressiondata[,c(candidate.genes,"group")]

gg.temp <- VlnPlot(dt,features.plot = "GSTM4",ident.include = epithelial.cells,do.return = T,point.size.use = -1, x.lab.rot = 45)

gg.temp <- gg.temp + theme(panel.border = element_blank(),
                 axis.line = element_line(),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 20),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20))

gg.temp








## Venndiagram 
library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    A=sc.can.deg,
    B=CAN.deg,
    C = gc.prediction
  ),
  filename = "figure2d-EGC.tiff",
  col = "transparent",
  category.names=c("SC","GSE2669","Prediction"),
  fill = c("cornflowerblue", "darkorchid1","yellow"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col =  c("dodgerblue", "seagreen3","darkorange1"),
  cat.cex = 1.5,
  cat.dist=c(0.1, 0.1, 0.05),
  cat.fontfamily = "serif",
  rotation.degree = 0,
  margin = 0.2,
  reverse=FALSE
);

intersect(intersect(sc.can.deg,CAN.deg),gc.prediction)



library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    A= sc.im.deg,
    B= IM.deg,
    C = gs.prediction
  ),
  filename = "figure2d-IM.tiff",
  col = "transparent",
  category.names=c("IM","GSE2669","CIPHER"),
  fill = c("cornflowerblue", "darkorchid1","yellow"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col =  c("dodgerblue", "seagreen3","darkorange1"),
  cat.cex = 1.5,
  cat.dist=c(0.1, 0.1, 0.1),
  cat.fontfamily = "serif",
  rotation.degree = 0,
  margin = 0.2,
  reverse=FALSE
);

sc.chg.deg <- union(sc.chg.deg,sc.im.deg)
chg.deg <- union(chg.deg,IM.deg)

library(VennDiagram)
venn.plot <- venn.diagram(
  x = list(
    A= sc.chg.deg,
    B= chg.deg,
    C = gs.prediction
  ),
  filename = "figure2d-CAG.tiff",
  col = "transparent",
  category.names=c("SC","GSE2669","Prediction"),
  fill = c("cornflowerblue", "darkorchid1","yellow"),
  alpha = 0.50,
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col =  c("dodgerblue", "seagreen3","darkorange1"),
  cat.cex = 1.5,
  cat.dist=c(0.1, 0.1, 0.05),
  cat.fontfamily = "serif",
  rotation.degree = 0,
  margin = 0.2,
  reverse=FALSE
);

intersect(intersect(sc.chg.deg,chg.deg),gs.prediction)


}


### Figure 2f
# Network --unique genes # Outputs -- node info & edges
# node info including celltypes, high-risk genes? receptor or ligands or 
{
# import netwrok data
load("ppi.string.symbol.RData")
network <- ppi_string_symbol[which(ppi_string_symbol$combined_score>700),c(1,2)]
network <- network[-which(is.na(network$protein1) | is.na(network$protein2)),]
network <- t(data.frame(apply(network,1,function(x){as.vector(sort(x))})))
network.temp <- apply(network,1,function(x){paste(x,collapse = "_")})
network <- data.frame(network[-which(duplicated(network.temp)),])
colnames(network) <- c("protein1","protein2")

epithelial.cells.markers <- Epithelial.cells.markers
epithelial.cells.markers <- epithelial.cells.markers[-which(as.character(epithelial.cells.markers$gene) %in% "UBC"),]
epithelial.cells.markers <- epithelial.cells.markers[which(epithelial.cells.markers$avg_logFC > log2(1.5)),]
epithelial.cells.markers.uniq <- names((table(epithelial.cells.markers$gene)))[which(table(epithelial.cells.markers$gene)<3)]
epithelial.cells.markers <- epithelial.cells.markers[which(epithelial.cells.markers$gene %in% epithelial.cells.markers.uniq),]

# epithelial.cells.markers <- epithelial.cells.markers[-which(is.na(epithelial.cells.markers$gene)),]

# rownames(epithelial.cells.markers) <- epithelial.cells.markers$gene

# functions
network_analysis <- function(network,key_words,epithelial.cells.markers){
  
 epithelial.cells.markers <- epithelial.cells.markers[which(epithelial.cells.markers$cluster %in% key_words),]
 epithelial.cells.markers.genes <- epithelial.cells.markers$gene 
   
 network.temp <- network[which((as.character(network$protein1) %in% epithelial.cells.markers.genes)&(as.character(network$protein2) %in% epithelial.cells.markers.genes)),]
 network.node <- network.temp
 
 valid.genes <- union(network.temp$protein1,network.temp$protein2)
 node.info <- epithelial.cells.markers[which(epithelial.cells.markers$gene %in% valid.genes),]
 node.info$celltype <-  sapply(strsplit(as.character(node.info$cluster),split = " "),function(x){paste(x[-1],collapse =  " ")})
 node.info <- node.info[-which(duplicated(node.info$gene)),]
 rownames(node.info) <- node.info$gene
 # 
 # node.info[,c(1,7)] <-  node.info[,c(7,1)]
 # colnames(node.info)[c(1,7)] <- colnames(node.info)[c(7,1)]
 
 node.info$cipher <- 0
 node.info$cipher[which(node.info$gene %in% cipher.prediction)] <- 1
 
 final.results <- list(network.node,node.info)  
 
 return(final.results)
 
}

epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.1 <- epithelial.cells.markers.temp[grep("CAG GMC",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.1 <- epithelial.cells.markers.temp.1$gene 

epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.2 <- epithelial.cells.markers.temp[grep("IM GMC",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.2 <- epithelial.cells.markers.temp.2$gene 

epithelial.cells.markers.genes.1 <- union(epithelial.cells.markers.genes.1,epithelial.cells.markers.genes.2)


epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.3 <- epithelial.cells.markers.temp[grep("CAG PMC",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.3 <- epithelial.cells.markers.temp.3$gene 

epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.4 <- epithelial.cells.markers.temp[grep("IM PMC",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.4 <- epithelial.cells.markers.temp.4$gene 

epithelial.cells.markers.genes.4  <- union(epithelial.cells.markers.genes.3,epithelial.cells.markers.genes.4)



epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.5 <- epithelial.cells.markers.temp[grep("IM Goblet cell",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.5 <- epithelial.cells.markers.temp.5$gene 


epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.5 <- epithelial.cells.markers.temp[grep("IM Enterocyte",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.6 <- epithelial.cells.markers.temp.5$gene 



epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.5 <- epithelial.cells.markers.temp[grep("Stem-like cell",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.7 <- epithelial.cells.markers.temp.5$gene 


epithelial.cells.markers.temp <- epithelial.cells.markers
epithelial.cells.markers.temp.5 <- epithelial.cells.markers.temp[grep("CAN Cancer cell",epithelial.cells.markers.temp$cluster),]
epithelial.cells.markers.genes.8 <- epithelial.cells.markers.temp.5$gene 


epithelial.cells.markers.genes <- Reduce(union,c(epithelial.cells.markers.genes.1,epithelial.cells.markers.genes.4, epithelial.cells.markers.genes.5,epithelial.cells.markers.genes.6,epithelial.cells.markers.genes.7,epithelial.cells.markers.genes.8))


 network.temp <- network[which((as.character(network$protein1) %in% epithelial.cells.markers.genes)&(as.character(network$protein2) %in% epithelial.cells.markers.genes)),]
 network.node <- network.temp
 
 valid.genes <- union(network.temp$protein1,network.temp$protein2)
 node.info <- epithelial.cells.markers[which(epithelial.cells.markers$gene %in% valid.genes),]
 node.info$celltype <-  sapply(strsplit(as.character(node.info$cluster),split = " "),function(x){paste(x[-1],collapse =  " ")})
 node.info <- node.info[-which(duplicated(node.info$gene)),]
 rownames(node.info) <- node.info$gene

 
 # node.info[,c(1,7)] <-  node.info[,c(7,1)]
 # colnames(node.info)[c(1,7)] <- colnames(node.info)[c(7,1)]
 
 node.info$cipher <- 0
 node.info$cipher[which(node.info$gene %in% cipher.prediction)] <- 1
 
 
 node.info$GMC <- "GMC"
 node.info$GMC[which(node.info$gene %in% epithelial.cells.markers.genes.4)] <- "PMC"
 node.info$GMC[which(node.info$gene %in% epithelial.cells.markers.genes.5)] <- "Goblet"
 node.info$GMC[which(node.info$gene %in% epithelial.cells.markers.genes.6)] <- "Enterocyte"
 node.info$GMC[which(node.info$gene %in% epithelial.cells.markers.genes.7)] <- "MSC"
 node.info$GMC[which(node.info$gene %in% epithelial.cells.markers.genes.8)] <- "EGC"
 
 node.info <- node.info[,c("GMC","cipher")]
 
 write.table(network.temp,file = "PMC.network.edge.txt",sep = "\t",quote = F,row.names = F)
 write.table(node.info,file = "PMC.network.nodeinfo.txt",sep = "\t",quote = F,row.names = F)

 
 
 
 
# 
#  
# NAG.network  <- list(network.node,node.info)  
# NAG.network.edge <- NAG.network[[1]]
# NAG.network.nodeinfo <- NAG.network[[2]]

# 
# CAG.network <- network_analysis(network,"CAG",epithelial.cells.markers)
# CAG.network.edge <- CAG.network[[1]]
# CAG.network.nodeinfo <- CAG.network[[2]]
# write.table(CAG.network.edge,file = "CAG.network.edge.txt",sep = "\t",quote = F,row.names = F)
# write.table(CAG.network.nodeinfo,file = "CAG.network.nodeinfo.txt",sep = "\t",quote = F,row.names = F)
# 
# IM.network <- network_analysis(network,"IM",epithelial.cells.markers)
# IM.network.edge <- IM.network[[1]]
# IM.network.nodeinfo <- IM.network[[2]]
# write.table(IM.network.edge,file = "IM.network.edge.txt",sep = "\t",quote = F,row.names = F)
# write.table(IM.network.nodeinfo,file = "IM.network.nodeinfo.txt",sep = "\t",quote = F,row.names = F)
# 
# CAN.network <- network_analysis(network,"CAN",epithelial.cells.markers)
# CAN.network.edge <- CAN.network[[1]]
# CAN.network.nodeinfo <- CAN.network[[2]]
# write.table(CAN.network.edge,file = "CAN.network.edge.txt",sep = "\t",quote = F,row.names = F)
# write.table(CAN.network.nodeinfo,file = "CAN.network.nodeinfo.txt",sep = "\t",quote = F,row.names = F)
# 
# 
# NAG.CAG.network <- network_analysis(network,"AG",epithelial.cells.markers)
# NAG.CAG.network.edge <- NAG.CAG.network[[1]]
# NAG.CAG.nodeinfo <- NAG.CAG.network[[2]]
# write.table(NAG.CAG.network.edge,file = "NAG.CAG.network.edge.txt",sep = "\t",quote = F,row.names = F)
# write.table(NAG.CAG.nodeinfo,file = "NAG.CAG.network.nodeinfo.txt",sep = "\t",quote = F,row.names = F)


}


## Figure 2g
all.epithelial.genes <- Reduce(union,c(NAG.network.nodeinfo$gene,CAG.network.nodeinfo$gene,IM.network.nodeinfo$gene,CAN.network.nodeinfo$gene))

# microenvironment cell markers 
load("microenvironment.Rdata")
microenvironment_marker_genelist <- microenvironment_marker$Gene

Epithelial.cells.markers <- read.csv("Epithelial.cells.marker.csv",stringsAsFactors = F)

# ligand-receptor database
ligand.receptor <- read.delim("ligand-receptors.txt",sep = "\t",stringsAsFactors = F)
ligand.receptor.temp <- ligand.receptor[which(((ligand.receptor$Ligand.ApprovedSymbol %in% all.epithelial.genes) & (ligand.receptor$Receptor.ApprovedSymbol %in% microenvironment_marker_genelist))| ((ligand.receptor$Receptor.ApprovedSymbol %in% all.epithelial.genes) & (ligand.receptor$Ligand.ApprovedSymbol %in% microenvironment_marker_genelist))),]

ligand.receptor.temp <- ligand.receptor.temp[,c("Ligand.ApprovedSymbol","Receptor.ApprovedSymbol")]

Epithelial.cells.markers <- Epithelial.cells.markers[-which(duplicated(Epithelial.cells.markers$gene)),]
rownames(Epithelial.cells.markers) <- Epithelial.cells.markers$gene

microenvironment_marker <- microenvironment_marker[-which(duplicated(microenvironment_marker$Gene)),]
rownames(microenvironment_marker) <- microenvironment_marker$Gene

ligand.receptor.temp$ligand.celltype.epi <- Epithelial.cells.markers[ligand.receptor.temp$Ligand.ApprovedSymbol,"cluster"]
ligand.receptor.temp$ligand.celltype.nonepi <- microenvironment_marker[ligand.receptor.temp$Ligand.ApprovedSymbol,"cell.type"]

ligand.receptor.temp$receptor.celltype.epi <- Epithelial.cells.markers[ligand.receptor.temp$Receptor.ApprovedSymbol,"cluster"]
ligand.receptor.temp$receptor.celltype.nonepi <- microenvironment_marker[ligand.receptor.temp$Receptor.ApprovedSymbol,"cell.type"]

write.table(ligand.receptor.temp,file = "ligand.receptor.subset.txt",sep = "\t",quote = F,row.names = F)



```



## Figure3 ##
## Focusing the mucous-secreting cells
```{r considering the Cancer cells}

library(dplyr)
library(reshape2)
library(Seurat)

load("data.temp.res.2.annotation.Rdata")

{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}


## Figure 5a -- Highlight the cancer cell cluster

{
Mucous.cells <- SubsetData(dt,ident.use = c("PMC","GMC"),do.center = T,do.scale = T)
Mucous.cells.list <- Mucous.cells@cell.names
CAG.IM.cells.list <- dt@cell.names[-which(dt@meta.data$batch %in% c("CAN2"))]
valid.cells <- intersect(Mucous.cells.list,CAG.IM.cells.list)

pdf(file = "Figure3a.pdf", width = 8,height = 8)
data.raw <- data.frame(dt@dr$tsne@cell.embeddings)
data.raw$Mucous.cells <- 1
data.raw$Mucous.cells[which(dt@cell.names %in% valid.cells)] <- 0
ggplot(data.raw, aes(x=tSNE_1, y=tSNE_2, color=factor(Mucous.cells))) + geom_point() + scale_colour_manual(values = c("#F8766D","#C4C4C4")) + theme(axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25),axis.text.x = element_text(size = 25),axis.text.y = element_text(size = 25))
dev.off()
}

## Figure 3b-C -- MUC6 & MUC5AC as markers
{
pdf(file = "Figure3b.pdf", width = 5, height =5)
gg <- FeaturePlot(object = dt, features.plot = c("MUC6"), cols.use = c("grey","#F8766D"), reduction.use = "tsne",pt.size = 1,do.return = T,no.legend = F) + theme(legend.text = element_text(size = 40),
                 legend.position = 'bottom',
                 legend.key.height = unit(1.4,"cm"),
                 legend.key.width = unit(1.4,"cm"),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 20))
gg 
dev.off()

pdf(file = "Figure3c.pdf", width = 5, height = 5)
gg <- FeaturePlot(object = dt, features.plot = c("MUC5AC"), cols.use = c("grey","#F8766D"), reduction.use = "tsne",pt.size = 1,do.return = T,no.legend = F) + theme(legend.text = element_text(size = 40),
                 legend.position = 'bottom',
                 legend.key.height = unit(1.4,"cm"),
                 legend.key.width = unit(1.4,"cm"),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 20))
gg 
dev.off()

}

## Figure 3d
{
 PMC.GMC.marker <- FindMarkers(dt,ident.1 = c("GMC"),ident.2 = c("PMC"),logfc.threshold = log2(1.5)) 
 write.table(PMC.GMC.marker,file = "PMC.GMC.marker.csv")
 # Figure 3D was used in the file saved as "C:\Users\Peng-Zhang\Desktop\胃炎癌-多部位研究\论文终版\图&论文\Figures\Final_result_figure_v1003\Figure3/Figure3d.pdf"
 
}


# Supplementary_Figure 3e
{
gastric.mucous.cells <- SubsetData(dt,ident.use = c("PMC","GMC"),do.center = T,do.scale = T)
gastric.mucous.cells  <- FindVariableGenes(object = gastric.mucous.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = gastric.mucous.cells@var.genes)
gastric.mucous.cells  <- RunPCA(object =gastric.mucous.cells,pc.genes = gastric.mucous.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
gastric.mucous.cells <- FindClusters(object = gastric.mucous.cells, reduction.type = "pca", dims.use = 1:20,
    resolution = c(1), print.output = 0, save.SNN = TRUE)
gastric.mucous.cells  <- RunTSNE(object = gastric.mucous.cells, dims.use = 1:20, do.fast = TRUE)

current.cluster.ids <- c(8,1,3,5,14,4,9,0,2,10,7,12,11,6,13)
# new.cluster.ids <- c(rep("PMC_NAG",2),rep("PMC_CAG",7),rep("PMC_IM",2),rep("PMC_EGC",1),rep("GMC_IM_1",1),rep("GMC_CAG",1),rep("GMC_IM_2",1))
new.cluster.ids <- c(rep("NAG",2),rep("CAG",7),rep("IM",2),rep("EGC",1),rep("IM",1),rep("CAG",1),rep("IM",1))
gastric.mucous.cells@ident <- plyr::mapvalues(x = gastric.mucous.cells@ident, from = current.cluster.ids, to = new.cluster.ids)

Pit.cells <- SubsetData(gastric.mucous.cells,ident.use = c("NAG","CAG","IM","EGC"),do.scale = T,do.center = T)
Pit.cells.markers <- FindAllMarkers(Pit.cells,logfc.threshold = log2(1.5),only.pos = T)
Pit.cells.order <- order(match(Pit.cells.markers$cluster, c("NAG","CAG","IM","EGC")))
Pit.cells.markers <- Pit.cells.markers[Pit.cells.order,]
Pit.cells.markers <- Pit.cells.markers %>% group_by(cluster) %>% top_n(20,avg_logFC)

pdf("Supplementary_Figure3e.pdf",width = 18,height = 20)
DoHeatmap(Pit.cells,genes.use = as.character(Pit.cells.markers$gene),slim.col.label = TRUE,remove.key = F,group.order = c("NAG","CAG","IM","EGC"),group.spacing = 0.4,cex.row = 20,cex.col = 40,col.low = "blue",col.mid = "white",col.high = "red",group.cex = 30,group.label.rot = T)
dev.off()
}


## Figure 3e
genes <- c("GAST","PGC","FABP1","TFF3","REG4","CDH17")
pdf(file = "Figure3e.pdf",width = 6,height = 10)
gg <- VlnPlot(object = gastric.mucous.cells, features.plot = genes, remove.legend = T,x.lab.rot = 45, point.size.use = 0,do.return = T,do.sort = F,size.x.use = 25,size.title.use = 30,size.y.use = 20,nCol = 2,same.y.lims = T) 
gg
dev.off()



## Figure 3f-h
# Figure 3f
{
gastric.mucous.cells <- SubsetData(dt,ident.use = c("GMC"),do.center = T,do.scale = T)
gastric.mucous.cells@meta.data$batch[gastric.mucous.cells@meta.data$batch %in% c("NAG1","NAG2")] <- "NAG"
gastric.mucous.cells@meta.data$batch[gastric.mucous.cells@meta.data$batch %in% c("CAG1","CAG2","CAG3")] <- "CAG"
gastric.mucous.cells@meta.data$batch[gastric.mucous.cells@meta.data$batch %in% c("IMS1","IMS2","IMS3","IMS4","IMW1","IMW2")] <- "IM"
gastric.mucous.cells@meta.data$batch[gastric.mucous.cells@meta.data$batch %in% c("EGC")] <- "EGC"

gastric.mucous.cells  <- FindVariableGenes(object = gastric.mucous.cells, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.3, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = gastric.mucous.cells@var.genes)
gastric.mucous.cells  <- RunPCA(object =gastric.mucous.cells,pc.genes = gastric.mucous.cells@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
gastric.mucous.cells <- FindClusters(object = gastric.mucous.cells, reduction.type = "pca", dims.use = 1:20,
    resolution = c(3), print.output = 0, save.SNN = TRUE)
gastric.mucous.cells  <- RunTSNE(object = gastric.mucous.cells, dims.use = 1:20, do.fast = TRUE)


pdf("figure3f.pdf",width = 10,height = 8)
gg <- TSNEPlot(object = gastric.mucous.cells,do.label = F,pt.size = 3,do.return = T,group.by = "batch") 
gg <- gg + theme(legend.text = element_text(size = 25),
                 legend.position = 'right',
                 legend.key.height = unit(1.4,"cm"),
                 legend.key.width = unit(1.4,"cm"),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 20),
                 axis.title.x = element_text(size = 20),
                 axis.title.y = element_text(size = 20))
gg
dev.off()
}


## Supplementary_Figure3g-1 
## Supplementary_Figure3g-2
{
valid.cells <- gastric.mucous.cells@cell.names[grep("IM",gastric.mucous.cells@meta.data$batch)]
gland.IM.cell <- SubsetData(gastric.mucous.cells,cells.use = valid.cells,do.center = T,do.scale = T)
gland.IM.cell <- FindVariableGenes(object = gland.IM.cell, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 1,do.plot = T)
length(x = gland.IM.cell@var.genes)
gland.IM.cell <- RunPCA(object = gland.IM.cell,pc.genes = gland.IM.cell@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40) 
gland.IM.cell <- FindClusters(object = gland.IM.cell, reduction.type = "pca", dims.use = 1:10,
    resolution = c(1), print.output = 0, save.SNN = TRUE)
gland.IM.cell  <- RunTSNE(object = gland.IM.cell, dims.use = 1:10, do.fast = TRUE)

pca.gland.IM.cell <- data.frame(gland.IM.cell@dr$pca@cell.embeddings[,1:4])
pca.gland.IM.cell$original.cluster <- as.numeric(gland.IM.cell@ident)
pca.gland.IM.cell$samples <- gland.IM.cell@meta.data$batch
  
pca.gland.IM.cell$OLFM4 <- gland.IM.cell@data["OLFM4",]
pca.gland.IM.cell$ODAM <- gland.IM.cell@data["ODAM",]
pca.gland.IM.cell$MUC6 <- gland.IM.cell@data["MUC6",]
pca.gland.IM.cell$LEFTY1 <- gland.IM.cell@data["LEFTY1",]
pca.gland.IM.cell$TFF3 <- gland.IM.cell@data["TFF3",]

pdf("Supplementary_figure3g_2.pdf",width = 8,height = 8)
ggplot()+ geom_point(data = pca.gland.IM.cell,aes(x = PC1,y = PC2, colour = OLFM4,size = 20)) + scale_color_gradient(low = "gray",high = "red",na.value = "yellow")
dev.off()
pdf("Supplementary_figure3g_3.pdf",width = 8,height = 8)
ggplot()+ geom_point(data = pca.gland.IM.cell,aes(x = PC1,y = PC2, colour = MUC2,size = 20)) + scale_color_gradient(low = "gray",high = "red",na.value = "yellow")
dev.off()
pdf("Supplementary_figure3g_4.pdf",width = 8,height = 8)
ggplot()+ geom_point(data = pca.gland.IM.cell,aes(x = PC1,y = PC2, colour = LEFTY1,size = 20)) + scale_color_gradient(low = "gray",high = "red",na.value = "yellow")
dev.off()
pdf("Supplementary_figure3g_5.pdf",width = 8,height = 8)
ggplot()+ geom_point(data = pca.gland.IM.cell,aes(x = PC1,y = PC2, colour = ODAM,size = 20)) + scale_color_gradient(low = "gray",high = "red",na.value = "yellow")
dev.off()

kmeans.result <- kmeans(pca.gland.IM.cell[,c(1,2)],centers = 2)
pca.gland.IM.cell$cluster <- kmeans.result$cluster
d <- dist(pca.gland.IM.cell[,c(1,2)], method = "euclidean")
fit <- hclust(d, method="ward") 
groups <- cutree(fit, k=2)
pca.gland.IM.cell$cluster.new <- groups
 
pca.gland.IM.cell$col <- "blue"
pca.gland.IM.cell$col[which(pca.gland.IM.cell$cluster.new == 2)] <- "coral"
  
pdf("Supplementary_figure3g_1.pdf",width = 8,height = 8)
ggplot()+ geom_point(data = pca.gland.IM.cell,aes(x = PC1,y = PC2, colour = col,size = 30)) + theme(legend.text = element_blank(),
                 axis.text.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 25))
dev.off()

}


## figure3g
{
res <- FindMarkers(gland.IM.cell,ident.1 = c(1,4),ident.2 = c(0,3,2),test.use = "bimod",min.pct = 0.05)
res$Gene <- rownames(res)
# res["TFF3","p_val_adj"] <- 2.11456*exp(-85)
# res["LEFTY1","p_val_adj"] <- 4.606887*exp(-45)
pdf(file = "figure3g.pdf",width = 8,height = 8)
# Make a basic volcano plot
with(res, plot(avg_logFC, -log10(p_val_adj),cex.axis = 1.5,cex.lab = 1.5,  pch=16, cex = 2,xlab = "LogFC (Cluster 1 vs Cluster 2)",ylab = "-log10 (p_value_adjust)",xlim=c(-4.5,3.5)))

# Add colored points: red if padj<0.05, orange of log2FC>1, green if both)
with(subset(res, (avg_logFC)>1), points(avg_logFC, -log10(p_val_adj), pch=19, col="#FF8D8D"))
with(subset(res, (avg_logFC) < -1), points(avg_logFC, -log10(p_val_adj), pch=19, col="#789FCF"))

# Label points with the textxy function from the calibrate plot
library(calibrate)
res$symbol <- 0
res$symbol[which(rownames(res) %in% c("PGC","REG3A","LEFTY1","PHLDA1","ODAM","TFF3","TNFRSF12A","OLFM4"))] <- 1
with(subset(res, symbol == 1), textxy(avg_logFC, -log10(p_val_adj), labs=Gene, cex=1.5))
dev.off()

}


## Figure 3h
# co-expressing markers with MUC2 ,OLFM4 and LEFTY1
{
MUC6.cells <- dt@data["MUC6",]
MUC6.cells <- MUC6.cells[which(as.numeric(MUC6.cells)>0)]
MUC6.cell.names <- names(MUC6.cells)[-order(MUC6.cells)[1:round(0.1*length(MUC6.cells))]]

OLFM4.cells <- dt@data["OLFM4",]
OLFM4.cells <- OLFM4.cells[which(as.numeric(OLFM4.cells)>0)]
OLFM4.cell.names <- names(OLFM4.cells)[-order(OLFM4.cells)[1:round(0.1*length(OLFM4.cells))]]

LEFTY1.cells <- dt@data["LEFTY1",]
LEFTY1.cells <- LEFTY1.cells[which(as.numeric(LEFTY1.cells)>0)]
LEFTY1.cell.names <- names(LEFTY1.cells)[-order(LEFTY1.cells)[1:round(0.1*length(LEFTY1.cells))]]

valid.cells <- Reduce(union,c(MUC6.cell.names,OLFM4.cell.names,LEFTY1.cell.names))
# valid.cells <- intersect(valid.cells,dt@cell.names[which(dt@meta.data$batch %in% c("IMS3","IMS4","CAN2"))])
 
co.expressed.cells <- intersect(union(OLFM4.cell.names,LEFTY1.cell.names),MUC6.cell.names)
only.OLFM4.cells <- setdiff(OLFM4.cell.names,co.expressed.cells)
only.LEFTY1.cells <- setdiff(LEFTY1.cell.names,co.expressed.cells)
only.OLFM4.LEFTY1.cells <- only.OLFM4.cells 
only.MUC6.cells <- setdiff(MUC6.cell.names,co.expressed.cells)

result.matrix <- matrix(0,nrow = 5,ncol = 3)
rownames(result.matrix) <- c("NAG","CAG","IMW","IMS","EGC")
colnames(result.matrix) <- c("only.MUC6","only.OLFM4.LEFTY1","OLFM4.LEFTY1.MUC6")

NAG.cells <- intersect(dt@cell.names[grep("NAG",dt@meta.data$batch)],valid.cells)
result.matrix[1,] <- c(length(intersect(NAG.cells,only.MUC6.cells))/length(NAG.cells),length(intersect(NAG.cells,only.OLFM4.LEFTY1.cells))/length(NAG.cells),length(intersect(NAG.cells,co.expressed.cells))/length(NAG.cells))


CAG.cells <- intersect(dt@cell.names[grep("CAG",dt@meta.data$batch)],valid.cells)
result.matrix[2,] <- c(length(intersect(CAG.cells,only.MUC6.cells))/length(CAG.cells),length(intersect(CAG.cells,only.OLFM4.LEFTY1.cells))/length(CAG.cells),length(intersect(CAG.cells,co.expressed.cells))/length(CAG.cells))

IMW.cells <- intersect(dt@cell.names[grep("IMW",dt@meta.data$batch)],valid.cells)
result.matrix[3,] <- c(length(intersect(IMW.cells,only.MUC6.cells))/length(IMW.cells),length(intersect(IMW.cells,only.OLFM4.LEFTY1.cells))/length(IMW.cells),length(intersect(IMW.cells,co.expressed.cells))/length(IMW.cells))

IMS.cells <- intersect(dt@cell.names[grep("IMS",dt@meta.data$batch)],valid.cells)
result.matrix[4,] <- c(length(intersect(IMS.cells,only.MUC6.cells))/length(IMS.cells),length(intersect(IMS.cells,only.OLFM4.LEFTY1.cells))/length(IMS.cells),length(intersect(IMS.cells,co.expressed.cells))/length(IMS.cells))

EGS.cells <- intersect(dt@cell.names[grep("CAN",dt@meta.data$batch)],valid.cells)
result.matrix[5,] <- c(length(intersect(EGS.cells,only.MUC6.cells))/length(EGS.cells),length(intersect(EGS.cells,only.OLFM4.LEFTY1.cells))/length(EGS.cells),length(intersect(EGS.cells,co.expressed.cells))/length(EGS.cells))

 # MUC6 and OLFM4
result.matrix[,c(2,3)] <- result.matrix[,c(3,2)]
colnames(result.matrix) <- c("MUC6","MUC6+OLFM4","OLFM4")
result.matrix[,c(1,3)] <- result.matrix[,c(3,1)]
colnames(result.matrix)[c(1,3)] <- colnames(result.matrix)[c(3,1)]
library(reshape2)
dt.data <- melt(result.matrix)
colnames(dt.data)[2] <- "Class"
pdf("Figure3h.pdf",width = 10,height = 6)
  g <- ggplot(dt.data, aes(fill=Class, y=value, x=Var1)) +
    geom_bar( stat="identity") +  scale_fill_manual(values=c(rainbow(36)[22], "#E69F00", "springgreen"))
  g <- g+theme(
    panel.background=element_rect(colour = NA, fill = "white"),
    axis.title.y = element_text(size=15,face = "bold"),
    axis.title.x = element_blank(),
    axis.ticks.length=unit(.4,"lines"),
    axis.line=element_line(colour="black"),
    axis.text = element_text(face = "bold"),
    axis.text.x=element_text(angle = 45,vjust = 0.5,size = 20),
    axis.text.y=element_text(size = 12),
    legend.key = element_rect(colour = NA,fill = "white"),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.title.align=0.2,
    strip.background =element_rect(fill="beige"),
    strip.text.x = element_text(size = 15,face = "bold")
  )
  g
dev.off()


# gastric mucous cells in pit & gland
LGR5.cells <- dt@data["LGR5",]
LGR5.cells <- LGR5.cells[which(as.numeric(LGR5.cells)>0)]
LGR5.cell.names <- names(LGR5.cells)[-order(LGR5.cells)[1:round(0.1*length(LGR5.cells))]]

proportion.LGR5.CAG.cells <- length(intersect(dt@cell.names[grep("CAG",dt@meta.data$batch)],LGR5.cell.names))/length(dt@cell.names[grep("CAG",dt@meta.data$batch)])
proportion.LGR5.IMW.cells <- length(intersect(dt@cell.names[grep("IMW",dt@meta.data$batch)],LGR5.cell.names))/length(dt@cell.names[grep("IMW",dt@meta.data$batch)])
proportion.LGR5.IMS.cells <- length(intersect(dt@cell.names[grep("IMS",dt@meta.data$batch)],LGR5.cell.names))/length(dt@cell.names[grep("IMS",dt@meta.data$batch)])
proportion.LGR5.EGC.cells <- length(intersect(dt@cell.names[grep("CAN",dt@meta.data$batch)],LGR5.cell.names))/length(dt@cell.names[grep("CAN",dt@meta.data$batch)])

}


## Supplementary_Figure 3h

# assessing the co-expression of MUC6 and OLFM4 in dividual cells
# MUC6 VS OLFM4 VS CDX2
{
MUC6.cells <- dt@data["HES6",]
MUC6.cells <- MUC6.cells[which(as.numeric(MUC6.cells)>0)]
MUC6.cell.names <- names(MUC6.cells)[-order(MUC6.cells)[1:round(0.1*length(MUC6.cells))]]

OLFM4.cells <- dt@data["SOX4",]
OLFM4.cells <- OLFM4.cells[which(as.numeric(OLFM4.cells)>0)]
OLFM4.cell.names <- names(OLFM4.cells)[-order(OLFM4.cells)[1:round(0.1*length(OLFM4.cells))]]

gastric.mucous.cells <- SubsetData(dt,ident.use = c("GMC"),do.center = T,do.scale = T)
IM.cells <- gastric.mucous.cells@cell.names[grep("IM",gastric.mucous.cells@meta.data$batch)]
# IM.cells <- gastric.mucous.cells@cell.names
gland.IMS.cell.data <- SubsetData(dt,cells.use = IM.cells,do.center = T,do.scale = T)
gland.IMS.cell.data <- as.data.frame(as.matrix(gland.IMS.cell.data@data))
gland.IMS.cell.data <- t(gland.IMS.cell.data)

valid.MUC6 <- which(as.numeric(gland.IMS.cell.data[,"MUC6"]) != 0)
valid.MUC6.cells <- intersect(rownames(gland.IMS.cell.data)[valid.MUC6],MUC6.cell.names)

valid.OLFM4 <- which(as.numeric(gland.IMS.cell.data[,"OLFM4"]) != 0)
valid.OLFM4.cells <- intersect(rownames(gland.IMS.cell.data)[valid.OLFM4],OLFM4.cell.names)
valid.CDX2 <- which(as.numeric(gland.IMS.cell.data[,"CDX2"]) != 0)

valid.cells <- intersect(valid.MUC6,valid.OLFM4)
gland.IM.data.temp <- gland.IMS.cell.data[valid.cells,]
gland.IM.data.temp <- gland.IM.data.temp[-which(gland.IM.data.temp[,"MUC6"]>5),]
gland.IM.data.temp <- gland.IM.data.temp[-which(gland.IM.data.temp[,"MUC6"]<2),]
gland.IM.data.temp <- gland.IM.data.temp[-which(gland.IM.data.temp[,"OLFM4"]>4),]
gland.IM.data.temp <- as.data.frame(gland.IM.data.temp)


pdf("Supplementary_Figure3h_HES6+SOX4.pdf",width = 8,height = 8)
par(mar = c(2,3,4,4))
p <- ggplot(gland.IM.data.temp, aes(MUC6,OLFM4)) + geom_point(size = 5) + geom_smooth(method = "lm", se = FALSE) + theme(
                 panel.background = element_blank(),
                 axis.text.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 25))
p
dev.off()
dd <- (summary(lm(OLFM4~MUC6,data = gland.IM.data.temp)))
cor.test(gland.IM.data.temp[,"SOX4"],gland.IM.data.temp[,"HES6"])

stem.cells.data <- SubsetData(dt,ident.use = c("Stem-like cell"),do.center = T,do.scale = T)
  
CDX2.cells <- stem.cells.data@data["CDX2",]
CDX2.cells <- CDX2.cells[which(as.numeric(CDX2.cells)>0)]
CDX2.cell.names <- names(CDX2.cells)[-order(CDX2.cells)[1:round(0.1*length(CDX2.cells))]]

OLFM4.cells <- stem.cells.data@data["OLFM4",]
OLFM4.cells <- OLFM4.cells[which(as.numeric(OLFM4.cells)>0)]
OLFM4.cell.names <- names(OLFM4.cells)[-order(OLFM4.cells)[1:round(0.1*length(OLFM4.cells))]]
 
valid.cells <- intersect(CDX2.cell.names,OLFM4.cell.names)
data.temp <- SubsetData(stem.cells.data,cells.use = valid.cells,do.center = T,do.scale = T) 
  
data.temp.1 <- t(data.frame(as.matrix(data.temp@data)))
cor.test(as.numeric(data.temp.1[,"CDX2"]),as.numeric(data.temp.1[,"OLFM4"]))
  
data.temp.1 <- data.frame(data.temp.1)
pdf("Supplementary_Figure3h_CDX2+OLFM4.pdf",width = 8,height = 8)
par(mar = c(2,3,4,4))
p <- ggplot(data.temp.1, aes(CDX2,OLFM4)) + geom_point(size = 5) + geom_smooth(method = "lm", se = FALSE) + theme(
                 panel.background = element_blank(),
                 axis.text.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 25))
p
dev.off()
cor.test(data.temp.1[,"CDX2"],data.temp.1[,"OLFM4"])
}
```




## Figure 4 ##
# Focusing Enteroendocrine cells
```{r}

library(dplyr)
library(reshape2)
library(Seurat)

load("data.temp.res.2.annotation.Rdata")

{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}

## Figure 4a
pdf(file = "Figure4a.pdf", width = 8,height = 8)
data.raw <- data.frame(dt@dr$tsne@cell.embeddings)
data.raw$Cancer <- 1
data.raw$Cancer[which(dt@ident == "Enteroendocrine")] <- 0
ggplot(data.raw, aes(x=tSNE_1, y=tSNE_2, color=factor(Cancer))) + geom_point() + scale_colour_manual(values = c("#FF8AB1","#C4C4C4")) + theme(
  legend.position = "none",
  axis.title.x = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 30))
dev.off()



## Figure 4b
#preliminary work:find cluster
# load("/Users/zhangyiding/Documents/singlecell/SecondEdition/data.temp.res.2.annotation.Rdata")
#####extract edocrine cell
endoc<- SubsetData(dt,ident.use ="Enteroendocrine",do.scale = T,do.center = T)

endoc<- FindVariableGenes(object = endoc, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.2, x.high.cutoff = 5, y.cutoff = 1,do.plot = T)

length(x = endoc@var.genes)
# cancer.cell <- ScaleData(object = cancer.cell, do.scale = T,do.center = T,vars.to.regress = c("nUMI", "percent.mito"))
endoc<- RunPCA(object =endoc, pc.genes = endoc@var.genes, do.print = TRUE, pcs.print = 1:5,genes.print = 5,pcs.compute = 120)

endoc <- FindClusters(object = endoc, reduction.type = "pca", dims.use = 1:20,
    resolution = 5, print.output = 0, save.SNN = TRUE)

endoc<- RunTSNE(object =endoc, dims.use = 1:20, do.fast = TRUE)

library(RColorBrewer)
color<-c(brewer.pal(8,"Paired"),brewer.pal(11,"RdGy")[1:3],brewer.pal(11,"Set2")[1:5],brewer.pal(8,"Dark2")[1:8])
barplot(rep(1,24),col=color)

pdf("endoc.pdf")
TSNEPlot(object = endoc,do.label = T,label.size = 5,colors.use = color)
dev.off()

TSNEPlot(object = endoc,do.label = T,group.by="lesion",label.size = 5,colors.use = color)



# 
# saveRDS(endoc,file = "./endoc_res5.rds")
endoc<-readRDS("endoc_res5.rds")
# endoc@meta.data$batch[which(endoc@meta.data$batch=="CAN2")]<-"EGC"
#stomach endocrine cell marker
FeaturePlot(object = endoc,features.plot = "GAST",cols.use = c("grey", "blue"), reduction.use = "tsne")#G cell
FeaturePlot(object = endoc,features.plot = "GHRL",cols.use = c("grey", "blue"), reduction.use = "tsne")#X cell
FeaturePlot(object = endoc,features.plot = "SST",cols.use = c("grey", "blue"), reduction.use = "tsne")#D cell
FeaturePlot(object = endoc,features.plot = "DCLK1",cols.use = c("grey", "blue"), reduction.use = "tsne")#Tuft cell
FeaturePlot(object = endoc,features.plot = "SLC6A4",cols.use = c("grey", "blue"), reduction.use = "tsne")#EC cell:Enterochromaffin (EC) cells

#enteroendocrine subtype marker
FeaturePlot(object = endoc,features.plot = "SCT",cols.use = c("grey", "blue"), reduction.use = "tsne")#S cell
FeaturePlot(object = endoc,features.plot = "CCK",cols.use = c("grey", "blue"), reduction.use = "tsne")#I cell
FeaturePlot(object = endoc,features.plot = "GCG",cols.use = c("grey", "blue"), reduction.use = "tsne")#L cell
FeaturePlot(object = endoc,features.plot = "GIP",cols.use = c("grey", "blue"), reduction.use = "tsne")#K cell
FeaturePlot(object = endoc,features.plot = "SST",cols.use = c("grey", "blue"), reduction.use = "tsne")#D cell
FeaturePlot(object = endoc,features.plot = "NTS",cols.use = c("grey", "blue"), reduction.use = "tsne")#N cell
FeaturePlot(object = endoc,features.plot = "GHRL",cols.use = c("grey", "blue"), reduction.use = "tsne")#A cell
FeaturePlot(object = endoc,features.plot = "TPH1",cols.use = c("grey", "blue"), reduction.use = "tsne")#EC cell_enterochromaffin cells
TSNEPlot(object = endoc,do.label = T,label.size = 5,colors.use = color)
refactor<-function(name,from,to){
  output<-as.character(name)
  output[which(output %in% from)]<-to
  output<-as.factor(output)
  names(output)<-names(name)
  return (output)
}

###relabel

cell<-refactor(endoc@ident,c(17,8,23),"DA/X")
cell<-refactor(cell,18,"DILK")
cell<-refactor(cell,c(11,15),"NA")
cell<-refactor(cell,0,"L")
cell<-refactor(cell,c(9,3),"D")
cell<-refactor(cell,c(14,1),"G-EC")
cell<-refactor(cell,c(16,10,12,13,2,20,21,22,4,5,6,7),"G")
cell<-refactor(cell,19,"GD")

endoc@meta.data$subtype<-cell
endoc@ident<-cell

library(RColorBrewer)
color<-c(brewer.pal(8,"Paired"),brewer.pal(11,"RdGy")[1:3],brewer.pal(11,"Set2")[1:5],brewer.pal(8,"Dark2")[1:8])
pdf("endoc-1.pdf")
TSNEPlot(object = endoc,do.label = T,label.size = 5,colors.use = color)
dev.off()

pdf("endoc-2.pdf")
TSNEPlot(object = endoc,do.label = T,group.by="batch",colors.use = color)
dev.off()

current.cluster.ids <- c("D", "DA/X" ,"DILK", "G", "G-EC", "GD" ,"L", "NA" )
 new.cluster.ids <- c( "C1","C2","C3","C4","C5","C6","C7","C8")
 endoc@ident <- plyr::mapvalues(x = endoc@ident, from = current.cluster.ids, to = new.cluster.ids)

marker.use<-c("GAST","SST","TPH1","GHRL","GIP","SCT","CCK","NTS","GCG")
matrix<-(x = FetchData(object = endoc, vars.all =marker.use,use.imputed = FALSE))

#get num of each cell type
data.plot <-data.frame()
for (cluster in levels(endoc@ident)){
  ind<-which(endoc@ident == cluster)
  num.total<-length(ind)
  data.use<-matrix[ind,]
  data.use[data.use<3]<-0
  left.ind<-which(rowSums(data.use)!=0)
  data.use<-data.frame(data.use[left.ind,])
  order<-apply(data.use,1,function(x) which(x==max(x)))
  factor<-factor(marker.use[as.numeric(unlist(order))],levels = marker.use)
  data<-data.frame(table(factor))
  data[,2]<-data[,2]/num.total
  colnames(data)<-c("Marker","Proportion")
  data$cell<-rep(cluster,length(marker.use))
  data.plot<-rbind(data.plot,data)
}

# p<-ggplot(data=data.plot, aes(x=cell, y=Fre, fill=Marker))+
#   geom_bar(stat="identity")+
#   theme(legend.position="bottom")+
#   guides(fill = guide_legend(direction = "horizontal"))+
#    scale_fill_manual(values=brewer.pal(12,"Paired")[1:9])

pdf("fig4b.pdf",width=8, height=8)

p<-TSNEPlot(object = endoc,do.label = T,label.size = 12,pt.size = 2,do.return = TRUE)+theme(
  legend.position = "none",
  axis.title.x = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 30))
p
dev.off()


## Figure 4c
library(easyGgplot2)
marker.use<-c("GAST","SST","CCK") 
#marker.use<-c("SST","TPH1","GHRL","GIP","SCT","CCK","NTS","GCG") 
matrix<-as.matrix(x = FetchData(object = endoc, vars.all =marker.use,,use.imputed = FALSE))

cell.1<-which(endoc@meta.data$batch %in% "NAG1")
cell.2<-which(endoc@meta.data$batch %in% "NAG2")
cell.3<-which(endoc@meta.data$batch %in% "CAG3")
cell.4<-which(endoc@meta.data$batch %in% "CAG2")
cell.5<-which(endoc@meta.data$batch %in% "CAG1")
cell.6<-which(endoc@meta.data$batch %in% "IMW2")
cell.7<-which(endoc@meta.data$batch %in% "IMW1")
cell.8<-which(endoc@meta.data$batch %in% "IMS1")
cell.9<-which(endoc@meta.data$batch %in% "IMS2")
cell.10<-which(endoc@meta.data$batch %in% "IMS3")
cell.11<-which(endoc@meta.data$batch %in% "IMS4")
cell.12<-which(endoc@meta.data$batch %in% "CAN1")

cell.stage<-list(cell.1=cell.1,cell.2=cell.2,cell.3=cell.3,cell.4=cell.4,cell.5=cell.5,cell.6=cell.6,cell.7=cell.7,cell.8=cell.8,cell.9=cell.9,cell.10=cell.10,cell.11=cell.11,cell.12=cell.12)

Mean<-matrix(0,ncol=length(marker.use),nrow=length(cell.stage))
rownames(Mean)<-c("NAG1","NAG2","CAG3","CAG2","CAG1","IMW2","IMW1","IMS1","IMS2","IMS3","IMS4","CAN1")
colnames(Mean)<-marker.use

for (i in 1:length(marker.use)){
  for (j in 1:length(cell.stage))
 Mean[j,i]<- apply( data.frame( matrix[cell.stage[[j]],marker.use[i] ]),2,mean)
}

 
gene<-factor(rep(marker.use,each=length(cell.stage)),levels =c("GAST","SST","CCK") )
#gene<-factor(rep(marker.use,each=length(cell.stage)),levels =c("SST","TPH1","GHRL","GIP","SCT","CCK","NTS","GCG") )

stage<-factor(rep(c("NAG1","NAG2","CAG3","CAG2","CAG1","IMW2","IMW1","IMS1","IMS2","IMS3","IMS4","CAN1"),length(marker.use)),level=c("NAG1","NAG2","CAG3","CAG2","CAG1","IMW2","IMW1","IMS1","IMS2","IMS3","IMS4","CAN1"))
Expression<-as.vector(Mean)


library(RColorBrewer)


c<-c(brewer.pal(9,"Blues")[c(7,2)],brewer.pal(9,"Spectral")[2])
#c<-c(brewer.pal(9,"Blues")[c(2)],brewer.pal(9,"BrBG")[4],brewer.pal(9,"Spectral")[5:1],brewer.pal(9,"RdGy")[1])
color<-rep(c,each=length(cell.stage))

df<-data.frame(Gene=gene,Stage=stage,Expression=Expression,color=color)


### set breaks number to 3 for each facets
require(grid)
# defining the breaks function, 
# s is the scaling factor (cf. multiplicative expand)
equal_breaks <- function(n = 3, s = 0.05, ...){
  function(x){
    # rescaling
    d <- s * diff(range(x)) / (1+2*s)
    round(seq(min(x)+d, max(x)-d, length=n), 1)
  }
}

g<-ggplot(df,aes(x=Stage,y=Expression))
g<-g+geom_bar(stat="identity",fill=color, colour="gray65")+
 scale_y_continuous(breaks=equal_breaks(n=3, s=0.05),expand = c(0, 0))+
 facet_grid(Gene~.)+ 
 theme(title = element_text(size=20),strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 45, hjust = 1,size=15),axis.text.y = element_text( color = "black", size = 15),strip.text = element_text(face="bold", size=15,lineheight=10.0))
#facet_grid(Gene~., scales="free_y"

#set background color for each facet
g <- ggplot_gtable(ggplot_build(g))
stripr <- which(grepl('strip-r', g$layout$name))
fills <- c("#2171B5","#6BAED6" ,"#DEEBF7"  ,"#FDAE61" )
#fills <- c("#6BAED6" ,"#DEEBF7" ,"#F6E8C3" ,"#FFFFBF", "#FEE08B" ,"#FDAE61","#F46D43" ,"#D53E4F", "#B2182B")

k <- 1
for (i in stripr) {
j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
k <- k+1
}


pdf("./fig4c_1.pdf",width = 7,height = 7)
#pdf("./fig4c_2.pdf",width = 7,height = 7)
grid.draw(g)#need run in console
dev.off()

## Figure 4d-1/4d-2

pdf("Figure4d_1.pdf",width=8, height=6)
colnames(data.plot)[3] <- "Cells"
p<-ggplot(data=data.plot, aes(x=Cells, y=Proportion, fill=Marker))+ geom_bar(stat="identity") + theme(legend.position="top",legend.key.size =  unit(0.4, "in"),title = element_text(size=15),axis.text.x = element_text(size=25),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25),legend.text=element_text(size=20),legend.title = element_text(size = 20)) +  guides(fill = guide_legend(direction = "horizontal"))+ scale_fill_manual(values=c("#6BAED6" ,"#DEEBF7" ,"#F6E8C3" ,"#FFFFBF", "#FEE08B" ,"#FDAE61","#F46D43" ,"#D53E4F", "#B2182B"))
p
dev.off()


 # theme(title = element_text(size=20),strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 45, hjust = 1,size=15),axis.text.y = element_text( color = "black", size = 15),strip.text = element_text(face="bold", size=15,lineheight=10.0))
current.cluster.ids <- c("D", "DA/X" ,"DILK", "G", "G-EC", "GD" ,"L", "NA" )
 new.cluster.ids <- c( "C1","C2","C3","C4","C5","C6","C7","C8")
 endoc@ident <- plyr::mapvalues(x = endoc@ident, from = current.cluster.ids, to = new.cluster.ids)
C1<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C1")]))
C1$Freq<-C1$Freq/sum(C1$Freq)
C1$Cluster<-"C1"
C2<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C2")]))
C2$Freq<-C2$Freq/sum(C2$Freq)
C2$Cluster<-"C2"
C3<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C3")]))
C3$Freq<-C3$Freq/sum(C3$Freq)
C3$Cluster<-"C3"
C4<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C4")]))
C4$Freq<-C4$Freq/sum(C4$Freq)
C4$Cluster<-"C4"
C5<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C5")]))
C5$Freq<-C5$Freq/sum(C5$Freq)
C5$Cluster<-"C5"
C6<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C6")]))
C6$Freq<-C6$Freq/sum(C6$Freq)
C6$Cluster<-"C6"
C7<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C7")]))
C7$Freq<-C7$Freq/sum(C7$Freq)
C7$Cluster<-"C7"
C8<-data.frame(table(endoc@meta.data$batch[which(endoc@ident=="C8")]))
C8$Freq<-C8$Freq/sum(C8$Freq)
C8$Cluster<-"C8"
data.cluster<-rbind(C1,C2,C3,C4,C5,C6,C7,C8)
colnames(data.cluster)<-c("Sample","Proportion","Cluster")

data.cluster$Sample <- gsub("CAN2","EGC",as.character(data.cluster$Sample))

data.cluster$Sample[grep("IMW",data.cluster$Sample)] <- "IMW"
data.cluster$Sample[grep("IMS",data.cluster$Sample)] <- "IMS"
data.cluster$Sample[grep("NAG",data.cluster$Sample)] <- "NAG"
data.cluster$Sample[grep("CAG",data.cluster$Sample)] <- "CAG"

data.cluster$Sample <-factor(data.cluster$Sample,levels=c("NAG","CAG","IMW","IMS","EGC"))

color <- c("#2B61AD","#AD5FA2","#F39E5B","#F4A27F","#DB4623")

pdf("Figure4d_2.pdf",width=8, height=6)
p<-ggplot(data=data.cluster, aes(x=Cluster, y=Proportion, fill=Sample))+
  geom_bar(stat="identity")+
  theme(legend.position="bottom",legend.key.size =  unit(0.4, "in"),title = element_text(size= 25),axis.text.x = element_text(size= 25),axis.text.y = element_text(size = 25),legend.text=element_text(size=20))+  guides(fill = guide_legend(direction = "horizontal"))+scale_fill_manual(values=color)
p
dev.off()

pdf("Figure4d_1.pdf",width=8, height=6)
colnames(data.plot)[3] <- "Cells"
p<-ggplot(data=data.plot, aes(x=Cells, y=Proportion, fill=Marker))+ geom_bar(stat="identity") + theme(legend.position="top",legend.key.size =  unit(0.4, "in"),title = element_text(size=25),axis.text.x = element_text(size=25),axis.title.x = element_text(size = 25),axis.title.y = element_text(size = 25), axis.text.y = element_text(size = 25),legend.text=element_text(size=20)) +  guides(fill = guide_legend(direction = "horizontal"))+ scale_fill_manual(values=c("#6BAED6" ,"#DEEBF7" ,"#F6E8C3" ,"#FFFFBF", "#FEE08B" ,"#FDAE61","#F46D43" ,"#D53E4F", "#B2182B"))
p
dev.off()

## Figure 4e
#OR51E1
pdf("4e.pdf",width = 4,height=4)
FeaturePlot(endoc,features.plot = "OR51E1",cols.use = c("grey","red"),pt.size = 1.5)
dev.off()

```





## Figure 5 ##
## Focusing the goblet cells
```{r}

library(dplyr)
library(reshape2)
library(Seurat)

load("data.temp.res.2.annotation.Rdata")

{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}

## Figure 5a

{
Goblet.cells <- SubsetData(dt,ident.use = c("Goblet cell"),do.center = T,do.scale = T)
Goblet.cells.list <- Goblet.cells@cell.names
# CAG.IM.cells.list <- dt@cell.names[-which(dt@meta.data$batch %in% c("CAN2"))]
valid.cells <- Goblet.cells.list

pdf(file = "Figure5a-1.pdf", width = 8,height = 8)
data.raw <- data.frame(dt@dr$tsne@cell.embeddings)
data.raw$Goblet.cells <- 1
data.raw$Goblet.cells[which(dt@cell.names %in% valid.cells)] <- 0
ggplot(data.raw, aes(x=tSNE_1, y=tSNE_2, color=factor(Goblet.cells))) + geom_point() + scale_colour_manual(values = c("#45B500","#C4C4C4")) + theme(
  legend.position = "none",
  axis.title.x = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 30))

dev.off()
}

{
  
pdf(file = "Figure5a-2.pdf", width = 8, height = 8)
gg <- TSNEPlot(object = goblet.cell,do.label = T,label.size = 25,pt.size = 4,do.return = T) 
gg <- gg + theme(
  # panel.border = element_blank(),
  # axis.line = element_line(),
  legend.text = element_text(size = 25),
                 legend.position = 'bottom',
                 legend.key.height = unit(1.4,"cm"),
                 legend.key.width = unit(1.4,"cm"),
                 axis.text.x = element_text(size = 30),
                 axis.text.y = element_text(size = 30),
                 axis.title.x = element_text(size = 30),
                 axis.title.y = element_text(size = 30))
gg
dev.off()  
}


## Supplementary figures 
{
goblet.batch <- table(Goblet.cells@meta.data$batch)[c("NAG2","CAG2","CAG3","IMW1","IMW2","IMS1","IMS2","IMS3","IMS4","EGC")]
overall.batch <- table(dt@meta.data$batch)[names(goblet.batch)]
propor.cells <- goblet.batch/overall.batch
# rotate_x <- function(data, column_to_plot, labels_vec, rot_angle) {
#     plt <- barplot(data[[column_to_plot]], col='steelblue', xaxt="n")
#     text(plt, par("usr")[3], labels = labels_vec, srt = rot_angle, adj = c(1.1,1.1), xpd = TRUE, cex=0.6) 
# }
end_point <- 0.5 + 2*length(propor.cells) -1
pdf("Supplementary_Figure5a-Goblet-cells.pdf",width = 6,height = 6)
barplot(propor.cells,ylab = "Proportion of goblet cells",ylim = c(0,0.08),xlab = "",space = 1, main = "")
text(seq(0.5,end_point,by = 2), 
     srt = 70, adj = 1, xpd = TRUE,
     labels = names(propor.cells), cex=0.65)
dev.off() 
}


# Figure 5b
library(ggcorrplot)
{
marker.goblet.genes <- FindAllMarkers(goblet.cell,only.pos = T,test.use = "roc")
goblet.cell.data.raw <- t(data.frame(as.matrix(goblet.cell@data[unique(marker.goblet.genes$gene),])))
p.mat <- round(cor(goblet.cell.data.raw), 1)
p.mat <- data.frame(p.mat); 
colnames(p.mat) <- rownames(p.mat) <- unique(marker.goblet.genes$gene)
pdf("Figure5b.pdf",width = 18,height = 18)
ggcorrplot(p.mat, hc.order = TRUE, outline.col = "white",tl.cex = 2)
dev.off()  
}



## Figure 5c
marker.goblet.genes <- FindMarkers(goblet.cell,ident.1 = c(0),ident.2 = c(1,2,3,4))
write.csv(marker.goblet.genes,file = "marker.goblet.genes.csv")
cluster1 <- data.frame(as.vector(c(3.2E-6,6.3E-4,1.4E-3,9.2E-3,2E-2)))
cluster1 <- -log10(cluster1)
rownames(cluster1) <- c("Metabolic pathways","Oxidative phosphorylation","Fat digestion and absorption","Terpenoid backbone biosynthesis","Amino sugar and nucleotide sugar metabolism")
colnames(cluster1) <- "pvalue"

pdf("Figure5c-1.pdf",width = 20,height = 20)
par(mar=c(2,8,3,3))
barplot(cluster1$pvalue,cex.axis = 6,cex.names = 2,ylim = c(0,6),axis.lty = 2,beside = T)
dev.off()

cluster0 <- data.frame(as.vector(c(3.2E-27,6.3E-20,1.4E-7,9.2E-6,2E-4)))
cluster0 <- -log10(cluster0)
rownames(cluster1) <- c("Metabolic pathways","Oxidative phosphorylation","Fat digestion and absorption","Terpenoid backbone biosynthesis","Amino sugar and nucleotide sugar metabolism")
colnames(cluster0) <- "pvalue"


pdf("Figure5c-2.pdf",width = 20,height = 20)
par(mar=c(2,8,3,3))
barplot(cluster0$pvalue,cex.axis = 6,cex.names = 2,ylim = c(0,30),axis.lty = 2,beside = T)
dev.off()



## Figure 5d
genes <- c("MUC2")
pdf("Figure4e-a.pdf",width = 8,height = 8)
gg <- VlnPlot(goblet.cell,genes,x.lab.rot = 90,remove.legend = T,point.size.use = 0,size.title.use = 30,do.return = T)
gg <- gg + theme(axis.text.x = element_text(size = 35),axis.text.y = element_text(size = 35),axis.ticks.x = element_line(size = 30),axis.title.x = element_text(size = 30),title = element_text(size = 40),axis.title = element_text(size = 0))
gg
dev.off()


genes <- c("HES6")
pdf("Figure4e-b.pdf",width = 8,height = 8)
gg <- VlnPlot(goblet.cell,genes,x.lab.rot = 90,remove.legend = T,point.size.use = 0,size.title.use = 30,do.return = T)
gg <- gg + theme(axis.text.x = element_text(size = 35),axis.text.y = element_text(size = 35),axis.ticks.x = element_line(size = 30),axis.title.x = element_text(size = 30),title = element_text(size = 40),axis.title = element_text(size = 0))
gg
dev.off()


## Figure 5e
total.goblet.cells <- goblet.cell@cell.names

low.10.percet.cells <- as.numeric(as.vector(goblet.cell@data))
low.10.percet.value <- quantile(low.10.percet.cells[which(low.10.percet.cells>0)],0.1)

MUC2.cells <- goblet.cell@cell.names[which(goblet.cell@data["MUC2",]>low.10.percet.value)]
HES6.cells <- goblet.cell@cell.names[which(goblet.cell@data["HES6",]>low.10.percet.value)]
HEPACAM2.cells <- goblet.cell@cell.names[which(goblet.cell@data["HEPACAM2",]>low.10.percet.value)]
KLK12.cells <- goblet.cell@cell.names[which(goblet.cell@data["KLK12",]>low.10.percet.value)]
MUC2.HEPACAM2.HES6.cells <- goblet.cell@cell.names[which(colSums(as.matrix(goblet.cell@data)[c("HEPACAM2","HES6","MUC2","KLK12"),])>low.10.percet.value)]

a1 <- length(MUC2.cells)/length(total.goblet.cells)
a2 <- length(HES6.cells)/length(total.goblet.cells)
# a3 <- length(HEPACAM2.cells)/length(total.goblet.cells)
# a4 <- length(KLK12.cells)/length(total.goblet.cells)
a5 <- length(MUC2.HEPACAM2.HES6.cells )/length(total.goblet.cells)

final.data <- as.vector(c(a1,a2,a5))
names(final.data) <- c("MUC2","HES6","MUC2+HES6")
final.data <- final.data[c("MUC2+HES6","HES6","MUC2")]
pdf(file = "figure5e.pdf",width = 8,height = 8)
barplot(final.data, main = NULL, horiz =T, xlab = "The proportion of cells",col = c("#FF9999","gray","gray","gray","gray"),cex.axis = 2,cex.names = 4,axisnames = F,xlim = c(0,1),width = 0.1)
dev.off()

## Supplementary figure SOX4 + HES6
HES6.cells <- dt@data["HES6",]
HES6.cells <- HES6.cells[which(as.numeric(HES6.cells)>0)]
HES6.cell.names <- names(HES6.cells)[-order(HES6.cells)[1:round(0.1*length(HES6.cells))]]

SOX4.cells <- dt@data["SOX4",]
SOX4.cells <- SOX4.cells[which(as.numeric(SOX4.cells)>0)]
SOX4.cell.names <- names(SOX4.cells)[-order(SOX4.cells)[1:round(0.1*length(SOX4.cells))]]

gastric.mucous.cells <- SubsetData(dt,ident.use = c("Goblet cell"),do.center = T,do.scale = T)
IM.cells <- gastric.mucous.cells@cell.names[grep("IM",gastric.mucous.cells@meta.data$batch)]

# IM.cells <- gastric.mucous.cells@cell.names
gland.IMS.cell.data <- SubsetData(dt,cells.use = IM.cells,do.center = T,do.scale = T)
gland.IMS.cell.data <- as.data.frame(as.matrix(gland.IMS.cell.data@data))
gland.IMS.cell.data <- t(gland.IMS.cell.data)

valid.HES6 <- which(as.numeric(gland.IMS.cell.data[,"HES6"]) != 0)
valid.HES6.cells <- intersect(rownames(gland.IMS.cell.data)[valid.HES6],HES6.cell.names)

valid.SOX4 <- which(as.numeric(gland.IMS.cell.data[,"SOX4"]) != 0)
valid.SOX4.cells <- intersect(rownames(gland.IMS.cell.data)[valid.SOX4],SOX4.cell.names)

valid.cells <- intersect(valid.HES6,valid.SOX4)
gland.IM.data.temp <- gland.IMS.cell.data[valid.cells,]
# gland.IM.data.temp <- gland.IM.data.temp[-which(gland.IM.data.temp[,"HES6"]>5),]
# gland.IM.data.temp <- gland.IM.data.temp[-which(gland.IM.data.temp[,"HES6"]<2),]
# gland.IM.data.temp <- gland.IM.data.temp[-which(gland.IM.data.temp[,"SOX4"]>4),]
gland.IM.data.temp <- as.data.frame(gland.IM.data.temp)
pdf("Supplementary_Figure3h_HES6+SOX4.pdf",width = 8,height = 8)
par(mar = c(2,3,4,4))
p <- ggplot(gland.IM.data.temp, aes(HES6,SOX4)) + geom_point(size = 5) + geom_smooth(method = "lm", se = FALSE) + theme(
                 panel.background = element_blank(),
                 axis.text.x = element_text(size = 25),
                 axis.text.y = element_text(size = 25),
                 axis.title.x = element_text(size = 25),
                 axis.title.y = element_text(size = 25))
p
dev.off()
dd <- (summary(lm(SOX4~HES6,data = gland.IM.data.temp)))
cor.test(gland.IM.data.temp[,"SOX4"],gland.IM.data.temp[,"HES6"])

```







## Figure 6 ##
# Focusing on the Cancer cell

```{r }

library(dplyr)
library(reshape2)
library(Seurat)
library(pheatmap)

load("data.temp.res.2.annotation.Rdata")
{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}

{
TSNEPlot(object = Cancer.cells,do.label = T,label.size = 5,pt.size = 2)
TSNEPlot(object = Cancer.cells,group.by = "batch",pt.size = 2)

}


## Figure 6a

pdf(file = "Figure6a.pdf", width = 8,height = 8)
data.raw <- data.frame(dt@dr$tsne@cell.embeddings)
data.raw$Cancer <- 1
data.raw$Cancer[which(dt@ident == "Cancer cell")] <- 0
ggplot(data.raw, aes(x=tSNE_1, y=tSNE_2, color=factor(Cancer))) + geom_point() + scale_colour_manual(values = c("#71C7F5","#C4C4C4")) + theme(
  legend.position = "none",
  axis.title.x = element_text(size = 30),
  axis.title.y = element_text(size = 30),
  axis.text.x = element_text(size = 30),
  axis.text.y = element_text(size = 30))
dev.off()


## Figure 6b -- Comparing the cancer and adjacent IM mucosa
library(ggplot2)

{

cancer.markers <- c("CEACAM5","CEACAM6","BAX","CCND2","SERPINB5")

Epithelial.cells.all.data <- (data.frame(as.matrix(Epithelial.cells.all@data[cancer.markers,])))
Epithelial.cells.all.data <- t(Epithelial.cells.all.data)
Epithelial.cells.all.data <- data.frame(Epithelial.cells.all.data)
Epithelial.cells.all.data$ident <- as.vector(as.character(Epithelial.cells.all@ident))

pdf("Figure6b-1.pdf",width = 16,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, CEACAM6,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60),axis.title.y = element_text(size = 50), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6b-2.pdf",width = 16,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, CCND2,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60),axis.title.y = element_text(size = 50), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6b-3.pdf",width = 16,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, BAX,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60),axis.title.y = element_text(size = 50), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()


pdf("Figure6b-4.pdf",width = 16,height = 10)
ggplot(Epithelial.cells.all.data, aes(ident, BAX,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_text(angle = 60, hjust = 1, size = 50),axis.text.y = element_text(size = 60),axis.title.y = element_text(size = 50), legend.key.size  = unit(1,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

}


## Figure 6c --- Comparing the overlapping gene set
{
   # cancer cell markers
  cancer.cell.markers <- FindMarkers(dt,ident.1 = c("Cancer cell"),only.pos = T)
  # GSE29272
  {
  GSE29272 <- read.table("GSE29272_DEGs.txt",sep = "\t",stringsAsFactors = F,header = T)
  GSE29272 <- GSE29272[-which(duplicated(GSE29272$Gene.symbol)),]
  rownames(GSE29272) <- as.character(GSE29272$Gene.symbol)
  GSE29272.genelist <- rownames(GSE29272)[which((GSE29272$adj.P.Val < 0.05) & (GSE29272$logFC) > log2(1.5))]
  } 
  length(intersect(GSE29272.genelist,rownames(cancer.cell.markers)))
}


## Figure 6d

load("Medoids.clusters.Rdata")

M <- cor(medoids.cluster)
# M1 <- dist(t(medoids.cluster))
# M <- as.data.frame(as.matrix(M1))
rownames(M) <- colnames(M) <- as.character(unique(dt@ident))
network.nodes <- as.data.frame(rownames(M))
network.nodes[,1] <- as.character(network.nodes$`rownames(M)`)
rownames(network.nodes) <- rownames(M)
network.node.num <- table(dt@ident)
network.nodes$num <- as.vector(as.numeric(network.node.num[rownames(network.nodes)]))

network.edges <- as.data.frame(matrix(0,nrow = nrow(M)^2,ncol = 3))
for(i in 1:nrow(M)){
  cluster_i <- rownames(M)[i]
  for (j in 1:ncol(M)){
   cluster_j <- rownames(M)[j] 
   row_k <- (i-1)*nrow(M) + j
    network.edges[row_k,] <- c(cluster_i,cluster_j,M[i,j])
  }
}

network.edges <- network.edges[-which(duplicated(network.edges$V3)),]


write.table(network.edges,file = "cellcluster.network.edge.Eu.txt",sep = "\t",quote = F)
write.table(network.nodes,file = "cellcluster.network.node.txt",sep = "\t",quote = F)



## Figure 6e
library(ggplot2)

medoids.cluster <- data.frame(medoids.cluster)
colnames(medoids.cluster) <- rownames(M)

pdf("Figure6e-1.pdf",width = 5,height = 5)
colnames(medoids.cluster)[13] <- "Cancer.cell"
colnames(medoids.cluster)[11] <- "Stem.like.cell"
par(mar = c(2,3,4,4))
p <- ggplot(medoids.cluster, aes(Stem.like.cell,Cancer.cell)) + geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE) + theme(
                 panel.background = element_blank(),
                 axis.text.x = element_text(size = 30),
                 axis.text.y = element_text(size = 30),
                 axis.title.x = element_text(size = 30),
                 axis.title.y = element_text(size = 30))
p
dev.off()

str(summary(lm(Cancer.cell~Stem.like.cell,data = medoids.cluster)))

# R^2 = 0.69

pdf("Figure6e-2.pdf",width = 5,height = 5)
colnames(medoids.cluster)[13] <- "Cancer.cell"
colnames(medoids.cluster)[11] <- "Enterocyte"
par(mar = c(2,3,4,4))
p <- ggplot(medoids.cluster, aes(Enterocyte,Cancer.cell)) + geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE) + theme(        
                 panel.background = element_blank(),
                 axis.text.x = element_text(size = 30),
                 axis.text.y = element_text(size = 30),
                 axis.title.x = element_text(size = 30),
                 axis.title.y = element_text(size = 0))
p
dev.off()

str(summary(lm(Cancer.cell~Enterocyte,data = medoids.cluster)))

# R^2 = 0.68


## Figure 5f

load("data.temp.res.2.annotation.Rdata")
{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}
dt.temp <- SubsetData(dt,ident.use = c("PMC","GMC","Enterocyte","Cancer cell","Stem-like cell"),do.scale = T,do.center = T)
refactor<-function(name,from,to){
  output<-as.character(name)
  output[which(output %in% from)]<-to
  output<-as.factor(output)
  names(output)<-names(name)
  return (output)
}
# dt.temp@ident <- refactor(dt.temp@ident, "NA?", "Enterocytes")
dt.temp@ident <- refactor(dt.temp@ident, c("PMC","GMC"), "Gastric mucous cells")
dt.temp@ident <- refactor(dt.temp@ident, c("Cancer cell"), "Cancer cells")


## Figure 6f
{

classical.markers.3 <- c("FABP1","CDH17","GRN")

Epithelial.cells.all.data <- (data.frame(as.matrix(dt.temp@data[classical.markers.3,])))
Epithelial.cells.all.data <- t(Epithelial.cells.all.data)
Epithelial.cells.all.data <- data.frame(Epithelial.cells.all.data)
Epithelial.cells.all.data$ident <- as.vector(as.character(dt.temp@ident))

pdf("Figure6f-1.pdf",width = 6,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, FABP1,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6f-2.pdf",width = 6,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, CDH17,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6f-3.pdf",width = 6,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, GRN,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf(file = "Figure6f--legend.pdf",width = 14,height = 6)
ggplot(Epithelial.cells.all.data, aes(ident, GRN,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "bottom",legend.key.width = unit(2.2,"cm"),axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0),legend.text = element_text(size=20),legend.spacing.x = unit(2,"lines"))
dev.off()

}




## Figure 6g
{
library(pheatmap)
load("data.temp.res.2.annotation.Rdata")

{
current.cluster.ids <- c("NA")
new.cluster.ids <- c("Fibroblast")
dt@ident <- plyr::mapvalues(x = dt@ident, from = current.cluster.ids, to = new.cluster.ids)
}

cell.included <- c("Cancer cell", "Enterocyte", "GMC","Goblet cell","PC","PMC","Stem-like cell","Enteroendocrine")
data.temp.used <- SubsetData(dt,ident.use = cell.included,do.scale = T,do.center = T)


EGC.genes.sepcific <- FindMarkers(dt,ident.1 = "Cancer cell",ident.2 = "Enterocyte",only.pos = T)
EGC.genes.sepcific.2 <- FindMarkers(dt,ident.1 = "Cancer cell",ident.2 = "Stem-like cell",only.pos = T)
EGC.genes.sepcific.3 <- FindMarkers(dt,ident.1 = "Cancer cell",ident.2 = "PMC",only.pos = T)
EGC.marker.list <- FindMarkers(dt,ident.1 = "Cancer cell",only.pos = T)


EGC.marker.list.list <- rownames(EGC.marker.list)[which((EGC.marker.list$avg_logFC > 0) & (EGC.marker.list$p_val_adj < 0.01) & (EGC.marker.list$pct.2<0.05)& (EGC.marker.list$pct.1 > 0.25))]

EGC.genes.sepcific.list <- rownames(EGC.genes.sepcific)[which((EGC.genes.sepcific$avg_logFC > 0) & (EGC.genes.sepcific$p_val_adj < 0.01))]
EGC.genes.sepcific.list.2 <- rownames(EGC.genes.sepcific.2)[which((EGC.genes.sepcific.2$avg_logFC > 0) & (EGC.genes.sepcific.2$p_val_adj < 0.01))]
EGC.genes.sepcific.list.3 <- rownames(EGC.genes.sepcific.3)[which((EGC.genes.sepcific.3$avg_logFC > 0) & (EGC.genes.sepcific.3$p_val_adj < 0.01))]

EGC.genes.specific.list.all <- intersect(intersect(EGC.genes.sepcific.list,EGC.genes.sepcific.list.2),EGC.genes.sepcific.list.3)


EGC.genes.sepcific.list <- intersect(EGC.genes.specific.list.all,EGC.marker.list.list) # up-regualted in EGC and down-regulated in enterocytes and other cell lineages
EGC.genes.sepcific.list <- intersect(EGC.marker.list.list,EGC.genes.specific.list.all)

classical.markers.3 <- EGC.genes.sepcific.list

Epithelial.cells.all.data <- (data.frame(as.matrix(data.temp.used@data[classical.markers.3,])))
Epithelial.cells.all.data <- t(Epithelial.cells.all.data)
Epithelial.cells.all.data <- data.frame(Epithelial.cells.all.data)
Epithelial.cells.all.data$ident <- as.vector(as.character(data.temp.used@ident))



pdf("Figure6g-1.pdf",width = 8,height = 4)
ggplot(Epithelial.cells.all.data, aes(ident, KLK10,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6g-2.pdf",width = 8,height = 4)
ggplot(Epithelial.cells.all.data, aes(ident, SLC11A2,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6g-3.pdf",width = 8,height = 4)
ggplot(Epithelial.cells.all.data, aes(ident,SULT2B1,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()


pdf("Figure6g-4.pdf",width = 8,height = 4)
ggplot(Epithelial.cells.all.data, aes(ident,KLK7,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6g-5.pdf",width = 8,height = 4)
ggplot(Epithelial.cells.all.data, aes(ident,KRT16,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(size = 60,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6g-6.pdf",width = 12,height = 10)
ggplot(Epithelial.cells.all.data, aes(ident, LMTK3,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_text(size = 60,angle = 60,hjust = 1),axis.title.x = element_blank(),axis.text.y = element_text(size = 90,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

pdf("Figure6g-7.pdf",width = 12,height = 10)
ggplot(Epithelial.cells.all.data, aes(ident, ECM1,fill = factor(ident),las = 2)) + geom_violin(scale = "width") + theme(legend.position = "none",axis.text.x = element_text(size = 60,angle = 60, hjust = 1),axis.title.x = element_blank(),axis.text.y = element_text(size = 90,face = "plain"),axis.title.y = element_text(size = 0), legend.key.size  = unit(1.5,"lines"),legend.text = element_text(size=25),legend.spacing.y = unit(2,"lines"))
dev.off()

rm(data.temp)
## Supplementary figure 6h

{

library(plotly)

label = c("NAG", "CAG","IMW","IMS","EGC")
dt@meta.data$batch <- gsub("CAN2","EGC",dt@meta.data$batch)
dt@meta.data$batch[grep("IMW",dt@meta.data$batch)] <- "IMW"
dt@meta.data$batch[grep("IMS",dt@meta.data$batch)] <- "IMS"
dt@meta.data$batch[grep("CAG",dt@meta.data$batch)] <- "CAG"
dt@meta.data$batch[grep("NAG",dt@meta.data$batch)] <- "NAG"

cancer.cell.info <- dt@meta.data[which(dt@ident %in% c("Cancer cell")),]
stem.cell.info <- dt@meta.data[which(dt@ident %in% c("Stem-like cell")),]
cancer.cell.samples <-  table(cancer.cell.info$batch)[label]
# cancer.cell.samples <- c(0,cancer.cell.samples)
stem.cell.samples <-  table(stem.cell.info$batch)[label]

col.total <- rainbow(16)


p <- plot_ly(
    type = "sankey",
    orientation = "h",
    width = 300,
    sizes = c(50,10),
    node = list(
      label = c("NAG", "CAG","IMW","IMS","EGC","Stem-like cell","Cancer cell"),
      color = c(rainbow(40)[39], "green", col.total[c(1:7,11:15)]),
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      )
    ),

    link = list(
      source = c(0:5,0:5),
      target = c(rep(6,6),rep(7,6)),
      value =  c(stem.cell.samples,cancer.cell.samples)
    )
  ) %>%
  layout(
    title = NULL)



tiff(file = "Supplementary_Figure6h.tiff", width = 200, height = 200, res = 200, unit = "mm")
p
dev.off()

}

}


## Supplementary figure 6i
valid.cells <- dt@cell.names[which(dt@meta.data$batch %in% c("IMS3","IMS4","CAN2"))]
Epithelial.cells.IIC <- SubsetData(dt,cell.use = valid.cells,do.center = T,do.scale = T)
Epithelial.cells.IIC@meta.data$orig.ident <- Epithelial.cells.IIC@ident
Epithelial.cells.IIC <- FindVariableGenes(object = Epithelial.cells.IIC, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.1, x.high.cutoff = 5, y.cutoff = 0.5,do.plot = T)
length(x = Epithelial.cells.IIC@var.genes)
Epithelial.cells.IIC <- RunPCA(object = Epithelial.cells.IIC,pc.genes = Epithelial.cells.IIC@var.genes, do.print = TRUE, pcs.print = 1:10,genes.print = 20,pcs.compute = 40)
# Epithelial.cells.IIC <- FindClusters(object = Epithelial.cells.IIC, reduction.type = "pca", dims.use = 1:20,
#     resolution = c(2), print.output = 0, save.SNN = TRUE)
Epithelial.cells.IIC <- RunTSNE(object = Epithelial.cells.IIC, dims.use = 1:20, do.fast = TRUE)

pdf(file = "Supplementary_Figure6i.pdf", width = 8, height = 8)
gg <- TSNEPlot(object = Epithelial.cells.IIC,do.label = T,label.size = 30,pt.size = 4,do.return = T) 
gg <- gg + theme(legend.text = element_text(size = 45),
                 legend.position = 'bottom',
                 legend.key.height = unit(1.4,"cm"),
                 legend.key.width = unit(1.4,"cm"),
                 axis.text.x = element_text(size = 35),
                 axis.text.y = element_text(size = 35),
                 axis.title.x = element_text(size = 35),
                 axis.title.y = element_text(size = 35))
gg
dev.off()  
```

### end


