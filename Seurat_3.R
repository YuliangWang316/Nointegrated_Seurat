library(Seurat)
library(cowplot)
KO.data <- Read10X(data.dir = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Treg_KO/outs/filtered_feature_bc_matrix")
WT.data <- Read10X(data.dir = "D:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/WT/outs/filtered_feature_bc_matrix")

KO.data <- as.data.frame(KO.data)
WT.data <- as.data.frame(WT.data)

for (i in 1:4105) {
   colnames(KO.data)[i] <- paste(colnames(KO.data)[i],"KO",i,sep = "-")  
}

for (i in 1:3905) {
  colnames(WT.data)[i] <- paste(colnames(WT.data)[i],"WT",i,sep = "-")  
}

KO.metadata<-data.frame(colnames(KO.data),rep("KO",4105))
WT.metadata<-data.frame(colnames(WT.data),rep("WT",3905))
colnames(KO.metadata)<-c("barcode","group")
colnames(WT.metadata)<-c("barcode","group")
pbmc.metadata<-rbind(WT.metadata,KO.metadata)
rownames(pbmc.metadata)<-pbmc.metadata[,1]
pbmc.data<-cbind(WT.data,KO.data)

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",meta.data = pbmc.metadata,min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 0.8)
pbmc <- RunUMAP(pbmc, dims = 1:20)
pbmc <- RunTSNE(pbmc, dims = 1:20)
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "umap",split.by = "group")
DimPlot(pbmc, reduction = "tsne")
DimPlot(pbmc, reduction = "tsne",split.by = "group")

Idents(pbmc)<-pbmc@meta.data$group
DimPlot(pbmc, reduction = "umap")
DimPlot(pbmc, reduction = "tsne")



