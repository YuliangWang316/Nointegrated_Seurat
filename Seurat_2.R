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
pbmc.markers<-FindAllMarkers(pbmc,logfc.threshold = 0,min.pct = 0,only.pos = TRUE)
write.table(pbmc.markers,file = "c:/Users/xjmik/Desktop/pbmc.markers.txt",sep = "\t")
features<-c("Pdcd1","Nrp1","Tnfrsf18","Icos","Ctla4","Lag3","Tnfrsf4","Havcr2","Entpd1","Il2ra")#Nt5e(CD73) delete

library(RColorBrewer)
library(ggplot2)
DoHeatmap(pbmc, features = features, size = 4,
          angle = 90) + scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(10))
DotPlot(pbmc, features = features) + RotatedAxis()
VlnPlot(pbmc,features = "Ifng",sort = TRUE)+stat_compare_means()
VlnPlot(pbmc,features = "Tbx21")
VlnPlot(pbmc,features = "Tigit")

library(ggplot2)
VlnPlot(pbmc,features = "Ctla4",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Lag3",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Pdcd1",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Havcr2",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Tnfrsf18",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Nrp1",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Icos",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))
VlnPlot(pbmc,features = "Tnfrsf4",pt.size = 0) + scale_fill_manual(values = c("#00BFC4","#F8766D"))

scale_fill_manual()
IL2_stat5<-read.table("g:/GSE98638/IL2_STAT5_symbol.txt",sep = "\t",header = TRUE)
IL6_stat3<-read.table("g:/GSE98638/IL6_STAT3_symbol.txt",sep = "\t",header = TRUE)

IL2_stat5<-IL2_stat5[2:200,]
IL6_stat3<-IL6_stat3[2:88,]

pbmc <- CellCycleScoring(pbmc, s.features = IL2_stat5, g2m.features = IL6_stat3, set.ident = FALSE)
IL2_stat5_list<-list(IL2_stat5)
IL6_stat3_list<-list(IL6_stat3)
pbmc <- AddModuleScore(pbmc,features = IL2_stat5_list,name = "IL2_stat5")
pbmc <- AddModuleScore(pbmc,features = IL6_stat3_list,name = "IL6_stat3")

A<-Idents(pbmc)
Idents(pbmc)<-pbmc$group
VlnPlot(pbmc,features = c("IL2_stat51","IL6_stat31","Jmjd1c","Nrp1","Pdcd1","Ifng"),pt.size = 0)
FeaturePlot(pbmc,features = c("IL2_stat51","IL6_stat31","Jmjd1c","Nrp1","Pdcd1","Ifng"))


#bulk
pbmc@meta.data$cell.cluster<-Idents(pbmc)
Idents(pbmc)<-pbmc@meta.data$group
library(ggplot2)
library(ggpubr)
VlnPlot(pbmc,features = "Nrp1") +stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Pdcd1")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Icos")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Tnfrsf4")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Havcr2")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Lag3")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Il2ra")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Ctla4")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Tnfrsf18")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Cd44")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Entpd1")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Nt5e")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Nt5e")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Scd1")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Fasn")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Acaca")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Acly")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Srebf1")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Srebf2")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Cd36")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Socs1")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Socs3")+stat_compare_means(method = "t.test",label.x.npc = "center")
VlnPlot(pbmc,features = "Irf1")+stat_compare_means(method = "t.test",label.x.npc = "center")

pbmc.markers<-FindMarkers(pbmc,ident.1 = "KO",ident.2 = "WT",logfc.threshold = 0,min.pct = 0)
write.table(pbmc.markers,file = "d:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/total.markers_KO_vs_WT.txt",sep = "\t")


gene_list=pbmc.markers[,c("avg_logFC","p_val_adj")]
#gene_list=log10(gene_list)
#gene_list[,"log2FoldChange"]=log2(gene_list[,"log2FoldChange"])
colnames(gene_list)=c("logFC","padj")
gene_list$change = ifelse(gene_list$padj < 0.000001 & abs(gene_list$logFC) >= 0.1, 
                          ifelse(gene_list$logFC> 0.1 ,'Up 512 genes','Down 573 genes'),
                          'no_change')
#colored_point<-gene_list[gene_list$threshold == "TRUE",]
library("ggplot2")
pdf("d:/Jmjd1c_Treg_Tumor/single_cell_RNAseq/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/volcanoplot.pdf")
g = ggplot(data=gene_list, aes(x=logFC, y=-log10(padj),color=change)) + geom_point(alpha=0.4, size=1.75)  + xlim(c(-2, 3)) + ylim(c(0, 310)) +xlab("log2 fold change ko_vs_wt") + ylab("-log10 p-value") + scale_color_manual(values=c("#191970","#000000","#B0171F")) + theme_set(theme_bw()) + theme(panel.grid.major=element_line(colour=NA))
#+geom_text(mapping=aes(label=rownames(colored_point)),data = colored_point,check_overlap = TRUE)
print(g)
dev.off()

library(dplyr)
library(Seurat)
library(patchwork)

top100 <- pbmc.markers %>% top_n(n = 100, wt = avg_logFC)
DoHeatmap(pbmc, features = top100$gene) #+ NoLegend()



avg.pbmc<-AverageExpression(pbmc)
write.table(avg.pbmc,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/avg.txt",sep = "\t")

library(dplyr)
Foxp3<-filter(as.data.frame(t(as.data.frame(pbmc@assays[["RNA"]]@scale.data))),Foxp3>0)
Foxp3_barcode<-rownames(Foxp3)

#1
Foxp3_t<-as.data.frame(t(Foxp3))
Foxp3_t_KO.data<-select(Foxp3_t,contains("KO"))
FOxp3_t_WT.data<-select(Foxp3_t,contains("WT"))
Foxp3_KO.metadata<-data.frame(colnames(Foxp3_t_KO.data),rep("KO",1095))
Foxp3_WT.metadata<-data.frame(colnames(FOxp3_t_WT.data),rep("WT",1073))
colnames(Foxp3_KO.metadata)<-c("barcode","group")
colnames(Foxp3_WT.metadata)<-c("barcode","group")
Foxp3.metadata<-rbind(Foxp3_KO.metadata,Foxp3_WT.metadata)
rownames(Foxp3.metadata)<-Foxp3.metadata[,1]
Foxp3.data<-cbind(Foxp3_t_KO.data,FOxp3_t_WT.data)#combine

#A directly combine ko wt matrix
Foxp3_Seurat <- CreateSeuratObject(counts = Foxp3.data, project = "pbmc3k",meta.data = Foxp3.metadata,min.cells = 3, min.features = 200)
Foxp3_Seurat[["percent.mt"]] <- PercentageFeatureSet(Foxp3_Seurat, pattern = "^mt-")
VlnPlot(Foxp3_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Foxp3_Seurat <- subset(Foxp3_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)#
#Foxp3_Seurat <-  NormalizeData(Foxp3_Seurat, normalization.method = "RC")
#Foxp3_Seurat <- FindVariableFeatures(Foxp3_Seurat, selection.method = "vst", nfeatures = 2000)
all.genes_2 <- rownames(Foxp3_Seurat)
Foxp3_Seurat <- ScaleData(Foxp3_Seurat, features = all.genes_2)
Foxp3_Seurat <- RunPCA(Foxp3_Seurat,features =all.genes_2)
ElbowPlot(Foxp3_Seurat)
Foxp3_Seurat <- FindNeighbors(Foxp3_Seurat, dims = 1:20)
Foxp3_Seurat <- FindClusters(Foxp3_Seurat, resolution = 0.8)
Foxp3_Seurat <- RunUMAP(Foxp3_Seurat, dims = 1:20)
Foxp3_Seurat <- RunTSNE(Foxp3_Seurat, dims = 1:20)
DimPlot(Foxp3_Seurat, reduction = "umap")
DimPlot(Foxp3_Seurat, reduction = "umap",split.by = "group")
DimPlot(Foxp3_Seurat, reduction = "tsne")
DimPlot(Foxp3_Seurat, reduction = "tsne",split.by = "group")
#bulk
Foxp3_Seurat@meta.data$cell.cluster<-Idents(Foxp3_Seurat)
Idents(Foxp3_Seurat)<-Foxp3_Seurat@meta.data$group
Foxp3_Seurat.markers<-FindMarkers(Foxp3_Seurat,ident.1 = "WT",ident.2 = "KO",logfc.threshold = 0,min.pct = 0)
write.table(Foxp3_Seurat.markers,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/Foxp3_directlycombine_wt_komatrix.txt",sep = "\t")
avg.Foxp3_Seurat<-AverageExpression(Foxp3_Seurat)
write.table(avg.Foxp3_Seurat,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/avg.Foxp3_Seurat.txt",sep = "\t")
#B no combine in fact as same as A 
#Foxp3_t_Seurat <- CreateSeuratObject(counts = Foxp3_t, project = "pbmc3k",meta.data = Foxp3.metadata,min.cells = 3, min.features = 200)
#Foxp3_t_Seurat[["percent.mt"]] <- PercentageFeatureSet(Foxp3_t_Seurat, pattern = "^mt-")
#VlnPlot(Foxp3_t_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Foxp3_t_Seurat <- subset(Foxp3_t_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
#Foxp3_t_Seurat <- NormalizeData(Foxp3_t_Seurat, normalization.method = "RC")
#Foxp3_t_Seurat <- FindVariableFeatures(Foxp3_t_Seurat, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(Foxp3_t_Seurat)
#Foxp3_t_Seurat <- ScaleData(Foxp3_t_Seurat, features = all.genes)
#Foxp3_t_Seurat <- RunPCA(Foxp3_t_Seurat, features = VariableFeatures(object = Foxp3_t_Seurat))
#ElbowPlot(Foxp3_t_Seurat)
#Foxp3_t_Seurat <- FindNeighbors(Foxp3_t_Seurat, dims = 1:20)
#Foxp3_t_Seurat <- FindClusters(Foxp3_t_Seurat, resolution = 0.8)
#Foxp3_t_Seurat <- RunUMAP(Foxp3_t_Seurat, dims = 1:20)
#Foxp3_t_Seurat <- RunTSNE(Foxp3_t_Seurat, dims = 1:20)
#DimPlot(Foxp3_t_Seurat, reduction = "umap")
#DimPlot(Foxp3_t_Seurat, reduction = "umap",split.by = "group")
#DimPlot(Foxp3_t_Seurat, reduction = "tsne")
#DimPlot(Foxp3_t_Seurat, reduction = "tsne",split.by = "group")
#bulk
#Foxp3_t_Seurat@meta.data$cell.cluster<-Idents(Foxp3_t_Seurat)
#Idents(Foxp3_t_Seurat)<-Foxp3_t_Seurat@meta.data$group
#Foxp3_t_Seurat.markers<-FindMarkers(Foxp3_t_Seurat,ident.1 = "WT",ident.2 = "KO",logfc.threshold = 0,min.pct = 0)
#write.table(Foxp3_t_Seurat.markers,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/Foxp3_t_Seurat_no_combine_directuse.txt",sep = "\t")
#avg.Foxp3_t_Seurat<-AverageExpression(Foxp3_t_Seurat)
#write.table(avg.Foxp3_t_Seurat,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/avg.Foxp3_t_Seurat.txt",sep = "\t")
#2
Foxp3_Sub<-SubsetData(pbmc,cells = Foxp3_barcode)
Foxp3_Sub[["percent.mt"]] <- PercentageFeatureSet(Foxp3_Sub, pattern = "^mt-")
VlnPlot(Foxp3_Sub, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Foxp3_Sub <- subset(Foxp3_Sub, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
Foxp3_Sub <- NormalizeData(Foxp3_Sub)#, normalization.method = "RC")
Foxp3_Sub <- FindVariableFeatures(Foxp3_Sub, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Foxp3_Sub)
Foxp3_Sub <- ScaleData(Foxp3_Sub, features = all.genes)
Foxp3_Sub <- RunPCA(Foxp3_Sub, features = VariableFeatures(object = Foxp3_Sub))
ElbowPlot(Foxp3_Sub)
Foxp3_Sub <- FindNeighbors(Foxp3_Sub, dims = 1:20)
Foxp3_Sub <- FindClusters(Foxp3_Sub, resolution = 0.8)
Foxp3_Sub <- RunUMAP(Foxp3_Sub, dims = 1:20)
Foxp3_Sub <- RunTSNE(Foxp3_Sub, dims = 1:20)
DimPlot(Foxp3_Sub, reduction = "umap")
DimPlot(Foxp3_Sub, reduction = "umap",split.by = "group")
DimPlot(Foxp3_Sub, reduction = "tsne")
DimPlot(Foxp3_Sub, reduction = "tsne",split.by = "group")
FeaturePlot(Foxp3_Sub,features = c("Xist","Jmjd1c"),split.by = "group")
#bulk
Foxp3_Sub@meta.data$cell.cluster<-Idents(Foxp3_Sub)
Idents(Foxp3_Sub)<-Foxp3_Sub@meta.data$group
Foxp3_Sub.markers<-FindMarkers(Foxp3_Sub,ident.1 = "WT",ident.2 = "KO",logfc.threshold = 0,min.pct = 0)
write.table(Foxp3_Sub.markers,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/Foxp3_Sub.markers.txt",sep = "\t")
avg.Foxp3_Sub<-AverageExpression(Foxp3_Sub)
write.table(avg.Foxp3_Sub,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/avg.Foxp3_Sub.txt",sep = "\t")

#Foxp3_Sub_RC<-SubsetData(pbmc,cells = Foxp3_barcode)
#Foxp3_Sub_RC[["percent.mt"]] <- PercentageFeatureSet(Foxp3_Sub_RC, pattern = "^mt-")
#VlnPlot(Foxp3_Sub_RC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#Foxp3_Sub_RC <- subset(Foxp3_Sub_RC, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 10)
#Foxp3_Sub_RC <- NormalizeData(Foxp3_Sub_RC,normalization.method = "RC")
#Foxp3_Sub_RC <- FindVariableFeatures(Foxp3_Sub_RC, selection.method = "vst", nfeatures = 2000)
#all.genes <- rownames(Foxp3_Sub_RC)
#Foxp3_Sub_RC <- ScaleData(Foxp3_Sub_RC, features = all.genes)
#Foxp3_Sub_RC <- RunPCA(Foxp3_Sub_RC, features = VariableFeatures(object = Foxp3_Sub_RC))
#ElbowPlot(Foxp3_Sub_RC)
#Foxp3_Sub_RC <- FindNeighbors(Foxp3_Sub_RC, dims = 1:20)
#Foxp3_Sub_RC <- FindClusters(Foxp3_Sub_RC, resolution = 0.8)
#Foxp3_Sub_RC <- RunUMAP(Foxp3_Sub_RC, dims = 1:20)
#Foxp3_Sub_RC <- RunTSNE(Foxp3_Sub_RC, dims = 1:20)
#DimPlot(Foxp3_Sub_RC, reduction = "umap")
#DimPlot(Foxp3_Sub_RC, reduction = "umap",split.by = "group")
#DimPlot(Foxp3_Sub_RC, reduction = "tsne")
#DimPlot(Foxp3_Sub_RC, reduction = "tsne",split.by = "group")

#Foxp3_Sub_RC@meta.data$cell.cluster<-Idents(Foxp3_Sub_RC)
#Idents(Foxp3_Sub_RC)<-Foxp3_Sub_RC@meta.data$group
#Foxp3_Sub_RC.markers<-FindMarkers(Foxp3_Sub_RC,ident.1 = "WT",ident.2 = "KO",logfc.threshold = 0,min.pct = 0)
#write.table(Foxp3_Sub_RC.markers,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/Foxp3_Sub_RC.txt",sep = "\t")
#avg.Foxp3_Sub_RC<-AverageExpression(Foxp3_Sub_RC)
#write.table(avg.Foxp3_Sub_RC,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/avg.Foxp3_Sub_RC.txt",sep = "\t")

a<-as.data.frame(pbmc@assays[["RNA"]]@scale.data)
library(dplyr)
a.ko.data<-select(a,contains("KO"))
a.wt.data<-select(a,contains("WT"))
K<-rep("NA",16956)
K<-as.data.frame(K)
B<-cbind(K,a)
D<-toupper(rownames(B))
E<-cbind(D,B)
colnames(E)[1]<-"Gene_Symbol"
colnames(E)[2]<-"Gene_Description"
F<-colnames(E)
G<-rbind(F,E)

write.table(G,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/pbmc_scaledata.txt",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)



c<-c(rep("KO",4037),rep("WT",3850))
c<-as.data.frame(c)
c<-t(c)
write.table(c,file = "d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/sample.txt",sep = "\t")

FeaturePlot(pbmc,features = "Foxp3")
VlnPlot(pbmc,features="Nrp1")




z<-as.data.frame(pbmc@assays[["RNA"]]@data)

VlnPlot(pbmc,features = "Cd4")
write.table(z,file="d:/201692_20210126_1_foxp3jmjd1c_mc205_tumour_wt_ko_treg/Seurat2_result/z.txt",sep="\t")
