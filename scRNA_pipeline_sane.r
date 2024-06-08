library(Seurat)
library(tidyverse)
library(future)
library(ggplot2)

#Read data matrix generated with RSEM
countMat=read.table('invitro_diff_scRNA_terminals_transcripts_reheader.mat',row.names=1, header=T)

#Create Seurat object
hEDT <- CreateSeuratObject(counts = countMat, min.cells = 3, min.features = 200, project = "hEDT", names.delim = "_", names.field = 1)
head(rownames(GetAssayData(hEDT)))

#Read mitochondrial and ribosomal transcripts
mitoread=read.table('mito_transcripts_list',row.names=1)
mitotrans=rownames(mitoread)
riboread=read.table('ribo_transcripts_list',row.names=1)
ribotrans=rownames(riboread) #,rownames(GetAssayData(hEDT)))
used_ribotrans=intersect(ribotrans,rownames(GetAssayData(hEDT)))


#Mark and check mitochondrial and ribosomal genes
hEDT[["percent.mito"]] <- PercentageFeatureSet(object = hEDT, features = mitotrans)
hEDT[["percent.ribo"]] <- PercentageFeatureSet(object = hEDT, features = used_ribotrans)


#Filter out abnormal mitochondrial or ribosomal gene expressing transcripts
hEDT <- subset(x = hEDT, subset = percent.mito < 25 & percent.ribo < 40 )

#Plot some filters
pdf('hEDT_nCount_RNA_post.pdf')
VlnPlot(object = hEDT, features = "nCount_RNA")
dev.off()

pdf('hEDT_nFeature_RNA_post.pdf')
VlnPlot(object = hEDT, features = "nFeature_RNA")
dev.off()


#Data normalization
hEDT <- NormalizeData(hEDT, normalization.method = "LogNormalize", scale.factor = 10000)
hEDT <- FindVariableFeatures(hEDT, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(hEDT)
hEDT <- ScaleData(hEDT, features = all.genes, vars.to.regress = c("nCount_RNA", "percent.mito", "percent.ribo")) # This is painfully slow. Consider using [plan("multisession", workers = 16)]


#saveRDS
saveRDS(hEDT, "scRNA/hEDT.rds")
#saveRDS
save(hEDT, file = "hEDT_filtered.RData")

#Run PCA
hEDT <- RunPCA (hEDT, features = VariableFeatures(object = hEDT), ndims.print = 1:2)
hEDT_pca12 <- hEDT

#Save and load later
save(hEDT_pca12, "hEDT_terminals_pca12.Robj")
hEDT<- load("hEDT_terminals_pca12.Robj", verbose=TRUE)

#Plot elbow plot
pdf('elbow_result.pdf', width=2.5, height=2.5)
ElbowPlot(object = hEDT, ndims = 50) 
dev.off()


#Run UMAP clustering with different dimensions
hEDT_umapd6 <- RunUMAP(object = hEDT_pca12, dims = 1:6)
hEDT_umapd10 <- RunUMAP(object = hEDT_pca12, dims = 1:10)
hEDT_umapd15 <- RunUMAP(object = hEDT_pca12, dims = 1:15)
hEDT_umapd20 <- RunUMAP(object = hEDT_pca12, dims = 1:20)



#Make dimplot
pdf('hEDT_umap_dimplotumapd6.pdf',width=4, height=4)
DimPlot(object = hEDT_umapd6, label = T, label.size = 3) + NoLegend()
dev.off()

pdf('hEDT_umap_dimplotumapd10.pdf',width=4, height=4)
DimPlot(object = hEDT_umapd10, label = T, label.size = 3) + NoLegend()
dev.off()


pdf('hEDT_umap_dimplotumapd15.pdf',width=4, height=4)
DimPlot(object = hEDT_umapd15, label = T, label.size = 3) + NoLegend()
dev.off()


pdf('hEDT_umap_dimplotumapd20.pdf',width=4, height=4)
DimPlot(object = hEDT_umapd20, label = T, label.size = 3) + NoLegend()
dev.off()


#Write data
write.table(hEDT_umapd6[["umap"]]@cell.embeddings, file = "hEDT_umapd6.tsv", sep="\t", col.names=NA, quote =F)





#Uniformly expressed
unifo=c('hEDT-00083255', 'hEDT-00038452', 'hEDT-00038886', 'hEDT-00015330', 'hEDT-00095535', 'hEDT-00061863', 'hEDT-00091739', 'hEDT-00014474', 'hEDT-00091870', 'hEDT-00046359')
pdf('hEDT_umap_dimplot6_uniform.pdf')
FeaturePlot(object = hEDT_umap6, features = unifo)
dev.off()

uniform1=c('hEDT-00083255','hEDT-00038452','hEDT-00029760','hEDT-00091739')
psc1=c('hEDT-00090876','hEDT-00065605','hEDT-00065702','hEDT-00024447')
endoderm1=c('hEDT-00003306','hEDT-00092327','hEDT-00088530','hEDT-00003308')
ectoderm1=c('hEDT-00027673','hEDT-00014647','hEDT-00083803','hEDT-00024275')
mesoderm1=c('hEDT-00074096','hEDT-00021389','hEDT-00052392','hEDT-00051834')
some_meso1=c('hEDT-00020555','hEDT-00075372','hEDT-00088817','hEDT-00089410')

uniform2=c('hEDT-00091870','hEDT-00091739','hEDT-00000934','hEDT-00002904')
psc2=c('hEDT-00065606','hEDT-00021227','hEDT-00058276','hEDT-00052011')
endoderm2=c('hEDT-00087437','hEDT-00066533','hEDT-00039117','hEDT-00079488')
ectoderm2=c('hEDT-00014648','hEDT-00069377','hEDT-00087183','hEDT-00069377')
mesoderm2=c('hEDT-00067907','hEDT-00024290','hEDT-00033940','hEDT-00055643')


ectoderm3=c('hEDT-00023130','hEDT-00053366','hEDT-00089867','hEDT-00072618')
ectoderm4=c('hEDT-00055883','hEDT-00072539','hEDT-00023130','hEDT-00087183')

#Uniformly expressed
pdf('hEDT_umap_dimplot6_uniform1.pdf')
FeaturePlot(object = hEDT_umapd6, features = uniform1)
dev.off()

pdf('hEDT_umap_dimplot6_uniform2.pdf')
FeaturePlot(object = hEDT_umapd6, features = uniform2)
dev.off()


pdf('hEDT_umap_dimplot6_PSC1.pdf')
FeaturePlot(object = hEDT_umapd6, features = psc1, order = T)
dev.off()

pdf('hEDT_umap_dimplot6_PSC2.pdf')
FeaturePlot(object = hEDT_umapd6, features = psc2, order = T)
dev.off()


pdf('hEDT_umap_dimplot6_endoderm1.pdf')
FeaturePlot(object = hEDT_umapd6, features = endoderm1, order = T)
dev.off()

pdf('hEDT_umap_dimplot6_endoderm2.pdf')
FeaturePlot(object = hEDT_umapd6, features = endoderm2, order = T)
dev.off()


pdf('hEDT_umap_dimplot6_ectoderm1.pdf')
FeaturePlot(object = hEDT_umapd6, features = ectoderm1, order = T)
dev.off()

pdf('hEDT_umap_dimplot6_ectoderm2.pdf')
FeaturePlot(object = hEDT_umapd6, features = ectoderm2, order = T)
dev.off()


pdf('hEDT_umap_dimplot6_ectoderm3.pdf')
FeaturePlot(object = hEDT_umapd6, features = ectoderm3, order = T)
dev.off()

pdf('hEDT_umap_dimplot6_ectoderm4.pdf')
FeaturePlot(object = hEDT_umapd6, features = ectoderm4, order = T)
dev.off()



pdf('hEDT_umap_dimplot6_mesoderm1.pdf')
FeaturePlot(object = hEDT_umapd6, features = mesoderm1, order = T)
dev.off()

pdf('hEDT_umap_dimplot6_mesoderm2.pdf')
FeaturePlot(object = hEDT_umapd6, features = mesoderm2, order = T)
dev.off()


pdf('hEDT_umap_dimplot6_some_meso.pdf')
FeaturePlot(object = hEDT_umapd6, features = some_meso, order = T)
dev.off()


pdf('hEDT_umap6_dotplot_mesoderm.pdf')
DotPlot(object = hEDT_umap, features = psc) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()




merged=c('hEDT-00083255', 'hEDT-00038452', 'hEDT-00091739', 'hEDT-00014474', 'hEDT-00091870', 'hEDT-00046359', 'hEDT-00091395', 'hEDT-00029760', 'hEDT-00000934', 'hEDT-00002904', 'hEDT-00024447', 'hEDT-00024356', 'hEDT-00074481', 'hEDT-00052011', 'hEDT-00065702', 'hEDT-00058276', 'hEDT-00019994', 'hEDT-00090876', 'hEDT-00065605', 'hEDT-00065606', 'hEDT-00098203', 'hEDT-00003306', 'hEDT-00092327', 'hEDT-00088530', 'hEDT-00003308', 'hEDT-00087437', 'hEDT-00092328', 'hEDT-00039117', 'hEDT-00075288', 'hEDT-00079488', 'hEDT-00075281', 'hEDT-00074096', 'hEDT-00021389', 'hEDT-00052392', 'hEDT-00051834', 'hEDT-00092760', 'hEDT-00074095', 'hEDT-00071997', 'hEDT-00083429', 'hEDT-00089132', 'hEDT-00001816', 'hEDT-00073863', 'hEDT-00001817', 'hEDT-00009683', 'hEDT-00089133', 'hEDT-00070884', 'hEDT-00087654', 'hEDT-00027673', 'hEDT-00014647', 'hEDT-00083803', 'hEDT-00024275', 'hEDT-00014648', 'hEDT-00053366', 'hEDT-00055588', 'hEDT-00087181', 'hEDT-00043906', 'hEDT-00057030', 'hEDT-00087183', 'hEDT-00089867', 'hEDT-00018842', 'hEDT-00043029', 'hEDT-00072618')
pdf('hEDT_umap6_dotplot_merged.pdf')
DotPlot(object = hEDT_umap, features = merged) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf('hEDT_umap6_dotplot_merged_flipped.pdf')
DotPlot(object = hEDT_umap, features = merged) + coord_flip() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


