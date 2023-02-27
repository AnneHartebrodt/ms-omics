setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(harmony)
require(tidyverse)
require(data.table)
library(cowplot)
library(rsvd)
library(Rtsne)
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
require(biovizBase)
require(wesanderson)
require(gtable)

library(magrittr)
library(liana)
require(SeuratObject)
require(stringr)
require(ComplexHeatmap)
require(stringr)
require(RColorBrewer)
require(colorspace)


filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

state.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/')


adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

adata@active.ident<-adata@meta.data$celltypes_annotated
Idents(adata)<-"celltypes_annotated"
adata@active.assay<-'RNA'

    
# Subset astrocytes
Idents(adata)<-"celltypes_annotated"


for(ct in c( 'Oligodendrocytes', 'Immune cells', 'OPCs', 'Astrocytes')){
# Remove weirdo from the plot
astrocytes<-subset(adata, subset = celltypes_annotated == ct )
astrocytes<-NormalizeData(astrocytes)
astrocytes <- ScaleData(astrocytes,scale.max = 10)
astrocytes <- RunPCA(astrocytes, npcs = 30)
astrocytes <- RunUMAP(astrocytes, reduction = "pca", dims = 1:30)
astrocytes <- FindNeighbors(astrocytes, reduction = "pca", dims = 1:30)
astrocytes <- FindClusters(astrocytes)
    

lesion.de.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/', ct)
de.genes<-fread(file.path(lesion.de.dir, paste0('DE_genes_', ct, '.tsv')), sep='\t', header=TRUE)


sub<-subset(astrocytes, features = unlist(de.genes$gene))
sub <- subset(sub, subset = nFeature_RNA > 0 )
mat<-sub@assays$RNA@scale.data

lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057')
color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057')
celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4])
sample_color = as.character(brewer.pal( 8,'Paired'))
sample_color = c('C1'=sample_color[1], 'C4'=sample_color[2], 'C3'=sample_color[8], 'MS1'=sample_color[3], 'MS3'=sample_color[4], 'MS6'=sample_color[5], 'MS7'=sample_color[6], 'MS9'=sample_color[7])
sex_color = c('male'='#3D4849', 'female'='#D7D9D9')

column_ha = HeatmapAnnotation(
  MS = sub$ms, Batch = sub$batch, Lesion=sub$lesion_type, Patient=sub$patient,
  col=list( Lesion=lesion_color, Patient=sample_color))

color<-divergingx_hcl(n = 300, palette = 'RdBu')[30:270]
q<-quantile(as.vector(mat), c(0.01, 0.99))

hr = hclust(dist(mat), method = "average")
clusters = dendextend::cutree(hr, k = 5)
hm<-Heatmap(mat, 
            show_column_names  = F, 
            top_annotation = column_ha,
            column_km  = 5,
            split = clusters,
            raster_by_magick = TRUE,
            
            col = circlize::colorRamp2(c(-2, 0, 6), colors= c('Red' , 'white' ,'Darkblue')),
) 

state.dir.ct<-file.path(file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/', ct))
dir.create(state.dir.ct)
pdf(file = file.path(state.dir.ct, paste0('all_edgR_heatmap_clustered.pdf')), width = 15, height = 10)
ht = draw(hm)
dev.off()

length_col<-lapply(column_order(ht), function(x) length(x))
length_col<-lapply(names(length_col), function(x) rep(x, length_col[x][[1]]))
annot<-data.table(rownames = colnames(mat)[unlist(column_order(ht))], hier.cluster = unlist(length_col))

astrocytes@meta.data$rownames<-rownames(astrocytes@meta.data)
astrocytes@meta.data<-merge(astrocytes@meta.data, annot, by='rownames')
rownames(astrocytes@meta.data)<-astrocytes@meta.data$rownames
Idents(astrocytes)<-'hier.cluster'

row_col<-lapply(row_order(ht), function(x) length(x))
row_col<-unlist(lapply(names(row_col), function(x) rep(x, row_col[[x]])))
names<-de.genes[unlist(row_order(ht))]$gene
df<-data.table(names, cluster=row_col)
fwrite(df, file=file.path(state.dir.ct, paste0('clustering_labels.tsv')), sep='\t')
}

