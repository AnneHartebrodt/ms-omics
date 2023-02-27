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
require(colorspace)

library(magrittr)
library(liana)
require(SeuratObject)
require(stringr)
require(ComplexHeatmap)
#require(leiden)
require(RColorBrewer)

filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

state.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/cellstates/')
dir.create(state.dir)


adata<- readRDS(file.path(filtering.dir, "annotated.rds"))
ct <- 'Oligodendrocytes'

ct.dir<-file.path(state.dir, ct)
dir.create(ct.dir)

adata@active.ident<-adata@meta.data$celltypes_annotated
Idents(adata)<-"celltypes_annotated"
adata@active.assay<-'RNA'

# Subset astrocytes
Idents(adata)<-"celltypes_annotated"
astrocytes<-subset(adata, subset = celltypes_annotated == ct)

#REcompute clusters
# astrocytes<-NormalizeData(astrocytes)
# astrocytes <- ScaleData(astrocytes,scale.max = 10)
# astrocytes <- RunPCA(astrocytes, npcs = 30)
# astrocytes<-RunHarmony(astrocytes, group.by.vars = 'batch')
# astrocytes <- RunUMAP(astrocytes, reduction = "harmony", dims = 1:30)
# astrocytes <- FindNeighbors(astrocytes, reduction = "harmony", dims = 1:30)
# astrocytes <- FindClusters(astrocytes, resolution = 0.5)

colooo <- rainbow(length(unique(astrocytes$seurat_clusters)))
names(colooo)<-unique( astrocytes$seurat_clusters)
# 
# p1<-DimPlot(astrocytes, group.by = 'seurat_clusters', cols = colooo)
# p2 <-DimPlot(astrocytes, group.by = 'batch')
# p4<-DimPlot(astrocytes, group.by = 'lesion_type')
# p<-p1 | p2 | p4
# p
# ggsave(p, file = file.path(ct.dir, paste0( 'umaps_harmony.pdf')), width = 20, height = 10, units = 'cm')

Idents(astrocytes)<-'seurat_clusters'
markers<-FindAllMarkers(astrocytes, group.by = 'seurat_clusters')
markers<-as.data.table(markers)
top<-unique(markers[p_val_adj<0.05 & (pct.1 >0.20 | pct.2>0.20) & abs(avg_log2FC)>1]$gene)
# 
# saveRDS(astrocytes, file=file.path(ct.dir, paste0(ct, '.rds')), compress = F)
fwrite(markers, file=file.path(ct.dir, paste0('markers.csv')), sep='\t')
# 

#de_oligo<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/celltypes/Oligodendrocytes_DE_genes.tsv')

lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057')
celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4])
sample_color = as.character(brewer.pal( 8,'Paired'))
sample_color = c('C1'=sample_color[1], 'C4'=sample_color[2], 'C3'=sample_color[8], 'MS1'=sample_color[3], 'MS3'=sample_color[4], 'MS6'=sample_color[5], 'MS7'=sample_color[6], 'MS9'=sample_color[7])
sex_color = c('male'='#3D4849', 'female'='#D7D9D9')

column_ha = HeatmapAnnotation(
  Clusters = astrocytes$seurat_clusters, Batch = astrocytes$batch, Lesion=astrocytes$lesion_type, Patient=astrocytes$patient,
  col=list(Clusters=colooo, Lesion=lesion_color, Patient=sample_color))


sub<-subset(astrocytes, features = top)
mat<-sub@assays$RNA@scale.data


Idents(sub)<-'seurat_clusters'
dp<-DotPlot(astrocytes, features = top, dot.min = 0.2)+theme(axis.text.x = element_text(angle=90, hjust = 1))

dp

hr = hclust(dist(mat), method = 'average')
k= 3
clusters = dendextend::cutree(hr, k = k)
cludi<-lapply(unique(clusters), function(x) names(clusters[clusters==x]))
names(cludi)<-paste0('State' , c(2,3,1))


# hm<-Heatmap(mat, show_column_names  = F, 
#             top_annotation = column_ha,
#             split = clusters,
#             raster_by_magick = TRUE,
#             column_split = sub$seurat_clusters,
#             heatmap_legend_param = list(
#               color_bar = "continuous"
#             ))
# 
# pdf(file = file.path(ct.dir, paste0( 'heatmap_clustered_seurat.pdf')), width = 15, height = 10)
# ht = draw(hm)
# dev.off()

k=3
pdf(file = file.path(ct.dir, paste0( 'heatmap_clustered_hierarchical.pdf')), width = 15, height = 10)
hm<-Heatmap(mat, 
            show_column_names  = F, 
            top_annotation = column_ha,
            column_km  = k,
            split = clusters,
            raster_by_magick = TRUE)
ht = draw(hm)
dev.off()

length_col<-lapply(column_order(ht), function(x) length(x))
length_col<-lapply(names(length_col), function(x) rep(x, length_col[x][[1]]))
annot<-data.table(rownames = colnames(mat)[unlist(column_order(ht))], hier.cluster = unlist(length_col))

astrocytes@meta.data$rownames<-rownames(astrocytes@meta.data)
astrocytes@meta.data<-merge(astrocytes@meta.data, annot, by='rownames')
rownames(astrocytes@meta.data)<-astrocytes@meta.data$rownames
Idents(astrocytes)<-'hier.cluster.y'

row_col<-lapply(row_order(ht), function(x) length(x))
row_col<-unlist(lapply(names(row_col), function(x) rep(x, row_col[[x]])))
names<-rownames(mat)[unlist(row_order(ht))]
df<-data.table(names, cluster=row_col)



dimplot<-DimPlot(astrocytes, group.by = 'hier.cluster')
dimplot 
ggsave(dimplot, filename = file.path( ct.dir, paste0( 'umap_with_hierarchical_labelling.pdf')))
ggsave(DimPlot(astrocytes, group.by = 'seurat_clusters'), filename = file.path( ct.dir, paste0(  'umap_with_hierarchical_labelling.pdf')))


Idents(astrocytes)<-'hier.cluster'
astrocytes$hier.cluster<-as.numeric(astrocytes$hier.cluster)
astrocytes$hier.cluster<-as.factor(astrocytes$hier.cluster)
p3<-DotPlot(
  astrocytes,
  cols = c('white', 'blue'),
  assay = 'RNA',
  features =cludi[c(3,1,2)],
  #group.by = 'hier.cluster',
  idents = NULL,
  split.by = NULL,
  cluster.idents = FALSE,
  scale = TRUE,
  scale.by = "radius",
  scale.max = NA,
  #dot.min = 0.1
)  
p3<-p3+
  theme(axis.text.x=element_text(size=18, angle=90,hjust = 1),axis.text.y=element_text(size=18), legend.box = 'horizontal') + ylab('Cluster')


ggsave(p3, file=file.path(ct.dir, paste0('dotplot.pdf')), width = 40, height = 10, units = 'cm')

df$cluster<-recode(df$cluster, '3'='1', '1'='2', '2'='3')
fwrite(df, file = file.path(ct.dir, 'row_labels.tsv'), sep='\t')
fwrite(annot, file = file.path(ct.dir, 'column_labels.tsv'), sep='\t')




