setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
#renv::init('final-r-env')
renv::activate('final-r-env')
renv::repair()
set.seed(11)

require(Seurat)
require(data.table)
require(SingleCellExperiment)
require(scds)
require(ggplot2)
require(dplyr)
require(tidyr)
require(devtools)
require(wesanderson)

require(biomaRt)
require('harmony')
require(DropletUtils)

base.dir<- file.path('//work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-rna/')


res.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results/ms-rna-with-schirmer')
dir.create(res.dir, recursive = T)

celltype.dir<-file.path(res.dir, 'celltypes')
dir.create(celltype.dir, recursive = T)

filtering.dir<-file.path(res.dir, 'filtering')

plot.dir<- file.path(res.dir, 'filtering/plots')
dir.create(plot.dir, recursive = T)



### STEP 1: read and combine the data sets
adatalist<-list()
cells<-list()
names<-c( 'C3_WM', 'C4_WM', 'MS1_AL','MS1_NAWM','MS3_CA','MS3_RL','MS7_AL','MS7_CA','MS7_RL', 'smpD', 'MS6_RL', 'MS6_CA', 'MS6_NAWM', 'C1_WM')
sex <-c('male', 'female', 'male', 'male', 'female', 'female', 'male', 'male', 'male', 'female', 'male', 'male', 'male', 'male')

for (fname in names){
  C3WM<-Read10X_h5(file.path(base.dir, fname, "/outs/raw_feature_bc_matrix.h5"))
  
  #create seurat object
  C3WM<-CreateSeuratObject(
    C3WM,
    project = "fname",
    assay = "RNA",
    min.cells = 15
  )
  
  C3WMss<-as.SingleCellExperiment(C3WM)
  e.out <- emptyDrops(counts(C3WMss))
  is.cell <- e.out$FDR <= 0.01

  C3WM[['is_cell']]<-is.cell
  C3WM <- subset(C3WM, subset = is_cell)
  
  C3WM[["percent.mt"]] <- PercentageFeatureSet(C3WM, pattern = "^MT-")
  C3WM <- subset(C3WM, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 1)
  
  ## Annotate doublets
  C3WMss<-as.SingleCellExperiment(C3WM)
  
  # run profile
  C3WMss = cxds(C3WMss,retRes = TRUE)
  C3WMss = bcds(C3WMss,retRes = TRUE,verb=TRUE)
  C3WMss = cxds_bcds_hybrid(C3WMss)
  
  C3WM[['hybrid']]<-C3WMss$hybrid_score
  
  adatalist[[fname]]<-C3WM
  cells[[fname]]<-C3WM@assays$RNA@counts@Dim[2]
}

schirmername<- c( 'SRR9123038', 'SRR9123039', 'SRR9123040', 'SRR9123041', 'SRR9123034', 'SRR9123035', 'SRR9123036', 'SRR9123037', 'SRR9123032')

for (fname in schirmername){
  print(fname)
  C3WM<-Read10X_h5(file.path('/work/AnneSophiaHartebrodt#6698/ms-project/ms-rna/data/schirmer/',paste0(fname, '_count'),'/outs/filtered_feature_bc_matrix.h5'))
  
  #create seurat object
  C3WM<-CreateSeuratObject(
    C3WM,
    project = "fname",
    assay = "RNA",
    min.cells = 15
  )
  
  # C3WMss<-as.SingleCellExperiment(C3WM, assay = 'RNA')
  # e.out <- emptyDrops(counts(C3WMss))
  # is.cell <- e.out$FDR <= 0.01

  # C3WM[['is_cell']]<-is.cell
  # C3WM <- subset(C3WM, subset = is_cell)
  # 
  C3WM[["percent.mt"]] <- PercentageFeatureSet(C3WM, pattern = "^MT-")
  C3WM <- subset(C3WM, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 1)
  
  ## Annotate doublets
  C3WMss<-as.SingleCellExperiment(C3WM)
  
  # run profile
  C3WMss = cxds(C3WMss,retRes = TRUE)
  C3WMss = bcds(C3WMss,retRes = TRUE,verb=TRUE)
  C3WMss = cxds_bcds_hybrid(C3WMss)
  
  C3WM[['hybrid']]<-C3WMss$hybrid_score
  
  adatalist[[fname]]<-C3WM
  cells[[fname]]<-C3WM@assays$RNA@counts@Dim[2]
  
}

schirmersex<-c('male','female','male','male','male','female','male','female','female')

adata<-merge(adatalist[['C3_WM']], adatalist[2:length(adatalist)], add.cell.ids =c(names, schirmername))
adata[['batch']]<-rep(c(names, schirmername), cells)
adata[['sex']] <- rep(c(sex, schirmersex), cells)

saveRDS(adata, file = file.path(filtering.dir, "before_filtering.rds"), compress=FALSE)
#adata<-readRDS(file.path(filtering.dir, "before_filtering.rds"))


## STEP 3: Normalize and Scale the data
adata <- NormalizeData(adata, normalization.method = "LogNormalize", scale.factor = 10000)
adata <- FindVariableFeatures(adata, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(adata)
adata <- ScaleData(adata, features = all.genes)
adata <- RunPCA(adata, features = VariableFeatures(object = adata))

# Integrate data and RUN harmony
adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.2)
adata <- RunUMAP(adata, dims = 1:10, n.components = 3)

adata<-RunHarmony(adata, group.by.vars = 'batch', reduction = 'pca' )
adata <- FindNeighbors(adata, reduction = "harmony", dims = 1:30) %>% FindClusters()
adata <- RunUMAP(adata, reduction = "harmony", dims = 1:30)

saveRDS(adata, file = file.path(filtering.dir, "transformed.rds"), compress=FALSE)
#adata<- readRDS(file.path(filtering.dir, "transformed.rds"))

### STEP 5: There is one weird cluster:
### What is going on?
adata[['is_doublet']]<- ifelse(adata$hybrid>quantile(adata$hybrid, .90)[1], 'doublet', 'singlet')
doublet1<-FeaturePlot(adata, features = 'hybrid')
doublet2<-DimPlot(adata, group.by = 'is_doublet')+ggtitle('Doublet predictions')
doublet1
doublet2

batch_before<-DimPlot(adata, group.by = 'batch')+ggtitle('Experimental Batch')
batch_before
ggsave(batch_before, file =file.path(plot.dir, paste0('batch_umap_before.pdf')) , width = 12, height = 10, units = 'cm')
ggsave(batch_before, file =file.path(plot.dir, paste0('batch_umap_before.png')) , width = 12, height = 10, units = 'cm')

feature_before<-DimPlot(adata, label = TRUE)+ggtitle('Clustering')
feature_before
ggsave(feature_before, file =file.path(plot.dir, paste0('feature_before_umap.pdf')) , width = 12, height = 10, units = 'cm')
ggsave(feature_before, file =file.path(plot.dir, paste0('feature_before_umap.png')) , width = 12, height = 10, units = 'cm')


ggsave(doublet1, file =file.path(plot.dir, paste0('scatter_hybrid_scores.pdf')), width = 12, height = 10, units = 'cm' )
ggsave(doublet1, file =file.path(plot.dir, paste0('scatter_hybrid_scores.png')), width = 12, height = 10, units = 'cm' )

ggsave(doublet2, file =file.path(plot.dir, paste0('binarized_doublet_scores.pdf')), width = 12, height = 10, units = 'cm' )
ggsave(doublet2, file =file.path(plot.dir, paste0('binarized_doublet_scores.png')) , width = 12, height = 10, units = 'cm')

# Ok this cluster is full of doublets and MT-RNA (debris)
DimPlot(adata, label = TRUE)
mtpercent.cluster.8 <- VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "seurat_clusters", pt.size = 0, combine = F)
mtpercent.cluster.8[[3]] | feature_before
ggsave(mtpercent.cluster.8[[3]], file = file.path(plot.dir, paste0('cluster13_is_suspicious.pdf')))
ggsave(mtpercent.cluster.8[[3]], file = file.path(plot.dir, paste0('cluster13_is_suspicious.png')))


# STEP 5.1. 
## This whole sample is frankly suspicious.
## STEP 6: remove the offending clusters
adata <- subset(adata, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 1)

adata<-subset(adata, subset = seurat_clusters !=13  & is_doublet == 'singlet')

saveRDS(adata, file = file.path(filtering.dir, 'filtered.rds'), compress=FALSE)
#adata<-readRDS(file.path(filtering.dir, 'filtered.rds'))

# STEP 5: Re-embed data
adata <- RunPCA(adata, features = VariableFeatures(object = adata))
adata <- FindNeighbors(adata, dims = 1:10)
adata <- FindClusters(adata, resolution = 0.2)
adata <- RunUMAP(adata, dims = 1:10, n.components = 2)
adata <-RunHarmony(adata, reduction = 'pca', group.by.vars = 'batch')
adata <- RunUMAP(adata, dims = 1:10, n.components = 2, reduction = 'harmony')

batch_after<-DimPlot(adata, group.by = 'batch')
batch_after
ggsave(batch_after, file =file.path(plot.dir, paste0('batch_umap_after.pdf')) , width = 12, height = 10, units = 'cm')
ggsave(batch_after, file =file.path(plot.dir, paste0('batch_umap_after.png')) , width = 12, height = 10, units = 'cm')

feature_after<-DimPlot(adata, label=TRUE)
feature_after
ggsave(feature_after, file =file.path(plot.dir, paste0('feature_after_umap.pdf')), width = 12, height = 10, units = 'cm')
ggsave(feature_after, file =file.path(plot.dir, paste0('feature_after_umap.png')), width = 12, height = 10, units = 'cm')


# STEP 7: make marker gene clusters
markers<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/annotation/marker_genes.txt')
astro<-FeaturePlot(adata, features = markers[Cell_type=='Astrocytes']$Gene_name)
astro
ggsave(astro, file =file.path(plot.dir, paste0('astrocytes_umap.pdf')) )
ggsave(astro, file =file.path(plot.dir, paste0('astrocytes_umap.png')) )

oligo<-FeaturePlot(adata, features = markers[Cell_type=='Oligodendrocytes']$Gene_name)
oligo
ggsave(oligo, file =file.path(plot.dir, paste0('oligodendrocytes_umap.pdf')) )
ggsave(oligo, file =file.path(plot.dir, paste0('oligodendrocytes_umap.png')) )

opcs<-FeaturePlot(adata, features = markers[Cell_type=='OPCs']$Gene_name)
opcs
ggsave(opcs, file =file.path(plot.dir, paste0('opcs_umap.png')) )

micro<-FeaturePlot(adata, features = markers[Cell_type=='Microglia-Macrophages']$Gene_name)
micro
ggsave(micro, file =file.path(plot.dir, paste0('microglia.png')) )

neurons<-FeaturePlot(adata, features = markers[Cell_type=='Neurons']$Gene_name)
neurons
ggsave(neurons, file =file.path(plot.dir, paste0('neurons.png')) )

bcell<-FeaturePlot(adata, features = markers[Cell_type=='B-cells']$Gene_name)
bcell
ggsave(bcell, file =file.path(plot.dir, paste0('bcells.png')) )

tcell<-FeaturePlot(adata, features = markers[Cell_type=='T-cells']$Gene_name)
tcell
ggsave(tcell, file =file.path(plot.dir, paste0('tcells.png')))

mrest<-FeaturePlot(adata, features = markers[Cell_type=='Microglia-resting']$Gene_name)
ggsave(mrest, file =file.path(plot.dir, paste0('mrest.png')))
mrest

endothelial<- FeaturePlot(adata, features = c('ERG', 'VWF', 'EMCN', 'CLDN5', 'EGFL7'), ncol = 3)
endothelial
ggsave(endothelial, file =file.path(plot.dir, paste0('endothelial.png')))
ggsave(endothelial, file =file.path(plot.dir, paste0('endothelial.pdf')))

peric<-FeaturePlot(adata, features = c('PDGFRB',  'NOTCH3', 'TBX18'))
peric
ggsave(peric, file =file.path(plot.dir, paste0('peric.png')))
ggsave(peric, file =file.path(plot.dir, paste0('peric.pdf')))

plot.dir.ind<-file.path(plot.dir, 'umaps')
dir.create(plot.dir.ind)
for(gene in c(markers$Gene_name,c('ERG', 'VWF', 'EMCN', 'CLDN5', 'EGFL7'))){
  tcell<-FeaturePlot(adata, features = gene)+NoLegend()
  ggsave(tcell, file =file.path(plot.dir.ind, paste0(gene, '.png')), width = 10, height = 10, units = 'cm')
  ggsave(tcell, file =file.path(plot.dir.ind, paste0(gene, '.pdf')), width = 10, height = 10, units = 'cm')
}




####
feature_after
coarse<-adata$seurat_clusters
coarse[coarse==1]<-0
coarse[coarse==6]<-0

coarse[coarse==7]<-3

coarse[coarse==10]<-2
coarse[coarse==11]<-2
coarse[coarse==8]<-2



coarse<-droplevels(coarse)

cellnames<-coarse
levels(cellnames)<-c('Oligodendrocytes', 'Astrocytes',  'Neurons','OPCs','Immune cells', 'Endothelial cells')
adata[['celltypes_annotated']]<-cellnames

celltype_color = c(as.character(wes_palette('GrandBudapest1', 4)), as.character(wes_palette('GrandBudapest2', 4)))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4],
                    'Neurons' = celltype_color[5], 'Endothelial cells' = celltype_color[6], 'Pericytes'=celltype_color[7], 'B-cells'=celltype_color[8])

cellt<- DimPlot(adata, group.by = 'celltypes_annotated', cols = celltype_color)+ggtitle('Annotation based on marker genes')
ggsave(cellt, file =file.path(plot.dir, paste0('annotated', '.png')), width = 13, height = 10, units = 'cm')
ggsave(cellt, file =file.path(plot.dir, paste0('annotated', '.pdf')), width = 13, height = 10, units = 'cm')

ggsave(cellt+NoLegend(), file =file.path(plot.dir, paste0('annotated_no_legend', '.png')), width = 10, height = 10, units = 'cm')
ggsave(cellt+NoLegend(), file =file.path(plot.dir, paste0('annotated_no_legend', '.pdf')), width = 10, height = 10, units = 'cm')


saveRDS(adata, file = file.path(filtering.dir, "annotated.rds"), compress=FALSE)
adata<-readRDS(file.path(filtering.dir, "annotated.rds"))

# Find marker genes for all groups
Idents(adata)<-'celltypes_annotated'
all.markers<- FindAllMarkers(adata, min.pct = 0.25, logfc.threshold = 0.25)
n<-rownames(all.markers)
all.markers<-as.data.table(all.markers)
all.markers$gene<-n
fwrite(all.markers, file=file.path(plot.dir, 'all_marker_genes.tsv'), sep='\t')

plot.dir<-
  all.markers$up_or_down<-ifelse(all.markers$avg_log2FC>0, 'UP', 'DOWN')
list_by_cluster<- all.markers[abs(avg_log2FC)>0.5 & p_val<0.05] %>% group_by(cluster) %>% group_split()

for (l in list_by_cluster){
  name <- l$cluster[1]
  fwrite(l, file.path(celltype.dir.dir, paste0(name, '_DE_genes.tsv')), sep='\t')
}

