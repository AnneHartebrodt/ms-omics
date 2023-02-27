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


filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

state.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/cellstates/')
dir.create(state.dir)


cellstates.astro<-list(c('PRKG1', 'AC012405.1','DPP10','AC073941.1','LSAMP','LINC00609', 'FBXL7'),
                 c('MT-CO3','MT-ND4','PLP1','MT-CO2','MT-CO1','IL1RAPL1','MT-ND1'),
                 c('AQP4','CLU','GJA1','PSAP','CPE','CST3','FTL'),
                 c('NLGN1','LINC01088','FAT3','KIAA1217','LRRTM4','PRUNE2','SORCS3'),
                 c('CSFAP299','SPAG17','ADGB','CFAP43','DNAH11','DNAH9','ZBBX'))
names(cellstates.astro)<-c('state0', 'state1', 'state2', 'state3', 'state4')


cellstates.immune<-list(c('RBFOX1', 'AC105402.3', 'CSMD1', 'MEG3', 'ROBO2', 'FAM155A', 'MTRNR2L8', 'OPCML', 'NRG1', 'RIMS2'),
           c('NAV2', 'DOCK8', 'SRGAP2', 'FRMD4A', 'APBB1IP', 'HS3ST4', 'SYNDIG1', 'ST6GAL1', 'CX3CR1', 'P2RY12'),
             c('CD163', 'ARGHGAP15', 'SRGAP1', 'F13A1', 'CD74', 'LRRK2', 'TRPS1', 'C1QA', 'DSE', 'UTRN'),
             c('FTL', 'APOE', 'CPM', 'GPNMB', 'ASAH1', 'C1QB', 'NHSL1', 'NLP', 'RPS19', 'CTSD'),
             c('IL1RAPL1', 'MAGI2', 'PCDH9', 'CTNNA3', 'DLG2', 'PTPRD', 'RNF220', 'ST18', 'SLC44A1', 'MAP7')
             )
names(cellstates.immune)<-c('state0', 'state1', 'state2', 'state3', 'state4')

cellstates.astro.new<-list(c('PRKG1', 'AC012405.1','DPP10','AC073941.1','LSAMP','LINC00609', 'FBXL7'),
                           c('MT-CO3','MT-ND4','PLP1','MT-CO2','MT-CO1','IL1RAPL1','MT-ND1'),
                           c('AQP4','CLU','GJA1','PSAP','CPE','CST3','FTL'),
                           c('NLGN1','LINC01088','FAT3','KIAA1217','LRRTM4','PRUNE2','SORCS3'),
                           c('CSFAP299','SPAG17','ADGB','CFAP43','DNAH11','DNAH9','ZBBX'),
                           c("CHI3L1", "MT1", "MT2", "M3", "SOD2", "VEGFA", "HIF1A", "CD44", "C3" ))
names(cellstates.astro.new)<-c('state0', 'state1', 'state2', 'state3', 'state4', 'state5')

cellstates.immune.new<-list(c('RBFOX1', 'AC105402.3', 'CSMD1', 'MEG3', 'ROBO2', 'FAM155A', 'MTRNR2L8', 'OPCML', 'NRG1', 'RIMS2'),
                            c('NAV2', 'DOCK8', 'SRGAP2', 'FRMD4A', 'APBB1IP', 'HS3ST4', 'SYNDIG1', 'ST6GAL1', 'CX3CR1', 'P2RY12'),
                            c('CD163', 'ARGHGAP15', 'SRGAP1', 'F13A1', 'CD74', 'LRRK2', 'TRPS1', 'C1QA', 'DSE', 'UTRN'),
                            c('FTL', 'APOE', 'CPM', 'GPNMB', 'ASAH1', 'C1QB', 'NHSL1', 'NLP', 'RPS19', 'CTSD'),
                            c('IL1RAPL1', 'MAGI2', 'PCDH9', 'CTNNA3', 'DLG2', 'PTPRD', 'RNF220', 'ST18', 'SLC44A1', 'MAP7'),
                            c("SPP1", "FTH1", "REL", "HIF1A", "C3" )
)
names(cellstates.immune.new)<-c('state0', 'state1', 'state2', 'state3', 'state4', 'state5')



cellstate_list<-list()
cellstate_list[['Astrocytes']]<-cellstates.astro
cellstate_list[['Immune cells']]<-cellstates.immune

cellstate_list_new<-list()
cellstate_list_new[['Astrocytes']]<-cellstates.astro.new
cellstate_list_new[['Immune cells']]<-cellstates.immune.new

csl<-list()
csl[['lerma']]<-cellstate_list
csl[['extra_genes']]<-cellstate_list_new

#for (ln in c('lerma', 'extra_genes')){
ct = 'Astrocytes'
ln = 'lerma'
print(ln)
cellstate_list<- csl[[ln]]
ct.dir<-file.path(state.dir, ct)
dir.create(ct.dir)
#for (vf in c(2000, 3000)){
vf <- 2000
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))
adata$batch<-sapply(adata$batch, function(x) gsub('smpD', 'MS9_NAWM', x))
adata$lesion_type<-sapply(adata$batch, function(x) str_split(x, '_')[[1]][2])
adata$patient<-sapply(adata$batch, function(x) str_split(x, '_')[[1]][1])

adata@active.ident<-adata@meta.data$celltypes_annotated
Idents(adata)<-"celltypes_annotated"
adata@active.assay<-'RNA'

cellstates<- cellstate_list[[ct]]
l<-lapply(cellstates, function(x) length(x))
l<-lapply(names(l), function(x) rep(x, l[x][[1]]))
cellstate.df<-data.table(genes = unlist(cellstates), state=unlist(l))
fwrite(cellstate.df, file=file.path(ct.dir, paste0(ln,'_', ct, '_cellstates.tsv')), sep='\t')

# Subset astrocytes
Idents(adata)<-"celltypes_annotated"

# Remove weirdo from the plot

astrocytes<-subset(adata, subset = celltypes_annotated == ct )

#REcompute clusters

astrocytes<-NormalizeData(astrocytes)
astrocytes <- ScaleData(astrocytes,scale.max = 10)
astrocytes <- RunPCA(astrocytes, npcs = 30)
astrocytes <- RunUMAP(astrocytes, reduction = "pca", dims = 1:30)
astrocytes <- FindNeighbors(astrocytes, reduction = "pca", dims = 1:30)
astrocytes <- FindClusters(astrocytes)

colooo <- rainbow(length(unique(astrocytes$seurat_clusters)))
names(colooo)<-unique( astrocytes$seurat_clusters)

p1<-DimPlot(astrocytes, group.by = 'seurat_clusters', cols = colooo)
p2 <-DimPlot(astrocytes, group.by = 'batch')
p4<-DimPlot(astrocytes, group.by = 'lesion_type')

p<-p1 | p2 | p4
p
ggsave(p, file = file.path(ct.dir, paste0(ln, vf, 'umaps.pdf')), width = 30, height = 10, units = 'cm')


sub<-subset(astrocytes, features = unlist(cellstates))
mat<-sub@assays$RNA@scale.data

lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057')
celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4])
sample_color = as.character(brewer.pal( 8,'Paired'))
sample_color = c('C1'=sample_color[1], 'C4'=sample_color[2], 'C3'=sample_color[8], 'MS1'=sample_color[3], 'MS3'=sample_color[4], 'MS6'=sample_color[5], 'MS7'=sample_color[6], 'MS9'=sample_color[7])
sex_color = c('male'='#3D4849', 'female'='#D7D9D9')

column_ha = HeatmapAnnotation(
  Clusters = sub$seurat_clusters, Batch = sub$batch, Lesion=sub$lesion_type, Patient=sub$patient,
  col=list(Clusters=colooo, Lesion=lesion_color, Patient=sample_color))

hr = hclust(dist(mat), method = "average")
clusters = dendextend::cutree(hr, k = 5)
cludi<-lapply(unique(clusters), function(x) names(clusters[clusters==x]))
names(cludi)<-paste0('State' , 1:5)


hm<-Heatmap(mat, show_column_names  = F, 
        top_annotation = column_ha,
         column_km  = 5,
        split = clusters) 
pdf(file = file.path(ct.dir, paste0(ln, vf, 'heatmap_clustered.pdf')), width = 15, height = 10)
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
names<-rownames(mat)[unlist(row_order(ht))]
df<-data.table(names, cluster=row_col)

fwrite(df, file = file.path(ct.dir, paste0(ln, vf, 'row_labels.tsv')), sep='\t')
fwrite(annot, file = file.path(ct.dir, paste0(ln, vf, 'column_labels.tsv')), sep='\t')

astrocytes$hier.cluster<-as.numeric(astrocytes$hier.cluster)
astrocytes$hier.cluster<-as.factor(astrocytes$hier.cluster)
Idents(astrocytes)<-'hier.cluster'
p3<-DotPlot(
  astrocytes,
  cols = c('white', 'blue'),
  assay = 'RNA',
  features =cludi,
  #group.by = 'hier.cluster',
  cluster.idents = FALSE,
  scale = TRUE,
  scale.max = NA,
)

p3<-p3+theme(axis.text.x=element_text(size=18, angle=90,hjust = 1),axis.text.y=element_text(size=18), legend.box = 'horizontal') +ylab('Cluster')
ggsave(p3, file=file.path(ct.dir, paste0(ln, vf, 'dotplot.pdf')), width = 40, height = 10, units = 'cm')


dimplot<-DimPlot(astrocytes, group.by = 'hier.cluster')
dimplot
ggsave(dimplot, filename = file.path( ct.dir, paste0(ln, vf, 'umap_with_hierarchical_labelling.pdf')))

 
allumaps<-p1|p2|dimplot                  
ggsave(allumaps, filename = file.path(ct.dir, paste0(ln, vf, '_umap_all.pdf')), width = 30, height = 10, units = 'cm')




ct <- 'Immune cells'
ln <- 'lerma'
  print(ln)
  cellstate_list<- csl[[ln]]
  ct.dir<-file.path(state.dir, ct)
  dir.create(ct.dir)
  #for (vf in c(2000, 3000)){
vf <- 2000
    adata<- readRDS(file.path(filtering.dir, "annotated.rds"))
    adata$batch<-sapply(adata$batch, function(x) gsub('smpD', 'MS9_NAWM', x))
    adata$lesion_type<-sapply(adata$batch, function(x) str_split(x, '_')[[1]][2])
    adata$patient<-sapply(adata$batch, function(x) str_split(x, '_')[[1]][1])
    
    adata@active.ident<-adata@meta.data$celltypes_annotated
    Idents(adata)<-"celltypes_annotated"
    adata@active.assay<-'RNA'
    
    cellstates<- cellstate_list[[ct]]
    l<-lapply(cellstates, function(x) length(x))
    l<-lapply(names(l), function(x) rep(x, l[x][[1]]))
    cellstate.df<-data.table(genes = unlist(cellstates), state=unlist(l))
    fwrite(cellstate.df, file=file.path(ct.dir, paste0(ln,'_', ct, '_cellstates.tsv')), sep='\t')
    
    # Subset astrocytes
    Idents(adata)<-"celltypes_annotated"
    
    # Remove weirdo from the plot
    
    astrocytes<-subset(adata, subset = celltypes_annotated == ct )
    
    #REcompute clusters
    
    astrocytes<-NormalizeData(astrocytes)
    if (vf==3000){
      astrocytes<-FindVariableFeatures(astrocytes, nfeatures = vf)
    }
    astrocytes <- ScaleData(astrocytes,scale.max = 10)
    astrocytes <- RunPCA(astrocytes, npcs = 30)
    astrocytes <- RunUMAP(astrocytes, reduction = "pca", dims = 1:30)
    astrocytes <- FindNeighbors(astrocytes, reduction = "pca", dims = 1:30)
    astrocytes <- FindClusters(astrocytes)
    
    colooo <- rainbow(length(unique(astrocytes$seurat_clusters)))
    names(colooo)<-unique( astrocytes$seurat_clusters)
    
    p1<-DimPlot(astrocytes, group.by = 'seurat_clusters', cols = colooo)
    p2 <-DimPlot(astrocytes, group.by = 'batch')
    p4<-DimPlot(astrocytes, group.by = 'lesion_type')
    
    p<-p1 | p2 | p4
    p
    ggsave(p, file = file.path(ct.dir, paste0(ln, vf, 'umaps.pdf')), width = 30, height = 10, units = 'cm')
    
   

    sub<-subset(astrocytes, features = unlist(cellstates))
    mat<-sub@assays$RNA@scale.data
    
    lesion_color  = c('NAWM' = '#00798c', 'AL'='#d1495b', 'CA'='#edae49', 'WM'='#66a182', 'RL'='#2e4057')
    celltype_color = as.character(wes_palette('GrandBudapest1', 4))
    celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4])
    sample_color = as.character(brewer.pal( 8,'Paired'))
    sample_color = c('C1'=sample_color[1], 'C4'=sample_color[2], 'C3'=sample_color[8], 'MS1'=sample_color[3], 'MS3'=sample_color[4], 'MS6'=sample_color[5], 'MS7'=sample_color[6], 'MS9'=sample_color[7])
    sex_color = c('male'='#3D4849', 'female'='#D7D9D9')
    
    column_ha = HeatmapAnnotation(
      Clusters = sub$seurat_clusters, Batch = sub$batch, Lesion=sub$lesion_type, Patient=sub$patient,
      col=list(Clusters=colooo, Lesion=lesion_color, Patient=sample_color))
    
    hr = hclust(dist(mat), method = "average")
    clusters = dendextend::cutree(hr, k = 5)
    cludi<-lapply(unique(clusters), function(x) names(clusters[clusters==x]))
    names(cludi)<-paste0('State' , 1:5)
    hm<-Heatmap(mat, show_column_names  = F, 
                top_annotation = column_ha,
                column_km  = 5,
                split = clusters) 
    pdf(file = file.path(ct.dir, paste0(ln, vf, 'heatmap_clustered.pdf')), width = 15, height = 10)
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
    names<-rownames(mat)[unlist(row_order(ht))]
    df<-data.table(names, cluster=row_col)
    
    fwrite(df, file = file.path(ct.dir, paste0(ln, vf, 'row_labels.tsv')), sep='\t')
    
    
    fwrite(annot, file = file.path(ct.dir, paste0(ln, vf, 'column_labels.tsv')), sep='\t')
    
    astrocytes$hier.cluster<-as.numeric(astrocytes$hier.cluster)
    astrocytes$hier.cluster<-as.factor(astrocytes$hier.cluster)
    Idents(astrocytes)<-'hier.cluster'
    p3<-DotPlot(
      astrocytes,
      cols = c('white', 'blue'),
      assay = 'RNA',
      features =cludi,
      #group.by = 'hier.cluster',
      cluster.idents = FALSE,
      scale = TRUE,
      scale.max = NA,
    )
    
    p3<-p3+theme(axis.text.x=element_text(size=18, angle=90,hjust = 1),axis.text.y=element_text(size=18), legend.box = 'horizontal') +ylab('Cluster')
    ggsave(p3, file=file.path(ct.dir, paste0(ln, vf, 'dotplot.pdf')), width = 40, height = 10, units = 'cm')
    
    
    
    dimplot<-DimPlot(astrocytes, group.by = 'hier.cluster')
    dimplot
    ggsave(dimplot, filename = file.path( ct.dir, paste0(ln, vf, 'umap_with_hierarchical_labelling.pdf')))
    
    
    allumaps<-p1|p2|dimplot                  
    ggsave(allumaps, filename = file.path(ct.dir, paste0(ln, vf, '_umap_all.pdf')), width = 30, height = 10, units = 'cm')
    
  


Idents(adata)<-'celltypes_annotated'
VlnPlot(adata, features = 'LRP2', idents = 'Oligodendrocytes', group.by = 'ms')
VlnPlot(adata, features = 'LRP2', idents = 'Oligodendrocytes', group.by = 'batch')

