setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(rlang)
require(Seurat)
library(anndata)
require(data.table)
require(SingleCellExperiment)
require(scds)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(tidyr)
require(devtools)
require(colorspace)
library(rlang)
require(ComplexHeatmap)
library(circlize)
library(wesanderson)
require(RColorBrewer)
require(magick)
require(ggvenn)

require(nVennR)
library(png)
library("jpeg")



out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna-spatial/umaps-top-de')
dir.create(out.dir, showWarnings = F, recursive = T)

de.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/plots')


filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering/')
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))


for(ds in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  de.dir.ds<- file.path(de.dir, ds)
  
  de_markers<-fread(file=file.path(de.dir.ds, paste0(ds, '_differential_genes.tsv')), sep='\t')
  
  top.de <-de_markers %>% group_by(cluster) %>% top_n(10, wt = abs(avg_log2FC))
  top.de<-as.data.table(top.de)
  f<- intersect(unique(top.de$gene) ,adata@assays$RNA@var.features)
  p5<-FeaturePlot(adata, features =f) 
  
  ggsave(p5, file=file.path(out.dir, paste0(ds, '_umap_top_de_genes_from_spatial.pdf')), width=25, height = 25)
  
}

