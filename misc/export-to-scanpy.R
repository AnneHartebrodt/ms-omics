setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)
require(Seurat)
require(data.table)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(tidyr)
require(colorspace)
library(circlize)
library(wesanderson)
require(RColorBrewer)
library(SeuratDisk)

base.dir<-"/work/AnneSophiaHartebrodt#6698//ms-project/submission/results-final/scanpy"
filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698//ms-project/submission/results-final/ms-rna/filtering/')
dir.create(base.dir)
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

SaveH5Seurat(adata, filename = file.path(base.dir, "msrna.h5Seurat"))
Convert(file.path(base.dir, "msrna.h5Seurat"), dest = "h5ad")



filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698//ms-project/submission/results-final/ms-atac/')
atac<- readRDS(file.path(filtering.dir, "integrated.rds"))

SaveH5Seurat(atac, filename = file.path(base.dir, "msatac.h5Seurat"))
Convert(file.path(base.dir, "msatac.h5Seurat"), dest = "h5ad")


