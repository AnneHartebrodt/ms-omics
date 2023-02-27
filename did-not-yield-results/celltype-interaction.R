setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

library(harmony)
require(data.table)
require(SingleCellExperiment)
require(scDblFinder)
require(scds)
library(cowplot)
library(rsvd)
library(Rtsne)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
require(biovizBase)
library(JASPAR2020)

library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
require(wesanderson)
require("readxl") # CRAN version
#BiocManager::install("GSEABase")
require(topGO)
require(GSEABase)
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x 

#axis labels (ex: ridgeplot)
library(ggplot2)  
library(magrittr)
library(liana)
require(SeuratObject)
require(liana)
require(dplyr)

filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/')
dir.create(filtering.dir)

liana.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/liana')
dir.create(liana.dir)

adata<- readRDS(file.path(filtering.dir, "integrated.rds"))


adata@active.assay<-'RNA'
adata@active.ident<- adata$batch
Idents(adata)<-'batch'
liana_objects<-list()
for (b in unique(adata$batch)){
  subbi <- subset(adata, subset = batch==b)
  Idents(subbi)<-'predicted_celltype'
  liana_test_subbi <- liana_wrap(subbi)
  liana_objects[[b]]<-liana_test_subbi
  saveRDS(liana_test_subbi, file = file.path(liana.dir, paste0(b, '.Rds')))
  
}


for (b in unique(adata$batch)){
  lo<-liana_objects[[b]] %>%  liana_aggregate()
  fwrite(lo, file=file.path(liana.dir, paste0(b, '.tsv')), sep='\t')
}

#saveRDS(liana_objects, file = file.path(liana.dir, 'all.Rds'))

#liana_objects<-readRDS(file.path(liana.dir, 'all.Rds'))

liana_objects<-list()
for (b in unique(adata$batch)){
  li<-readRDS(file = file.path(liana.dir, paste0(b, '.Rds')))
  liana_objects[[b]]<-li
  
}

lian<-list()
for(b in unique(adata$batch)){
  liana_test_C4_WM <- liana_objects[[b]] %>%
    liana_aggregate()
  liana_test_C4_WM$batch<-b
  liana_test_C4_WM$lesion<-stringr::str_split(b, '_')[[1]][2]
  lian[[b]]<-as.data.table(liana_test_C4_WM)
  
}


lian<-rbindlist(lian)

lian$lesion<-sapply(lian$lesion, function(x) gsub('NAWN', 'NAWM', x))
lian_sub<-lian %>% group_by(lesion, source, target, ligand.complex, receptor.complex) %>% 
  summarise(logFC= mean(sca.LRscore), rank=min(aggregate_rank), specificity = mean(natmi.edge_specificity, ))%>% 
  group_by(source, target) %>%slice_min(rank, n=10) 

lian_sub<-as.data.table(lian_sub)
celltypes<-c('Astrocytes', 'Oligodendrocytes', 'Immune cells', 'OPCs')
for(ct in celltypes){
  gp<-ggplot(lian_sub[source==ct], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
    facet_grid(c('source','lesion'),scales = 'free')+scale_color_viridis_b()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))+
    ylab("")+xlab('')
  ggsave(gp, file=file.path(liana.dir,  paste0(ct, '.pdf')), width = 30, height = 30, units = 'cm' )
}


