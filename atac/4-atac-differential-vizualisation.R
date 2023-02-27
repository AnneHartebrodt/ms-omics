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
library(JASPAR2022)

library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
require(wesanderson)
require(GenomicRanges)
require(dplyr)

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
vizualisation.dir<-file.path(out.dir, 'track_vizualisations_from_DE')
dir.create(vizualisation.dir)

da.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/differential_celltypes')

hm.integrated<-readRDS( file = file.path(out.dir, "integrated.rds"))

hm.integrated$meta<-as.ordered(paste0(hm.integrated$predicted_celltype, '.', hm.integrated$is_ms))
hm.integrated$celltype_short<-recode(hm.integrated$predicted_celltype, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPCs', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')


celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('Oligodendrocytes.MS'=celltype_color[1], 'Astrocytes.MS'=celltype_color[2], 'Immune cells.MS'=celltype_color[3], 'OPCs.MS'=celltype_color[4],
                    'Oligodendrocytes.C'=celltype_color[1], 'Astrocytes.C'=celltype_color[2], 'Immune cells.C'=celltype_color[3], 'OPCs.C'=celltype_color[4])


hm.integrated@active.assay<-'peaks'
Idents(hm.integrated)<-'meta'
annotations <- Annotation(object = hm.integrated)

for (ct in c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs')){
  viz.dir.ct<-file.path(vizualisation.dir, ct)
  dir.create(viz.dir.ct)
  
  regions<-fread(file.path(da.dir, paste0(ct, '_DE_accessible_region_all.tsv')))
for (chrom in regions$region){
    cov_plot<-CoveragePlot(
    object = hm.integrated,
    region = chrom,
    annotation = FALSE,
    peaks = FALSE,
    idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
  )+scale_fill_manual(values = celltype_color)
  
  
  
  gene_plot <- AnnotationPlot(
    object = hm.integrated,
    region = chrom,
  )
  
  
  p<-CombineTracks(
    plotlist = list(cov_plot, gene_plot),
  )
  
  ggsave(p, file=file.path(viz.dir.ct, paste0(chrom, '.pdf')), width = 20, height = 20, units = 'cm')
}
}
