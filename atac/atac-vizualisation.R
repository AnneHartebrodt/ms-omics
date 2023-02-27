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


out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
vizualisation.dir<-file.path(out.dir, 'track_vizualisations_candidates')
dir.create(vizualisation.dir)

hm.integrated<-readRDS( file = file.path(out.dir, "integrated.rds"))

hm.integrated$meta<-as.ordered(paste0(hm.integrated$predicted_celltype, '.', hm.integrated$is_ms))

hm.integrated$celltype_short<-recode('Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPCs', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')
celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('Oligodendrocytes.MS'=celltype_color[1], 'Astrocytes.MS'=celltype_color[2], 'Immune cells.MS'=celltype_color[3], 'OPCs.MS'=celltype_color[4],
                    'Oligodendrocytes.C'=celltype_color[1], 'Astrocytes.C'=celltype_color[2], 'Immune cells.C'=celltype_color[3], 'OPCs.C'=celltype_color[4])


hm.integrated@active.assay<-'peaks'
Idents(hm.integrated)<-'meta'
annotations <- Annotation(object = hm.integrated)

#genes<-c("SOX9", "SOX10", "HIF1A", "STAT1","STAT3", "REL", "RELA", "EGR1", "PLAG1", "CTCFL", "ZIC5", "EGR3","TF", "TFRC", "LRP2", "SPP1", "APOE", "C3", "C1QB", "CD44", "SPP1")
genes<-c("KLF2", "KLF4", "KLF7", "KLF9", "KLF10", "KLF15", "KLF16", "KLF17", "SP1", "SP2" ,"SP5", "SP8", "CEBPD" )
for (ge in genes){
  
isgene <- annotations$gene_name == ge
isgene <- !is.na(x = isgene) & isgene
annot.sub <- annotations[isgene]
gr <- GRanges(seqnames = as.character(x = seqnames(x = annot.sub))[[1]],
                ranges = IRanges(start = min(start(x = annot.sub)),
                                 end = max(end(x = annot.sub))))
  
cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = gr,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
)+scale_fill_manual(values = celltype_color)


expr_plot <- ExpressionPlot(
  object = hm.integrated,
  features = ge,
  assay = "RNA",
  idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
)+scale_fill_manual(values = celltype_color)


gene_plot <- AnnotationPlot(
  object = hm.integrated,
  region = gr,
)


p<-CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  expression.plot = expr_plot,
  widths = c(10, 2)
)

ggsave(p, file=file.path(vizualisation.dir, paste0(ge, '.pdf')), width = 20, height = 20, units = 'cm')
}
