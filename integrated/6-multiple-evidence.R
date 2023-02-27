setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5)
require(rhdf5filters)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell-new')
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
#library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
#require(biovizBase)
#library(JASPAR2022)

#library(TFBSTools)
#library(motifmatchr)
#library(BSgenome.Hsapiens.UCSC.hg38)
require(wesanderson)
require(GenomicRanges)
require(dplyr)
require(cowplot)

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
vizualisation.dir<-file.path(out.dir, 'track_vizualisations_from_DE')
dir.create(vizualisation.dir)

da.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/differential_celltypes')

hm.integrated<-readRDS( file = file.path(out.dir, "integrated.rds"))

hm.integrated$celltype_short<-recode(hm.integrated$predicted_celltype, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPCs', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')
hm.integrated$meta<-as.ordered(paste0(hm.integrated$celltype_short, '.', hm.integrated$is_ms))


celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('OL.MS'=celltype_color[1], 'AS.MS'=celltype_color[2], 'IM.MS'=celltype_color[3], 'OPCs.MS'=celltype_color[4],
                    'OL.C'=celltype_color[1], 'AS.C'=celltype_color[2], 'IM.C'=celltype_color[3], 'OPCs.C'=celltype_color[4])


hm.integrated@active.assay<-'peaks'
Idents(hm.integrated)<-'meta'
annotations <- Annotation(object = hm.integrated)


gene<-'SLC7A11'
ct <- 'Oligodendrocytes'
viz.dir.ct<-file.path(vizualisation.dir, ct, gene)
dir.create(viz.dir.ct)

regions<-fread(file.path(da.dir, 'Oligodendrocytes_DE_accessible_region_annot.tsv'))
regions<-regions[gene_name==gene]
chrom_reg<-'chr4-138185885-138242708'
ranges.show <- StringToGRanges(regions$region[1])
ranges.show.2 <- StringToGRanges(regions$region[2])
ranges.show.3 <- StringToGRanges(regions$region[3])
ranges.show.4 <- StringToGRanges(regions$region[4])
ranges.show$color <-'mediumblue'
ranges.show.2$color <-'mediumblue'
ranges.show.3$color <-'mediumblue'
ranges.show.4$color <-'mediumblue'

cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = chrom_reg,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('OL.MS','OL.C'),
  region.highlight = c(ranges.show, ranges.show.2, ranges.show.3, ranges.show.4))+
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())+scale_fill_manual(values = celltype_color)

cov_plot

#ggsave(cov_plot, file=file.path(viz.dir.ct, paste0( 'SCL7A11_overview.pdf')), width = 30, height = 5, units = 'cm')


gene_plot <- AnnotationPlot(
  object = hm.integrated,
  region = gene,
)
#ggsave(gene_plot, file=file.path(viz.dir.ct, paste0( 'SCL7A11_overview.pdf')), width = 30, height = 5, units = 'cm')

combo<-CombineTracks(list(cov_plot, gene_plot), heights = c(0.5, 0.5))
#ggsave(combo, file=file.path(viz.dir.ct, paste0( 'SCL7A11_overview.pdf')), width = 30, height = 5, units = 'cm')

#combo<-CombineTracks(list(cov_plot, gene_plot), heights = c(1, 0.5))
#ggsave(combo, file=file.path(viz.dir.ct, paste0( 'SCL7A11_overview_2.pdf')), width = 30, height = 5, units = 'cm')
filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/', 'filtering')
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

adata$celltype_short<-recode(adata$celltypes_annotated, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPCs', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')
adata$meta<-as.ordered(paste0(adata$celltype_short, '.', adata$ms))

Idents(adata)<-'meta'

expr_plot1 <- ExpressionPlot(
  object = adata,
  features = gene,
  
  assay = "RNA",
  idents = c('OL.MS','OL.C')
)+scale_fill_manual(values = celltype_color)+
  scale_x_continuous(limits=c(0,2))+
  ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+ggtitle('RNA\nexpression')

expr_plot1

expr_plot <- ExpressionPlot(
  object = hm.integrated,
  features = gene,
  assay = "RNA",
  idents = c('OL.MS','OL.C')
)+scale_fill_manual(values = celltype_color)+ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+ggtitle('predicted\nexpression')



pl<-plot_grid(expr_plot1, expr_plot)
pl<-plot_grid(pl, NULL, nrow = 2, rel_heights =c(0.9, 0.1) )+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p<-CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  heights = c(0.65, 0.35)
)
pp<-plot_grid(p, pl, rel_widths = c(0.77, 0.23))
ggsave(pp, file=file.path(viz.dir.ct, paste0( 'SCL7A11_overview_expr.pdf')), width = 30, height = 5, units = 'cm')


sub<-subset(adata, subset=celltype_short=='OL')
VlnPlot(sub, features = gene, group.by = 'batch')


gene = 'MTRNR2L1'
regions<-fread(file.path(da.dir, 'Oligodendrocytes_DE_accessible_region_annot.tsv'))
regions<-regions[gene_name==gene]
cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = regions$region,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('OL.MS','OL.C'))+
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

gene_plot <- AnnotationPlot(
  object = hm.integrated,
  region = gene,
)+ theme(axis.text.y  = element_blank(),
         axis.title.y = element_blank(),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank())

viz.dir.ct<-file.path(vizualisation.dir, ct, gene)
dir.create(viz.dir.ct)

expr_plot1 <- ExpressionPlot(
  object = adata,
  features = gene,
  assay = "RNA",
  idents = c('OL.MS','OL.C')
)+scale_fill_manual(values = celltype_color)+
  scale_x_continuous(limits=c(0,2))+
  ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+
  ggtitle('RNA\nexpression')

expr_plot1

expr_plot <- ExpressionPlot(
  object = hm.integrated,
  features = gene,
  assay = "RNA",
  idents = c('OL.MS','OL.C')
)+scale_fill_manual(values = celltype_color)+ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+
  ggtitle('predicted\nexpression')



pl<-plot_grid(expr_plot1, expr_plot)
pl<-plot_grid(pl, NULL, nrow = 2, rel_heights =c(0.9, 0.1) )+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p<-CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  heights = c(0.65, 0.35)
)
pp<-plot_grid(p, pl, rel_widths = c(0.6, 0.4))

ggsave(pp, file=file.path(viz.dir.ct, paste0( 'MTRNR2L1_overview_expr.pdf')), width = 15, height = 5, units = 'cm')



adata$celltype_short<-recode(adata$celltypes_annotated, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPCs', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')
adata$meta<-as.ordered(paste0(adata$celltype_short, '.', adata$ms))

Idents(adata)<-'meta'

gene = 'VIM'
regions<-fread(file.path(da.dir, 'Oligodendrocytes_DE_accessible_region_annot.tsv'))
regions<-regions[gene_name==gene]
cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = regions$region,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('OL.MS','OL.C'))+
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())

gene_plot <- AnnotationPlot(
  object = hm.integrated,
  region = gene,
)+ theme(axis.text.y  = element_blank(),
         axis.title.y = element_blank(),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank())

viz.dir.ct<-file.path(vizualisation.dir, ct, gene)
dir.create(viz.dir.ct)

expr_plot1 <- ExpressionPlot(
  object = adata,
  features = gene,
  assay = "RNA",
  idents = c('OL.MS','OL.C')
)+scale_fill_manual(values = celltype_color)+
  scale_x_continuous(limits=c(0,2))+
  ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+ggtitle('RNA\nexpression')

expr_plot1

expr_plot <- ExpressionPlot(
  object = hm.integrated,
  features = gene,
  assay = "RNA",
  idents = c('OL.MS','OL.C')
)+scale_fill_manual(values = celltype_color)+ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+
  ggtitle('predicted\nexpression')



pl<-plot_grid(expr_plot1, expr_plot)
pl<-plot_grid(pl, NULL, nrow = 2, rel_heights =c(0.9, 0.1) )+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p<-CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  heights = c(0.65, 0.35)
)
pp<-plot_grid(p, pl, rel_widths = c(0.6, 0.4))

ggsave(pp, file=file.path(viz.dir.ct, paste0( 'VIM_overview_expr.pdf')), width = 15, height = 5, units = 'cm')




gene = 'TLE4'
regions<-fread(file.path(da.dir, 'Oligodendrocytes_DE_accessible_region_annot.tsv'))
regions<-regions[gene_name==gene]
cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = regions$region,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('IM.MS','IM.C'))+
  theme(axis.text.y  = element_blank(),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank())+scale_fill_manual(values = celltype_color)

gene_plot <- AnnotationPlot(
  object = hm.integrated,
  region = gene,
)+ theme(axis.text.y  = element_blank(),
         axis.title.y = element_blank(),
         axis.line.y = element_blank(),
         axis.ticks.y = element_blank())

viz.dir.ct<-file.path(vizualisation.dir, ct, gene)
dir.create(viz.dir.ct)

expr_plot1 <- ExpressionPlot(
  object = adata,
  features = gene,
  assay = "RNA",
  idents = c('IM.MS','IM.C')
)+scale_fill_manual(values = celltype_color)+
  scale_x_continuous(limits=c(0,2))+
  ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+
  ggtitle('RNA\nexpression')

expr_plot1

expr_plot <- ExpressionPlot(
  object = hm.integrated,
  features = gene,
  assay = "RNA",
  idents = c('IM.MS','IM.C')
)+scale_fill_manual(values = celltype_color)+ylab("")+xlab("")+
  theme(axis.ticks.x.top = element_blank(), axis.text.x.top = element_blank())+
  ggtitle('predicted\nexpression')



pl<-plot_grid(expr_plot1, expr_plot)
pl<-plot_grid(pl, NULL, nrow = 2, rel_heights =c(0.9, 0.1) )+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

p<-CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  heights = c(0.65, 0.35)
)
pp<-plot_grid(p, pl, rel_widths = c(0.57, 0.43))

ggsave(pp, file=file.path(viz.dir.ct, paste0( 'TLE4_Immune_overview_expr.pdf')), width = 17, height = 5, units = 'cm')






