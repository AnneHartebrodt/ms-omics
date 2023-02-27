setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(data.table)
require(wesanderson)
library(ggplot2)         
require(colorspace)
require(stringr)
require(ggplot2)
require(Seurat)
require(stringr)

filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

state.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/specific_expression/')
dir.create(state.dir)

base.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/')


adata<- readRDS(file.path(filtering.dir, "annotated.rds"))
adata$celltype_short<-recode(adata$celltypes_annotated, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPC', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')
#saveRDS(adata, file = file.path(filtering.dir, "annotated.rds"), compress = F)

genes<-c('C1QA', 'C1QB', 'C1QC', 'APOE', 'APOD', 'C3', 'MTRNR2L1', 'SLC7A11', 'CRYAB', 'TF', 'TFRC', 'LRP2', 'NRCAM', 'NRXN3', 'ERBB4', 'LINC01608', 'CNTN2',
         'SPP1', "TF", 'SPP1', 'CD44', 'CHI3L1', "APOE", "C3", "C1QB", "C1QA", "C1QC", 'PLP1', 'CNP',
         "FOXP2", "LUZP2", "BCL6", "CEBPD", "FLI1", "ETV6", "JUN", "JUNB", "TCF7L1", "SOX6", "FOS", "ETS2", "IKZF1", "RUNX1", "EPAS1")
genes<-unique(genes)
for(gene in genes){
APOD<-VlnPlot(object = adata, features = gene, group.by = 'celltype_short', split.by = 'ms')+
  ylim(0.01, 4)+xlab('')+ylab("Expression Level\n(0 truncated)")+
  theme(axis.text.x = element_text(angle=0, hjust = 0.5, size = 12), legend.margin = margin(0,0,0,0))
APOD
ggsave(APOD, file=file.path(state.dir, paste0(gene, '.pdf')), width=5, height = 2.25, units = 'in')
}



spatial.dir<-file.path(state.dir, 'spatial_feature_plots')
dir.create(spatial.dir)
spatial<- c('CA_2A','CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')

features <- c('MTRNR2L1', 'SLC7A11', 'NFL', 'NFM', 'NPY')
features <- c( 'CXCL14')

features <- c( 'NEFL', 'NEFM')

features<-c('PLP1', 'CNP', 'MT-ND1', 'MT-CO1', 'SNAP25', 'NRGN', "SPP1", 'CNTN2', 'TF')


features<-c('TLE4', 'VIM')
for(ds in spatial){
  # READ AN PROCESS REQUIRED DATA
  seurat_obj<-readRDS(file.path(base.dir, ds, paste0(ds, '_processed_.rds')))
  for (f in features){
  dir.create(file.path(spatial.dir, f))
  sp<-SpatialFeaturePlot(seurat_obj, features = f)
  ggsave(sp, file=file.path(spatial.dir, f,  paste0(f,'_', ds, '.pdf')), width = 10, height = 10, units = 'cm')
  }
}





astrocytes.genes<-c("TEF", "PRR12", "MAN2A1", "LRRC69", "SCFD2")
immune.genes<-c("LPIN1")
oligo.genes<-c("TNR", "CRB2", "CDH23", "INPP5A", "RNH1", "OR4N2", "KIF19", "PIAS4", "NETO1", "GTPBP1", "KCNIP3", "KCNS3", "SLCO5A1", "KIF5C", "SLC7A11", "FAM167A", "GPAT4")


astro<-DotPlot(subset(adata, subset = celltypes_annotated=="Astrocytes"), 
        features = astrocytes.genes,
        group.by = 'lesion_type')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_fill_continuous_divergingx(palette = 'RdBu') +
  scale_color_continuous_divergingx(palette = 'RdBu')+xlab("")+ylab("")
ggsave(astro, file=file.path(state.dir, 'astrocytes.pdf'), width=15, height = 8, units = 'cm')


id<-subset(adata, subset = celltypes_annotated=="Immune cells")
immune<-DotPlot(id, 
        features = immune.genes,
        group.by = 'lesion_type')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_fill_continuous_divergingx(palette = 'RdBu') +
  scale_color_continuous_divergingx(palette = 'RdBu')+xlab("")+ylab("")

ggsave(immune, file=file.path(state.dir, 'immunecells.pdf'), width=15, height = 8, units = 'cm')


oligo<-DotPlot(subset(adata, subset = celltypes_annotated=="Oligodendrocytes"),
        features = oligo.genes,
        group.by = 'meta')+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_fill_continuous_divergingx(palette = 'RdBu') +
  scale_color_continuous_divergingx(palette = 'RdBu')+xlab("")+ylab("")
ggsave(oligo, file=file.path(state.dir, 'oligodendrocytes.pdf'), width=15, height = 8, units = 'cm')

