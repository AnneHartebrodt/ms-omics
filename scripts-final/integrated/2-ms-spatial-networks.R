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
library(proxy)
library(ComplexHeatmap)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)
require(harmony)
library(UCell)
library(dplyr)
library(tidyverse)
library(cowplot)
library(patchwork)
require(ggplot2)
library(proxy)
library(igraph)

# enable parallel processing for network analysis (optional)
enableWGCNAThreads(nThreads = 12)

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna-spatial/networks')
dir.create(out.dir, showWarnings = F, recursive = T)

de.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/plots')
base.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/')


spatial<- c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')
#spatial<-c('RL_1D')
for(ds in spatial){
  
out.dir.ds<-file.path(out.dir, ds)
dir.create(out.dir.ds, recursive = T)


# READ AN PROCESS REQUIRED DATA
seurat_obj<-readRDS(file.path(base.dir, ds, paste0(ds, '_processed_.rds')))
CARD<-readRDS(file.path(base.dir, ds, paste0(ds, '_card_object.rds')))

# make a dataframe containing the image coordinates for each sample
image_df <- do.call(rbind, lapply(names(seurat_obj@images), function(x){
  seurat_obj@images[[x]]@coordinates
}))

# merge the image_df with the Seurat metadata
new_meta <- merge(seurat_obj@meta.data, image_df, by='row.names')

maj_celltype<-CARD@Proportion_CARD
row.name<-rownames(maj_celltype)
maj_celltype<-as.data.table(maj_celltype)
maj_celltype$Row.names<-row.name
maj_celltype[, major_celltype:=  names(.SD)[max.col(.SD)], .SDcols = 1:4]
new_meta<-merge(new_meta, maj_celltype, by = 'Row.names', all.x = T)
new_meta[which(is.na(new_meta$major_celltype)),]$major_celltype<-'No prediction'

new_meta<-as.data.frame(new_meta)

# fix the row ordering to match the original seurat object
rownames(new_meta) <- new_meta$Row.names
ix <- match(as.character(colnames(seurat_obj)), as.character(rownames(new_meta)))
new_meta <- new_meta[ix,]

# add the new metadata to the seurat object
seurat_obj@meta.data <- new_meta

# Set the IDENT to the predicted major celltype
Idents(seurat_obj) <- 'major_celltype'
seurat_obj@active.assay<-'Spatial'

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "vis"
)


seurat_obj <- MetaspotsByGroups(
  seurat_obj,
  group.by = c("major_celltype"),
  ident.group = "major_celltype",
  assay = 'Spatial',
  slot = 'counts'
)

seurat_obj  <- NormalizeMetacells(seurat_obj)

# set up the expression matrix, set group.by and group_name to NULL to include all spots
seurat_obj  <- SetDatExpr(
  seurat_obj,
  group.by=NULL,
  group_name = NULL,
  assay = 'Spatial',
  slot = 'counts')

# test different soft power thresholds
seurat_obj <- TestSoftPowers(seurat_obj)
plot_list <- PlotSoftPowers(seurat_obj)

powerplot<-wrap_plots(plot_list, ncol=2)
ggsave(powerplot, filename = file.path(out.dir.ds, paste0(ds, 'power_plot.pdf')), width = 30, height = 20, units = 'cm')

# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name=ds,
  overwrite_tom=TRUE
)


# plot the dendrogram
pdf(file.path(out.dir.ds, paste0(ds, 'dendrogram.pdf')))
dendrogram<-PlotDendrogram(seurat_obj, main='Spatial hdWGCNA dendrogram')
dev.off()

seurat_obj <- ModuleEigengenes(seurat_obj)
seurat_obj <- ModuleConnectivity(seurat_obj)

seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(ds, "-")
)

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')


# get module eigengenes and gene-module assignment tables
MEs <- GetMEs(seurat_obj)
modules <- GetModules(seurat_obj)
modules<-as.data.table(modules)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add the MEs to the seurat metadata so we can plot it with Seurat functions
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
dp <- DotPlot(seurat_obj, features=mods, dot.min=0.1,  idents = c('Oligodendrocytes', 'Astrocytes'))
# flip the x/y axes, rotate the axis labels, and change color scheme:
dp <- dp +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue') +
  xlab('') + ylab('')

ggsave(dp, filename = file.path(out.dir.ds, paste0(ds, 'module_expression_dotplot.pdf')), width = 20, height = 10, units = 'cm')

spfe<-SpatialFeaturePlot(
  seurat_obj,
  features = mods,
  alpha = c(0.1, 1),
)

ggsave(spfe, filename = file.path(out.dir.ds, paste0(ds, 'module_feature_plot.pdf')), width = 20, height = 10, units = 'cm')


celltype_color = c(as.character(wes_palette('GrandBudapest1', 4)), as.character(wes_palette('GrandBudapest2', 4)))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4],
                    'Neurons' = celltype_color[5], 'Endothelial cells' = celltype_color[6], 'Pericytes'=celltype_color[7], 'B-cells'=celltype_color[8])

#Make plot with majority celltype prediction
deconv.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna-spatial/deconvolution', ds)
dir.create(deconv.dir, showWarnings = F, recursive = T)
sp<-SpatialDimPlot(seurat_obj, group.by = 'major_celltype', cols = celltype_color)+
  guides(fill=guide_legend(title="Majority celltype"))
ggsave(sp, file=file.path(deconv.dir, paste0(ds, 'majority_celltype.pdf')), width = 20, height = 20, units = 'cm')


p <- PlotKMEs(seurat_obj, ncol=5)
ggsave(p, file=file.path(out.dir.ds, paste0(ds, 'module_eigengene_plot.pdf')), width = 20, height = 20, units = 'cm')


plot_list<- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE# order so the points with highest hMEs are on top
)
# stitch together with patchwork
mfp<-wrap_plots(plot_list)
ggsave(mfp, file=file.path(out.dir.ds, paste0(ds, 'module_feature_plot.pdf')), width = 20, height = 20, units = 'cm')

# vlp <- VlnPlot(
#   seurat_obj,
#   features = 'CA_2A1',
#   group.by = 'major_celltype',
#   pt.size = 0,
#   idents = c('Oligodendrocytes', 'Astrocytes') # don't show actual data points
# ) + geom_boxplot(width=.25, fill='white')+
#   xlab('') + ylab('hME') + NoLegend()
# 
# vlp
# ggsave(mfp, file=file.path(deconv.dir, paste0(ds, 'module_feature_plot.pdf')), width = 20, height = 20, units = 'cm')

# hubgene network
pdf(file=file.path(out.dir.ds, paste0(ds, 'module_hubgenes_network.pdf')))
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=5,
  edge_prop = 0.9,
  mods = "all"
)
dev.off()

ModuleNetworkPlot(seurat_obj, outdir = out.dir.ds)

}

