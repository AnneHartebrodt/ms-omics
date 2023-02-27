setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)
require(Seurat)
library(tidyverse)
require(data.table)
require(ggrepel)
require(colorspace)
library(circlize)
library(wesanderson)
require(RColorBrewer)
require(harmony)
library(UCell)
library(dplyr)
library(cowplot)
library(patchwork)
require(ggplot2)
# load
library(MuSiC)
library(CARD)



## COPY FROM GITHUB TO allow for finetuning of the plot i.e change of point size
CARD.visualize.prop <- function(proportion,spatial_location,ct.visualize = ct.visualize,colors = c("lightblue","lightyellow","red"),NumCols){
  if(is.null(colors)){
    colors = c("lightblue","lightyellow","red")
  }else{
    colors = colors
  }
  res_CARD = as.data.frame(proportion)
  res_CARD = res_CARD[,order(colnames(res_CARD))]
  location = as.data.frame(spatial_location)
  if(sum(rownames(res_CARD)==rownames(location))!= nrow(res_CARD)){
    stop("The rownames of proportion data does not match with the rownames of spatial location data")
  }
  ct.select = ct.visualize
  res_CARD = res_CARD[,ct.select]
  res_CARD_scale = as.data.frame(apply(res_CARD,2,function(x){
    (x - min(x)) / (max(x) - min(x))
  } ))
  res_CARD_scale$x = as.numeric(location$x)
  res_CARD_scale$y = as.numeric(location$y)
  mData = melt(res_CARD_scale,id.vars = c("x","y"))
  colnames(mData)[3] <- "Cell_Type"
  b = c(0,1)
  p = suppressMessages(ggplot(mData, aes(x, y)) + 
                         geom_point(aes(colour = value),size = 1.0) +
                         scale_color_gradientn(colours = colors) + 
                         #scale_color_viridis_c(option = 2)+
                         scale_x_discrete(expand = c(0, 1)) + scale_y_discrete(expand = c(0,1))+ 
                         facet_wrap(~Cell_Type,ncol = NumCols)+ 
                         coord_fixed()+
                         theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
                               #legend.position=c(0.14,0.76),
                               panel.background = element_blank(),
                               plot.background = element_blank(),
                               panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
                               axis.text =element_blank(),
                               axis.ticks =element_blank(),
                               axis.title =element_blank(),
                               legend.title=element_text(size = 14,face="bold"),
                               legend.text=element_text(size = 11),
                               strip.text = element_text(size = 12,face="bold"),
                               legend.key = element_rect(colour = "transparent", fill = "white"),
                               legend.key.size = unit(0.45, 'cm')))
  return(p)
}

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna-spatial/deconvolution')
dir.create(out.dir, showWarnings = F, recursive = T)

de.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/plots')
base.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/')
filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

# READ annotated adata object
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

sc_count<-adata@assays$RNA@counts
sc_meta<-adata@meta.data

for(ds in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  out.dir.ds<-file.path(out.dir, ds)
  dir.create(out.dir.ds)
seurat_obj<-readRDS(file.path(base.dir, ds, paste0(ds, '_processed_.rds')))

spatial_count<-seurat_obj@assays$Spatial@counts

# make a dataframe containing the image coordinates for each sample
image_df <- do.call(rbind, lapply(names(seurat_obj@images), function(x){
  seurat_obj@images[[x]]@coordinates
}))

spatial_location<-image_df[, c(3,2)]
colnames(spatial_location)<-c('x', 'y')

#mirror the coordinates
spatial_location$y<-sapply(spatial_location$y, function(a) abs(max(spatial_location$y)-a))

CARD_obj = createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "celltypes_annotated",
  ct.select = c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs', 'Neurons', 'Endothelial cells', 'Pericytes', 'B-cells'),
  sample.varname = "batch",
  minCountGene = 100,
  minCountSpot = 1) 


CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)

celltype_color = c(as.character(wes_palette('GrandBudapest1', 4)), as.character(wes_palette('GrandBudapest2', 4)))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4],
                    'Neurons' = celltype_color[5], 'Endothelial cells' = celltype_color[6], 'Pericytes'=celltype_color[7], 'B-cells'=celltype_color[8])




ct.visualize = c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs', 'Neurons', 'Endothelial cells', 'Pericytes', 'B-cells')

## visualize the spatial distribution of the cell type proportion
p2 <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,        
  spatial_location = CARD_obj@spatial_location, 
  ct.visualize = ct.visualize,                 ### selected cell types to visualize
  colors = c("lightblue","lightyellow","red"), ### if not provide, we will use the default colors
  NumCols = 3)      ### number of columns in the figure panel

ggsave(p2, file = file.path(out.dir.ds, 'celltype_location_prediction.pdf'), width=30, height=20, units='cm')
saveRDS(CARD_obj, file= file.path(base.dir, ds, paste0(ds, '_card_object.rds')))
}


