setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)


require(rlang)
require(Seurat)
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

result.dir<- '/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/'
plot.dir<- file.path(result.dir, 'plots')
dir.create(plot.dir, recursive = T)

sobjects<-list()

# NAWM did not really work 'NAWM_2D'

make_venn<-function(top.de, foldchange=1, outf){
  
  venni <-list()
  for (c in unique(top.de$cluster)){
    venni[[c]]<- top.de[cluster==c]$gene
  }
  
  myV <- plotVenn(venni, setColors =DiscretePalette(length(unique(sdata$seurat_clusters))), 
                  outFile = outf)
  reg<-listVennRegions(myV)
  
  dflist<-list()
  for (i in 1:length(reg)){
    df<-data.table(strsplit( names(reg[i]), '\\(|\\)')[[1]][2], reg[[i]])
    dflist[[i]]<-df
  }
  df<-rbindlist(dflist)
  
  
  gv<-ggvenn(venni,  fill_color = DiscretePalette(length(unique(sdata$seurat_clusters))), 
             show_elements = TRUE, 
             label_sep = "\n",
             text_size = 2
  )
  gv2<-ggvenn(venni,  fill_color =DiscretePalette(length(unique(sdata$seurat_clusters))),
              label_sep = "\n",
              text_size = 3)
  return(list(gv2, gv, df))
}

for(ds in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  #Make a plot directory
  plot.dir.ds<- file.path(plot.dir, ds)
  dir.create(plot.dir.ds, showWarnings = F)
  
  img<- Read10X_Image(image.dir = file.path( '/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-spatial/',ds,'/outs/spatial/'),
                      image.name = 'tissue_lowres_image.png',
                      filter.matrix = TRUE)
  sdata<-Load10X_Spatial(
    data.dir = file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-spatial/', ds, 'outs'),
    filename = "filtered_feature_bc_matrix.h5",
    assay = "Spatial",
    slice = "slice1",
    filter.matrix = TRUE,
    to.upper = FALSE,
    image = img
  )
  
  dir.create(file.path(result.dir, ds), showWarnings = F, recursive = T)
  saveRDS(sdata, file = file.path(result.dir, ds, paste0(ds, '.rds')), compress = FALSE)
  
  sdata[["percent.mt"]] <- PercentageFeatureSet(sdata, pattern = "^MT-")
  p1<-VlnPlot(sdata, features = c( "nCount_Spatial", "percent.mt"), ncol = 3)
  ggsave(p1, file=file.path(plot.dir.ds, paste0(ds, '_violin_diagnostic.pdf')))
  p2<-SpatialFeaturePlot(sdata, features = "nCount_Spatial") + theme(legend.position = "right")
  ggsave(p2, file=file.path(plot.dir.ds, paste0(ds, '_umap_nCountSpatial.pdf')))

  # normalize and embedd
  if (ds %in% c('RL_1B', 'RL_1D')){
    # For some reason sc tranform fails, probably due to 0 UMIs (fails on log)
    sdata<- FindVariableFeatures(sdata)
    sdata<- NormalizeData(sdata)
    sdata<-ScaleData(sdata)
  }else{
  sdata <- SCTransform(sdata, assay = "Spatial", verbose = TRUE)
  }
  sdata <- RunPCA(sdata,  verbose = FALSE)
  sdata <- FindNeighbors(sdata, reduction = "pca", dims = 1:30)
  sdata <- FindClusters(sdata, verbose = FALSE, modularity.fxn = 1)
  sdata <- RunUMAP(sdata, reduction = "pca", dims = 1:30)
  cols<-DiscretePalette(length(unique(sdata$seurat_clusters)))
  names(cols)<-levels(Idents(sdata))
  p3 <- DimPlot(sdata, reduction = "umap", label = TRUE,cols = cols)

  p4 <- SpatialDimPlot(sdata, label = TRUE, label.size = 3, cols=cols)
  p34<-p3+p4

  ggsave(p34, file=file.path(plot.dir.ds, paste0(ds, '_umap_clusters.pdf')))

  de_markers <- FindAllMarkers(sdata)
  de_markers<- as.data.table(de_markers)
  fwrite(de_markers, file=file.path(plot.dir.ds, paste0(ds, '_differential_genes.tsv')), sep='\t')


  top.de <-de_markers %>% group_by(cluster) %>% top_n(10, wt = abs(avg_log2FC))
  top.de<-as.data.table(top.de)
  p5<-SpatialFeaturePlot(sdata, features =unique(top.de$gene), slot = 'scale.data')
  p5

  ggsave(p5, file=file.path(plot.dir.ds, paste0(ds, '_umap_top_de_genes.pdf')), width=25, height = 25)


  v<-make_venn(top.de, foldchange = 1, file.path(plot.dir.ds, paste0(ds, '_venn_top_de_genes.svg')))
  ggsave(v[[2]], file=file.path(plot.dir.ds, paste0(ds, '_venn_top_de_genes.pdf')))
  fwrite(v[[3]], file=file.path(plot.dir.ds, paste0(ds, '_overlap_top_de_cluster.tsv')), sep='\t')
  sobjects[[ds]]<-sdata
  saveRDS(sdata, file = file.path(result.dir, ds, paste0(ds, '_processed_.rds')), compress = FALSE)
  
}

saveRDS(sobjects, file=file.path(result.dir, 'spatial_objects.Rds'))



