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
require(clusterProfiler)

de.dir<-'/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/plots'
result.dir<- '/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna-spatial/'
result.dir<- file.path(result.dir, 'expression-of-spatial-in-RNA')
dir.create(result.dir, recursive = T)

filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/', 'filtering')


adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

for(ds in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  #Make a plot directory
  result.dir.ds<- file.path(result.dir, ds)
  dir.create(result.dir.ds, showWarnings = F)
  de.genes<-fread(file.path(de.dir, ds, paste0(ds, '_differential_genes.tsv')), sep='\t', header=TRUE)
  de.genes<-de.genes[abs(avg_log2FC)>0.5]
  
  de<-unique(de.genes[,.(gene, cluster)]) %>% pivot_wider(names_from = c('gene'), values_from = c('cluster'), values_fn = ~paste0(.x, collapse = ','))
  de<-t(de)
  rn<-de[,1]
  de<-rownames(de)
  names(de)<-rn
  
  reg.outdir<-file.path(result.dir, ds)
  dir.create(reg.outdir)
  
  
  if(length(de)>0){
  
  if(length(de)>400){
    
    for (clu in unique(names(de))){
      dp1<-DotPlot(adata ,features = de[names(de)==clu], 
                   scale = F, group.by = 'celltypes_annotated')+theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))+
        scale_color_continuous_divergingx('RdBu')+
        theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5))
      
      ggsave(dp1, file=file.path(reg.outdir, paste0(clu, 'cluster_expression.pdf')), width = 50, height = 10, units = 'cm')
    }    
  }
     
  
  else{  
  dp1<-DotPlot(adata ,features = de[sapply(names(de), function(x) nchar(x)==1)], 
               scale = F, group.by = 'celltypes_annotated')+theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))+
    scale_color_continuous_divergingx('RdBu')+
    theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5))
  
  ggsave(dp1, file=file.path(reg.outdir, paste0('cluster_specific_expression.pdf')), width = 50, height = 10, units = 'cm')
  
  
  dp2<-DotPlot(adata ,features = de[sapply(names(de), function(x) nchar(x)>1)], 
             scale = F, group.by = 'celltypes_annotated', 
  )+theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1),strip.text.x = element_text(angle=90, hjust = 0, vjust = 0.5))+
   scale_color_continuous_divergingx('RdBu')
  
  ggsave(dp2, file=file.path(reg.outdir, paste0('cluster_multiple_expression.pdf')), width = 50, height = 10, units = 'cm')
  }
  }
}

