setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(harmony)
require(tidyverse)
require(data.table)
library(cowplot)
library(rsvd)
library(Rtsne)
library(Signac)
library(Seurat)
library(ggplot2)
library(patchwork)
library(GenomicRanges)
library(future)
require(biovizBase)
require(wesanderson)
require(gtable)

library(magrittr)
library(liana)
require(SeuratObject)
require(stringr)
require(ComplexHeatmap)
require(stringr)
require(RColorBrewer)
require(colorspace)
require(CARD)
require(MuSiC)
library(scran)
library(scater)
require(edgeR)

filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

subcluster.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/subclustering/')
dir.create(subcluster.dir)
adata<-readRDS(file=file.path(subcluster.dir, 'fine_annotation.rds'))


# CHECK cluster composition
batch_comp<-ggplot(adata@meta.data, aes(clu.name, fill=batch))+geom_bar()+
  scale_fill_viridis_d()+theme_bw()

ggsave(batch_comp, filename = file.path(subcluster.dir, 'composition_of_subclusters.pdf'))

Idents(adata)<-'clu.name'
markers <- FindAllMarkers(object = adata)
genenames<-rownames(markers)
markers<-as.data.table(markers)
markers$gene<-genenames
fwrite(markers, file=file.path(subcluster.dir, 'marker_genes.tsv'), sep='\t')


top_markers<- markers %>% group_by(cluster) %>%top_n(n=20, wt = avg_log2FC) 
top_markers<-as.data.table(top_markers)
fwrite(top_markers, file=file.path(subcluster.dir, 'top_20_markers_per_cellstate.tsv'), sep='\t')



reslist<-list()
keeplist<-list()
for (ct in unique(adata$clu.name)[!unique(adata$clu.name) %in% c('IM-5', 'IM-3', 'IM-2', 'AS-4')]){
  print(ct)
  
  cell.dir<-file.path(subcluster.dir, ct)
  dir.create(cell.dir)
  
  list_of_de<-list()

    
    astro<-subset(adata, subset = clu.name==ct)
    
    Idents(astro)<-'ms'
    astro$batch<-as.factor(astro$batch)
    sse <- as.SingleCellExperiment(astro)
    #cdr<-calc.cdr(astro@assays$RNA@counts)
    #sse$cellular_detection_rate<-as.numeric(cdr)
    
    summed <- aggregateAcrossCells(sse, 
                                   id=colData(sse)[,c("clu.name", "batch")])
    
    
    current <- summed

    y <- DGEList(counts(current), samples=colData(current))
    y
    
    discarded <- current$ncells < 10
    y <- y[,!discarded]
    summary(discarded)
    
    
    keep <- filterByExpr(y, group=current$batch, min.prop=0.5)
    y <- y[keep, ,keep.lib.sizes=FALSE]
    summary(keep)
    keeplist[[ct]]<-keep
    
    
    y <- calcNormFactors(y)
    y$samples
    
    par(mfrow=c(2,3))
    for (i in seq_len(ncol(y))) {
      plotMD(y, column=i)
    }
    
    
    design <- model.matrix(~1+as.factor(sex)+as.factor(ms), y$samples)
    design
    
    y <- estimateDisp(y, design)
    
    summary(y$trended.dispersion)
    
    #plotBCV(y)
    
    fit<-glmFit(y, design)
    #fit <- glmQLFit(y, design,robust = TRUE)
    summary(fit$var.prior)
    
    #summary(fit$df.prior)
    #plotQLDisp(fit)
    
    res <- glmLRT(fit, coef=ncol(design))
    #res <- glmQLFTest(fit, coef=ncol(design))
    reslist[[ct]]<-res
    s<-summary(decideTests(res))
    
    ta<-topTags(res, n=1000)
    ta<-ta@.Data[[1]]
    na<-rownames(ta)
    ta$gene<-na
    ta$celltype<-ct
    ta<-as.data.table(ta)
    ta<-ta[FDR<0.01]

  
  
  
  astro<-subset(adata, subset = clu.name==ct)
  
  Idents(astro)<-'ms'
  
  degenes<-ta
  degenes<-as.data.table(degenes)
  degenes$up<-ifelse(degenes$logFC>0, 'upregulated', 'downregulated')
  degenes[,.N, by=c('gene', 'up')][order(N, decreasing = T)]
  fwrite(degenes, file=file.path(cell.dir, paste0('DE_genes_', ct, '.tsv')), sep='\t')
  
}

rl<-list()
for (ct in unique(adata$clu.name)[!unique(adata$clu.name) %in% c('IM-5', 'IM-3', 'IM-2', 'AS-4')]){
  rl[[ct]]<-reslist[[ct]]$samples[,c("batch.1", "ncells")]
}

df<-rl[[1]]
for(r in rl[2:length(rl)]){
  df<-merge(df, r, by= 'batch.1', all=T)
}
colnames(df)<-c('batch', names(rl))
df[is.na(df)]<-0
fwrite(df, file.path(subcluster.dir, 'cells_in_clusters_summary.tsv'), sep='\t')


broad<-list('OL'='Oligodendrocytes', 'AS'='Astrocytes', 'IM'="Immune cells")
for (cat in c('OL', 'AS', 'IM')){
delist<-list()
for (ct in paste0(cat, '-', 1:7)){
  cell.dir<-file.path(subcluster.dir, ct)
  dir.create(cell.dir)
  if (file.exists(file.path(cell.dir, paste0('DE_genes_', ct, '.tsv')))){
  degenes<-fread(file=file.path(cell.dir, paste0('DE_genes_', ct, '.tsv')), sep='\t')
  }
  delist[[ct]]<-degenes
}

degenes<-rbindlist(delist)
degenes[, .N, by = c('gene', 'up')][N>1]

cell.group.dir<-file.path(subcluster.dir, paste0('all-', cat))
dir.create(cell.group.dir)
sub<-subset(adata, subset = celltypes_annotated == broad[cat][[1]])

exp <- FetchData(sub, unique(degenes$gene))
exp<-exp>0
exp<-as.data.table(exp)
exp$ms<-sub$ms

sum<-exp[, lapply(.SD, mean), by=ms]
sum<-as.data.table(sum)
sum<-sum[,-c('ms')]
sumsum<-sum[, lapply(.SD, sum)]
sumsum<-sumsum/nrow(sum)
deg<-degenes[gene %in% colnames(sumsum[,which(sumsum>0.05), with=F])]



Idents(sub)<-'clu.name'
sub$meta<-paste0(sub$clu.name, sub$ms)
sub$meta<-as.ordered(sub$meta)
Idents(sub)<-'meta'
if (length(unique(deg[up=='upregulated']$gene)>0)){
multi_up<-DotPlot(sub, features =unique(deg[up=='upregulated']$gene), scale = F)+theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))+
  scale_color_continuous_divergingx('RdBu')+ggtitle('Upregulated genes in MS')
multi_up
ggsave(multi_up, file=file.path(cell.group.dir, paste0(cat, '_upregulated_in_MS.pdf')) ,width =10+(max(0, (length(unique(deg[up=='upregulated']$gene))-5))), height = 11, units='cm')
}
if (length(unique(deg[up=='downregulated']$gene)>0)){
multi_down<-DotPlot(sub, features =unique(deg[up=='downregulated']$gene), scale = F)+theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust = 1))+
  scale_color_continuous_divergingx('RdBu')+ggtitle('Downregulated genes in MS')
multi_down
ggsave(multi_down, file=file.path(cell.group.dir, paste0(cat, '_downregulated_in_MS.pdf')) ,width =10+(max(0, (length(unique(deg[up=='upregulated']$gene))-5))), height = 11, units='cm')
}
}





#### Look at the ones which did not work

reslist<-list()
keeplist<-list()
for (ct in c('IM-5', 'IM-3', 'IM-2', 'AS-4')){
  print(ct)
  
  cell.dir<-file.path(subcluster.dir, ct)
  dir.create(cell.dir)
  
  list_of_de<-list()
  
  
  astro<-subset(adata, subset = clu.name==ct)
  
  Idents(astro)<-'ms'
  astro$batch<-as.factor(astro$batch)
  sse <- as.SingleCellExperiment(astro)
  #cdr<-calc.cdr(astro@assays$RNA@counts)
  #sse$cellular_detection_rate<-as.numeric(cdr)
  
  summed <- aggregateAcrossCells(sse, 
                                 id=colData(sse)[,c("clu.name", "batch")])
  
  
  current <- summed
  
  y <- DGEList(counts(current), samples=colData(current))
  y
  
  discarded <- current$ncells < 10
  y <- y[,!discarded]
  summary(discarded)
  
  
  keep <- filterByExpr(y, group=current$batch, min.prop=0.5)
  y <- y[keep, ,keep.lib.sizes=FALSE]
  summary(keep)
  keeplist[[ct]]<-keep
  
  
  y <- calcNormFactors(y)
  y$samples
  
  reslist[[ct]]<-y
}

rl<-list()
for (ct in  c('IM-5', 'IM-3', 'IM-2', 'AS-4')){
  rl[[ct]]<-reslist[[ct]]$samples[,c("batch.1", "ncells")]
}

df<-rl[[1]]
for(r in rl[2:length(rl)]){
  df<-merge(df, r, by= 'batch.1', all=T)
}
colnames(df)<-c('batch', c('IM-5', 'IM-3', 'IM-2', 'AS-4'))
df[is.na(df)]<-0
fwrite(df, file.path(subcluster.dir, 'clusters_only_MS.tsv'), sep='\t')

