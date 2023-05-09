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


filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

liana.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/liana')
dir.create(liana.dir)

adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

adata@active.ident<-adata@meta.data$celltypes_annotated
adata@active.ident<- adata$batch
adata@active.assay<-'RNA'
Idents(adata)<-'batch'
liana_objects<-list()

# 1. subset the combined object to only include one batch at the time
# Cell-cell interactions can only manifest physically if they are in the same sample so using the full object
# would only yield false positives.
for (b in unique(adata$batch)){
  subbi <- subset(adata, subset = batch==b)
  subbi@active.ident<-subbi@meta.data$celltypes_annotated
  liana_test_subbi <- liana_wrap(subbi)
  liana_objects[[b]]<-liana_test_subbi
}


# 2. Aggregate liana results, using the predefined funtion
for (b in unique(adata$batch)){
  lo<-liana_objects[[b]] %>%  liana_aggregate()
  fwrite(lo, file=file.path(liana.dir, paste0(b, '.tsv')), sep='\t')
}

saveRDS(liana_objects, file = file.path(liana.dir, 'all.Rds'))

#liana_objects<-readRDS(file.path(liana.dir, 'all.Rds'))

lian<-list()



for(b in unique(adata$batch)){
  liana_test_C4_WM <- liana_objects[[b]] %>%
    liana_aggregate()
  liana_test_C4_WM$batch<-b
  liana_test_C4_WM$lesion<-stringr::str_split(b, '_')[[1]][2]
  lian[[b]]<-as.data.table(liana_test_C4_WM)

}

# Make all kinds of plots
require(dplyr)
lian<-rbindlist(lian)
lian_sub<-lian %>% group_by(lesion, source, target, ligand.complex, receptor.complex) %>% 
  summarise(logFC= mean(sca.LRscore), rank=min(aggregate_rank), specificity = mean(natmi.edge_specificity, ))%>% 
  group_by(source, target) %>%slice_min(rank, n=10) 

lian_sub<-as.data.table(lian_sub)

celltypes<-c('AS', 'OL', 'IM', 'OPCs')
lian_sub$source<-recode(lian_sub$source, Astrocytes='AS', Neurons='N', Oligodendrocytes='OL', "Immune cells"='IM', "B-cells"='BC', Pericytes="PC", "Endothelial cells"='EC')
lian_sub$target<-recode(lian_sub$target, Astrocytes='AS', Neurons = 'N', Oligodendrocytes='OL', "Immune cells"='IM', "B-cells"='BC', Pericytes="PC", "Endothelial cells"='EC')

for(ct in unique(lian_sub$source)){
gp<-ggplot(lian_sub[source==ct], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('source', 'lesion'),scales = 'free')+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  ylab("")+xlab('')
ggsave(gp, file=file.path(liana.dir,  paste0(ct, '.pdf')), width = 30, height = 20, units = 'cm' )
}

gp<-ggplot(lian_sub[source=='OL' & (ligand.complex %in% c('TF', 'SPP1', 'CNTN2', 'APOD') | receptor.complex %in%  c('TF', 'SPP1', 'CNTN2', 'APOD'))], 
           aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('source', 'lesion'),scales = 'free')+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8))
gp
ggsave(gp, file=file.path(liana.dir,  paste0('OL_TF_CNTN2_APOD', '.pdf')), width = 30, height = 7, units = 'cm' )

gp<-ggplot(lian_sub[source=='OL' & (ligand.complex %in% c('TF', 'SP1') | receptor.complex %in%  c('TF', 'SP1'))], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('source', 'lesion'),scales = 'free')+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8))+ylab("")+xlab('')
gp
ggsave(gp, file=file.path(liana.dir,  paste0('OL_TF', '.pdf')), width = 25, height = 6, units = 'cm' )



gp<-ggplot(lian_sub[source=='IM' & (ligand.complex == 'C3' | receptor.complex == 'C3')], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('source', 'lesion'),scales = 'free')+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')

ggsave(gp, file=file.path(liana.dir,  paste0('IM_C3.pdf')), width = 25, height = 7.5, units = 'cm' )




liana_test_C4_WM <- liana_objects[['C4_WM']] %>%
  liana_aggregate()

dp<-liana_test_C4_WM %>% liana_dotplot(source_groups = c("Oligodendrocytes"),
                             target_groups = c("Immune cells", "Astrocytes", "OPCs"), ntop=20)

liana_test_C4_WM %>% liana_dotplot(source_groups = c("Oligodendrocytes"),
                                   target_groups = c("Immune cells", "Astrocytes", "OPCs"), ntop=20)

liana_test_C4_WM <- liana_test_C4_WM %>%
  # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected



p <- chord_freq(liana_test_C4_WM,
                source_groups = c("Immune cells", "Astrocytes", "OPCs", 'Oligodendrocytes'),
                target_groups = c("Immune cells", "Astrocytes", "OPCs", 'Oligodendrocytes'))



