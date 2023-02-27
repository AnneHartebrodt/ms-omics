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
require(future)

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')

qc.plot.dir<-file.path(out.dir, 'qc', 'plots')

dir.create(qc.plot.dir, recursive = TRUE)

hm.integrated.new<-readRDS( file = file.path(out.dir, "integrated.rds"))

prov<-hm.integrated$batch
levels(prov)
levels(prov)<-c('C', 'C', 'MS', 'MS', 'MS', 'MS', 'MS', 'MS', 'MS')
hm.integrated[['is_ms']]<- prov


####### DIFFERENTIAL EXPRESSION ANALYSIS
# change back to working with peaks instead of gene activities
# Do Differential expression of MS vs. control
DefaultAssay(hm.integrated) <- 'peaks'
Idents(hm.integrated)<-'predicted_celltype'

celltype.de.dir<-file.path(out.dir, 'differential_celltypes')
dir.create(celltype.de.dir)

# differentially accessible regions between celltypes
da_peaks <- FindAllMarkers(
  object = hm.integrated,
  min.pct = 0.1,
  #test.use = 'LR',
  #latent.vars = 'peak_region_fragments',
  logfc.threshold = 0.25
)

da_peaks$region<-rownames(da_peaks)
da_peaks<-as.data.table(da_peaks)
da_peaks<-da_peaks[abs(avg_log2FC)>0.25 & p_val_adj < 0.05]
fwrite(da_peaks, file.path(celltype.de.dir, paste0( 'celltypes_DE_accessible_region_all.tsv')),sep='\t')
relevant<- da_peaks[cluster %in% c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs') & (pct.1>0.25 | pct.2>0.25)]
fwrite(relevant, file.path(celltype.de.dir, paste0( 'celltypes_DE_accessible_region_filtered.tsv')),sep='\t')


closest <- ClosestFeature(hm.integrated, relevant$region)
closest<-as.data.table(closest)
fwrite(closest, file.path(celltype.de.dir, paste0('celltypes_closest_genomic_region_relevant.tsv')),sep='\t')


all<-merge(closest, relevant, by.x='query_region', by.y='gene')

tyeps<-ggplot(all, aes(x=gene_biotype, fill=type))+geom_bar()+facet_wrap(~cluster, nrow = 1)+theme(axis.text.x = element_text(angle=90))
ggsave(tyeps, file=file.path(celltype.de.dir, 'biotypes_genomic_region_relevant.pdf'))


######## DIFFERENTIALLY ACCESSIBLE REGIONS BETWEEN MS AND CTRL for all celltypes ########
# used celltypes
celltypes<- c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs')
Idents(hm.integrated)<-'is_ms'

# use find marker function for DA
for(ct in  celltypes){
  sub<-subset(hm.integrated, subset = predicted_celltype==ct)
  
  da_peaks <- FindMarkers(
    object = sub,
    min.pct = 0.1,
    #test.use = 'LR',
    #latent.vars = 'peak_region_fragments',
    ident.1 = 'MS',
    logfc.threshold = 0.25
  )
  
  da_peaks$region<-rownames(da_peaks)
  da_peaks<-as.data.table(da_peaks)
  da_peaks<-da_peaks[abs(avg_log2FC)>0.25 & p_val_adj < 0.05]
  fwrite(da_peaks, file.path(celltype.de.dir, paste0(ct, '_DE_accessible_region_all.tsv')),sep='\t')
  
}

# save lists and annotations
list_of_de<-list()
list_of_closest<-list()
for(ct in  celltypes){
  sub<-subset(hm.integrated, subset = predicted_celltype==ct)
  de_genes<- fread(file.path(celltype.de.dir, paste0(ct, '_DE_accessible_region_all.tsv')),sep='\t')
  closestr <- ClosestFeature(sub, de_genes$region)
  de_genes<-merge(de_genes, closestr, by.x = 'region', by.y = 'query_region')
  list_of_de[[ct]]<-de_genes
  fwrite(de_genes, file.path(celltype.de.dir, paste0(ct, '_DE_accessible_region_annot.tsv')),sep='\t')
  
}
