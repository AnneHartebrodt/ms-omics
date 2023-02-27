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


out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
celltype.de.dir<-file.path(out.dir, 'differential_celltypes')


hm.integrated<-readRDS( file = file.path(out.dir, "integrated.rds"))

############## MOTIF ANALYSIS #########################

### PREPROCESSING: READ AND FETCH DATA
hm.intgrated<-RegionStats(hm.integrated, genome = BSgenome.Hsapiens.UCSC.hg38)
pwm <- getMatrixSet(
  x = JASPAR2022,
  opts = list(species = 9606, all_versions = FALSE)
)
# remove scaffolds
main.chroms <- standardChromosomes(BSgenome.Hsapiens.UCSC.hg38)
keep.peaks <- which(as.character(seqnames(granges(hm.integrated))) %in% main.chroms)
hm.integrated[["ATAC"]] <- subset(hm.integrated[["peaks"]], features = rownames(hm.integrated[["peaks"]])[keep.peaks])
hm.integrated <- AddMotifs(hm.integrated, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm, verbose = T, assay = 'ATAC')


### CANDIDATE TRANSCRIPTION FACTORS
tfs<- fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/annotation/transcription_factors.tsv', header = FALSE)
mytfs<- tfs[tfs$V1 %in% as.character(lapply(pwm@listData, function(x) x@name))]



#### CHECK THE CANDIDATE TFS FOR CELLTYPE SPECIFIC ACCESSIBILITY
### CELLTYPE SPECIFIC AND MS VS. CTRL FOR THE WHOLE DATA

## EXTRACT THE MOTIFS FROM THE PWM. NOT ALL CANDIDATES ARE CONTAINED IN THE PWM
n<-as.character(lapply(pwm@listData, function(x) x@ID))[which(as.character(lapply(pwm@listData, function(x) x@name)) %in% tfs$V1)]
m<-as.character(lapply(pwm@listData, function(x) x@name))[which(as.character(lapply(pwm@listData, function(x) x@name)) %in% tfs$V1)]
for (tf in n){
    hm.integrated <- Footprint(
    object = hm.integrated,
    assay = 'ATAC',
    motif.name = tf,
    genome = BSgenome.Hsapiens.UCSC.hg38,
    in.peaks=TRUE
  )
}

footprint.dir<-file.path(out.dir, 'footprints')
dir.create(footprint.dir)
for (tf in 1:length(n)){
  Idents(hm.integrated)<-'predicted_celltype'
  fp1 <- PlotFootprint(hm.integrated, features = n[tf], assay = 'ATAC', idents = c("Astrocytes", 
                                                                                'Oligodendrocytes', 'Immune cells', 'OPCs'))
  ggsave(fp1, file=file.path(footprint.dir, paste0(n[tf], '_', m[tf], '_celltypes.pdf')), width = 15, height = 10, units = 'cm')
  ggsave(fp1, file=file.path(footprint.dir, paste0(n[tf], '_', m[tf], '_celltypes.png')), width = 15, height = 10, units = 'cm')
  
  Idents(hm.integrated)<-'is_ms'
  fp2<- PlotFootprint(hm.integrated, features = n[tf], assay = 'ATAC')
  
  ggsave(fp2, file=file.path(footprint.dir, paste0(n[tf], '_', m[tf], '_ms_ctrl_all.pdf')), width = 15, height = 10, units = 'cm')
  ggsave(fp2, file=file.path(footprint.dir, paste0(n[tf], '_', m[tf], '_ms_ctrl_all.png')), width = 15, height = 10, units = 'cm')
}



## CHECK THE CANDIDATE FOOTPRINTS FOR ENRICHMENT WITHIN THE CELLTYPES
### CELLTYPE SPECIFIC MS. VS CONTROL
celltypes<- c( 'Astrocytes', 'Immune cells', 'OPCs','Oligodendrocytes')
Idents(hm.integrated)<-'is_ms'
for(ct in  celltypes){
    sub<-subset(hm.integrated, subset = predicted_celltype==ct)
    footprint.dir.ct<-file.path(footprint.dir, ct)
    dir.create(footprint.dir.ct)
    for (tf in 1:length(n)){
      sub <- Footprint(
      object = sub,
      assay = 'ATAC',
      motif.name = n[tf],
      genome = BSgenome.Hsapiens.UCSC.hg38,
      in.peaks=TRUE
    )
    
  fp2 <- PlotFootprint(sub, features = n[tf], assay = 'ATAC')
  ggsave(fp2, file=file.path(footprint.dir.ct, paste0(n[tf], '_', m[tf], '_ms_ctrl_all.pdf')), width = 15, height = 10, units = 'cm')
  ggsave(fp2, file=file.path(footprint.dir.ct, paste0(n[tf], '_', m[tf], '_ms_ctrl_all.png')), width = 15, height = 10, units = 'cm')
}
}
saveRDS(hm.integrated, file = file.path(out.dir, "integrated_tf.rds"))

hm.integrated <- readRDS(file.path(out.dir, "integrated_tf.rds"))


### USE THE DIFF ACC. REGIONS FROM THE DIFFERENTIAL EXPRESSION TO FIND MOTIFS
### CELLTYPE SPECIFIC MS. VS CONTROL
celltypes<-c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs')
DefaultAssay(hm.integrated)<-'ATAC'
for(ct in  celltypes){
  sub<-subset(hm.integrated, subset = predicted_celltype==ct)
  Idents(sub)<-'is_ms'
  de_genes<-fread(file.path(celltype.de.dir, paste0(ct, '_DE_accessible_region_annot.tsv')),sep='\t')
  enriched.motifs <- FindMotifs(
    object = sub,
    features = de_genes$region
  )
  enriched.motifs<-as.data.table(enriched.motifs)
  enriched.motifs<- enriched.motifs[p.adjust<0.05]
  fwrite(enriched.motifs, file.path(celltype.de.dir, paste0(ct, '_enriched_motifs.tsv')),sep='\t')
}


#### USE THE ENRICHED MOTIFS IDENTIFIED EALIER AND ANALYSE THE FOOTPRINTS
### CELLTYPE SPECIFIC MS. VS CONTROL
Idents(hm.integrated)<-'is_ms'
DefaultAssay(hm.integrated)<-'ATAC'
for(ct in  celltypes){
  sub<-subset(hm.integrated, subset = predicted_celltype==ct)
  Idents(sub)<-'is_ms'
  
  enriched.motifs<-fread(file.path(celltype.de.dir, paste0(ct, '_enriched_motifs.tsv')),sep='\t')
  celltype.de.dir.plots<-file.path(celltype.de.dir, ct)
  dir.create(celltype.de.dir.plots, showWarnings = FALSE)
  for (m in 1: length(enriched.motifs$motif)){
  sub <- Footprint(
    object = sub,
    assay = 'ATAC',
    motif.name = enriched.motifs$motif[m],
    genome = BSgenome.Hsapiens.UCSC.hg38,
    in.peaks=TRUE
  )
  fp1 <- PlotFootprint(sub, features = enriched.motifs$motif[m], assay = 'ATAC')
  ggsave(fp1, file=file.path(celltype.de.dir.plots, paste0(enriched.motifs$motif[m], '_', enriched.motifs$motif.name[m], '_ms_vs_ctrl.pdf')), width = 15, height = 10, units = 'cm')
  ggsave(fp1, file=file.path(celltype.de.dir.plots, paste0(enriched.motifs$motif[m], '_', enriched.motifs$motif.name[m], '_ms_vs_ctrl.png')), width = 15, height = 10, units = 'cm')
  }
}

hm.integrated <- readRDS(file.path(out.dir, "integrated_tf.rds"))
celltypes<- c( 'Astrocytes', 'Immune cells', 'OPCs','Oligodendrocytes')
footprint.dir<-file.path(out.dir, 'footprints-new-candidates')
not_found<-c('LUZP2',"SOX6","FOS", 'NFATCZ', 'RUNX1', "NFATCZ")
further_candidates<-c( "IKZF1","FOXP2", "BCL6", "CEBPD", "FLI1", "ETV6", "JUN", "JUNB", "TCF7L1", "ETS2", "NFATCZ")
Idents(hm.integrated)<-'is_ms'
DefaultAssay(hm.integrated)<-'ATAC'
for(ct in  celltypes){
  sub<-subset(hm.integrated, subset = predicted_celltype==ct)
  Idents(sub)<-'is_ms'
  
  enriched.motifs<-fread(file.path(celltype.de.dir, paste0(ct, '_enriched_motifs.tsv')),sep='\t')
  footprint.dir.ct<-file.path(footprint.dir, ct)
  dir.create(footprint.dir.ct, recursive = T)
  for (m in further_candidates){
    sub <- Footprint(
      object = sub,
      assay = 'ATAC',
      motif.name = m,
      genome = BSgenome.Hsapiens.UCSC.hg38,
      in.peaks=TRUE
    )
    fp1 <- PlotFootprint(sub, features = m, assay = 'ATAC')
    ggsave(fp1, file=file.path(footprint.dir.ct, paste0(m, '_ms_vs_ctrl.pdf')), width = 15, height = 10, units = 'cm')
    ggsave(fp1, file=file.path(footprint.dir.ct, paste0(m, '_ms_vs_ctrl.png')), width = 15, height = 10, units = 'cm')
  }
}

