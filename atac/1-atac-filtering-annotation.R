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
#BiocManager::install("plger/scDblFinder")
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


out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
plot.dir<-file.path(out.dir, 'plots')
dir.create(plot.dir)

qc.plot.dir<-file.path(out.dir, 'qc', 'plots')
dir.create(qc.plot.dir, recursive = TRUE)

initial.plot.dir<-file.path(out.dir, 'initial')
dir.create(initial.plot.dir, recursive = TRUE)

rna.filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering/')


## read the data, already aggregated from cellranger atac
counts <- Read10X_h5(filename = "/work/AnneSophiaHartebrodt#6698/ms-project/ms-atac/all_peaks/outs/filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/work/AnneSophiaHartebrodt#6698/ms-project/ms-atac/all_peaks/outs/singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = '/work/AnneSophiaHartebrodt#6698/ms-project/ms-atac/all_peaks/outs/fragments.tsv.gz',
  min.cells = 10,
  min.features = 50
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  project='ATAC',
  meta.data = metadata
)




# run profile
## Get annotations
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
#keepStandardChromosomes(annotations, 'coarse')

# add the gene information to the object
Annotation(pbmc) <- annotations

# Add sample name
sum<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/ms-atac/all_peaks/outs/aggregation_csv.csv')
prov<-as.character(sapply(names(pbmc$orig.ident), function(x) paste0(strsplit(x, '-')[[1]][2], collapse='_')))
prov<-as.factor(prov)
levels(prov)<-sum$library_id
pbmc[['batch']]<-prov

levels(prov)
levels(prov)<-c('C', 'C', 'MS', 'MS', 'MS', 'MS', 'MS', 'MS', 'MS')
pbmc[['is_ms']]<- prov


saveRDS(pbmc, file = file.path(out.dir, 'initial.rds'), compress=FALSE)

# manually checked with web_summary.html if samples are correctly assigned!
dbllist<-list()
for(b in unique(pbmc$batch)){
  sub <- subset(
    x = pbmc,
    subset = batch == b )
  sce<-as.SingleCellExperiment(sub)
  # identify doublets
  sce <- scDblFinder(sce, artificialDoublets=1, aggregateFeatures=TRUE, nfeatures=30, processing="normFeatures", iter = 3)
  df<- data.table(sce@colData@rownames, sce$scDblFinder.class, sce$scDblFinder.score)
  dbllist[[b]]<-df
}

## Annotate doublets
doublet_calls<- rbindlist(dbllist)
colnames(doublet_calls)<- c('barcode', 'prediction', 'score')
meta<-pbmc@meta.data
meta$barcode<-rownames(meta)
mm<-merge(meta, doublet_calls, by.x='barcode', by.y = 'barcode')
pbmc@meta.data<-mm
rownames(pbmc@meta.data)<-pbmc@meta.data$barcode


## Signac QC pipeline
# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(object = pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)


# not 
pbmc$blacklist_fraction <- FractionCountsInRegion(
  object = pbmc, 
  assay = 'peaks',
  regions = blacklist_hg19
)

# add blacklist ratio and fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments

# plot TSS enrichment for batch and overall
pbmc$high.tss <- ifelse(pbmc$TSS.enrichment > 2, 'High', 'Low')
tss<-TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
tss2<-TSSPlot(pbmc, group.by = 'batch') + NoLegend()
ggsave(tss, file= file.path(qc.plot.dir, paste0('tss_enrichment.png')))
ggsave(tss2, file= file.path(qc.plot.dir, paste0('tss_enrichment.batch.png')))

#tss
#tss2

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
fgm1<-FragmentHistogram(object = pbmc, group.by = 'batch')
fgm2<-FragmentHistogram(object = pbmc, group.by = 'nucleosome_group')
ggsave(fgm1, file= file.path(qc.plot.dir, paste0('fragment_histogram.batch.png')))
ggsave(fgm2, file= file.path(qc.plot.dir, paste0('fragment_histogram.nucleosome.group.png')))

#pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

#VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2, group.by = "batch", combine = FALSE)

diagnostic<-VlnPlot(
  object = pbmc,
  features = c('pct_reads_in_peaks', 'peak_region_fragments',
               'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5
)
features = c('pct_reads_in_peaks', 'peak_region_fragments',
             'TSS.enrichment', 'blacklist_fraction', 'nucleosome_signal')

for(i in 1:length(diagnostic)){
  ggsave(diagnostic[[i]], file= file.path(qc.plot.dir, paste0(features[i], '.png')))
  
}


pbmc <- subset(
  x = pbmc,
  subset = peak_region_fragments > 3000 &
    peak_region_fragments < 20000 &
    pct_reads_in_peaks > 15 &
    blacklist_fraction < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 2 &
    prediction == 'singlet'
)
pbmc

DefaultAssay(pbmc) <- 'peaks'

# Normalize, feature selection, SVD, UMAP, neighbors and clustering
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:30)
pbmc <- FindClusters(object = pbmc, verbose = TRUE)


batch<-pbmc[['batch']]$batch
levels(batch)<-sum$library_id
pbmc[['batch']]<-batch


# Look at the embedding
umap<-DimPlot(object = pbmc) 
umap
ggsave(umap, file= file.path(initial.plot.dir, paste0('initial_umap.png')))

### LOAD MARKER GENE LIST

markers<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/annotation/marker_genes.txt')
# Look at some gene activities to confirm that the batch effect is indeed very
# strong
gene.activities <- GeneActivity(pbmc)
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

FeaturePlot(
  object = pbmc,
  slot = 'counts',
  features = markers[Cell_type=='Oligodendrocytes']$Gene_name,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)


saveRDS(pbmc, file = file.path(out.dir, 'prep.rds'), compress=FALSE)




### INTEGRATION USING HARMONY:
hm.integrated <- RunHarmony(object = pbmc, group.by.vars = 'batch', reduction = 'lsi', assay.use = 'peaks', project.dim = FALSE)
hm.integrated <- RunUMAP(object = hm.integrated, dims = 2:30, reduction='harmony')
hm.integrated <- FindNeighbors(object = hm.integrated, reduction = 'lsi', dims = 3:30)
hm.integrated <- FindClusters(object = hm.integrated, verbose = FALSE, algorithm = 3)



## LOAD THE CLEANED AND ANNOTATED SINGLE RNA SEQ OBJECT
adata<- readRDS( file = file.path(rna.filtering.dir, "annotated.rds"))

## TRANSFER THE LABELS
transfer.anchors <- FindTransferAnchors(reference = adata, query = hm.integrated, features = VariableFeatures(object = adata),
                                        reference.assay = "RNA", query.assay = "RNA", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = adata$celltypes_annotated,
                                     weight.reduction = hm.integrated[['harmony']], dims = 2:30)

##ct.pred<-celltype.predictions$predicted.id
#names(ct.pred)<-rownames(celltype.predictions)
#hm.integrated[["RNA"]] <- celltype.predictions
hm.integrated$predicted_celltype <- celltype.predictions$predicted.id

celltype_color = c(as.character(wes_palette('GrandBudapest1', 4)), as.character(wes_palette('GrandBudapest2', 4)))
celltype_color  = c('Oligodendrocytes'=celltype_color[1], 'Astrocytes'=celltype_color[2], 'Immune cells'=celltype_color[3], 'OPCs'=celltype_color[4],
                    'Neurons' = celltype_color[5], 'Endothelial cells' = celltype_color[6], 'Pericytes'=celltype_color[7], 'B-cells'=celltype_color[8])

p1<-DimPlot(adata, group.by='celltypes_annotated', cols = celltype_color)+ggtitle("Celltype annotation (RNASeq)")+ NoLegend() 
p2 <- DimPlot(hm.integrated, group.by = "predicted_celltype", label = FALSE, cols=celltype_color) +  ggtitle("Predicted celltypes (ATACSeq)")
p3 <- DimPlot(adata, group.by = "batch", label = FALSE)  + ggtitle("Batches (RNASeq)")
p4 <- DimPlot(hm.integrated, group.by = "batch", label = FALSE)  + ggtitle("Batches (ATACSeq)")

pp<-p1+p2
ggsave(pp, file =file.path(plot.dir, paste0('celltypes_atac_rna', '.png')), width = 23, height = 10, units = 'cm')
ggsave(pp, file =file.path(plot.dir, paste0('celltypes_atac_rna', '.pdf')), width = 23, height = 10, units = 'cm')




astroatac<- FeaturePlot(
  object = hm.integrated,
  features = markers[Cell_type=='Astrocytes']$Gene_name,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
astroatac
ggsave(astroatac, file =file.path(plot.dir, paste0('astrocytes', '.png')), width = 23, height = 10, units = 'cm')
ggsave(astroatac, file =file.path(plot.dir, paste0('astrocytes', '.pdf')), width = 23, height = 10, units = 'cm')


olicoatac<-FeaturePlot(
  object = hm.integrated,
  features = markers[Cell_type=='Oligodendrocytes']$Gene_name,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
ggsave(olicoatac, file =file.path(plot.dir, paste0('oligodendrocytes', '.png')), width = 23, height = 10, units = 'cm')
ggsave(oligoatac, file =file.path(plot.dir, paste0('oligodendrocytes', '.pdf')), width = 23, height = 10, units = 'cm')



opcatac<-FeaturePlot(
  object = hm.integrated,
  slot = 'data',
  features = markers[Cell_type=='OPCs']$Gene_name,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
ggsave(opcatac, file =file.path(plot.dir, paste0('opcs', '.png')), width = 23, height = 10, units = 'cm')
ggsave(opcatac, file =file.path(plot.dir, paste0('opcs', '.pdf')), width = 23, height = 10, units = 'cm')


immuneatac<-FeaturePlot(
  object = hm.integrated,
  slot = 'counts',
  features = markers[Cell_type %in% c('Microglia-Macrophages', 'T-cells', 'B-cells', 'Microglia-resting')]$Gene_name,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
ggsave(immuneatac, file =file.path(plot.dir, paste0('immune', '.png')), width = 23, height = 10, units = 'cm')
ggsave(immuneatac, file =file.path(plot.dir, paste0('immune', '.pdf')), width = 23, height = 10, units = 'cm')

plot.dir.ind<-file.path(plot.dir, 'umaps')
dir.create(plot.dir.ind)
for(gene in c(markers$Gene_name,c('ERG', 'VWF', 'EMCN', 'CLDN5', 'EGFL7'))){
  tcell<-FeaturePlot(hm.integrated, features = gene)+NoLegend()
  ggsave(tcell, file =file.path(plot.dir.ind, paste0(gene, '.png')), width = 10, height = 10, units = 'cm')
  ggsave(tcell, file =file.path(plot.dir.ind, paste0(gene, '.pdf')), width = 10, height = 10, units = 'cm')
}



saveRDS(hm.integrated, file = file.path(out.dir, "integrated.rds"), compress=FALSE)
#hm.integrated<-readRDS( file = file.path(out.dir, "integrated.rds"))


