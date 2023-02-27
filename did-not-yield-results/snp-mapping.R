setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell-new')
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
#require(biovizBase)
require(wesanderson)
require(biomaRt)
require(gtable)
require(biomaRt)
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
require(GenomicRanges)

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
vizualisation.dir<-file.path(out.dir, 'track_viz_SNP')
dir.create(vizualisation.dir)

hm.integrated<-readRDS( file = file.path(out.dir, "integrated.rds"))

hm.integrated$meta<-as.ordered(paste0(hm.integrated$predicted_celltype, '.', hm.integrated$is_ms))

celltype_color = as.character(wes_palette('GrandBudapest1', 4))
celltype_color  = c('Oligodendrocytes.MS'=celltype_color[1], 'Astrocytes.MS'=celltype_color[2], 'Immune cells.MS'=celltype_color[3], 'OPCs.MS'=celltype_color[4],
                    'Oligodendrocytes.C'=celltype_color[1], 'Astrocytes.C'=celltype_color[2], 'Immune cells.C'=celltype_color[3], 'OPCs.C'=celltype_color[4])


hm.integrated@active.assay<-'peaks'
Idents(hm.integrated)<-'meta'
annotations <- Annotation(object = hm.integrated)


atac.de.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/differential_celltypes')
ct<-'Oligodendrocytes'
#de_acc<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-atac/all_peaks/outs/peaks.bed', sep = '\t', skip = '# prim', header=F)

### gRCH38
de_acc<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-atac/all_peaks/outs/peaks.bed', sep = '\t', skip = '# prim', header=F)

peak_ranges<-GRanges(seqnames = de_acc$V1,
                     ranges = IRanges(start = de_acc$V2,
                                      end = de_acc$V3,
                                      names = de_acc$V4))

snps_all<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-atac/snps_vep.txt', sep = '\t')
snps_all<-unique(snps_all)
snps_all<-as.data.table(snps_all)
snps_all[, c('chrom','start:end') := tstrsplit(Location, ':')]
snps_all[, c('start','end') := tstrsplit(`start:end`, '-')]
snps_all<-unique(snps_all[, c('#Uploaded_variation', 'chrom', 'start', 'end')])

snp_ranges<-GRanges(seqnames = paste0('chr', snps_all$chrom),
                    ranges = IRanges(start = as.numeric(snps_all$start),
                                     end = as.numeric(snps_all$end),
                                     names = snps_all$`#Uploaded_variation`))

grl<-list()
for (i in 1: length(snp_ranges)){
  rang<-GenomicRanges::intersect(peak_ranges, snp_ranges[i], ignore.strand = TRUE)
  rang<-as.data.table(rang)
  rang$mysnp<-snp_ranges[i]@ranges@NAMES
  grl[[snp_ranges[i]@ranges@NAMES]]<-rang
}

grl<-rbindlist(grl)
myrang2<-as.data.table(grl)



myrang2$seqnames<-stringr::str_replace(myrang2$seqnames, 'chr', '')
myrang2$seqnames<-as.numeric(myrang2$seqnames)

snps_all2<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-atac/snps_vep.txt', sep = '\t')
snps_all2<-unique(snps_all2)
snps_all2<-as.data.table(snps_all2)
snps_all2[, c('chrom','start:end') := tstrsplit(Location, ':')]
snps_all2[, c('start','end') := tstrsplit(`start:end`, '-')]

gn2<-unique(snps_all2[chrom %in% myrang2$seqnames & start %in% myrang2$start]$SYMBOL)

gn.sv<-data.table(gene=gn2)
fwrite(gn.sv, file = file.path(state.dir, 'open_gene.tsv'), sep= '\t')

adata$meta<-as.ordered(paste0(adata$celltypes_annotated, '.', adata$ms))
Idents(adata)<-'meta'

for (ge in gn2){
  print(ge)
  #'ELMO1',  

  if (ge %in% c('DIPK1A' ,'NA', '-' , 'YPEL3-DT', 'ZNF767P')){
    print(ge)
    print('NA')
    next
  }

  isgene <- annotations$gene_name == ge
  isgene <- !is.na(x = isgene) & isgene
  annot.sub <- annotations[isgene]
  gr <- GRanges(seqnames = as.character(x = seqnames(x = annot.sub))[[1]],
                ranges = IRanges(start = min(start(x = annot.sub)),
                                 end = max(end(x = annot.sub))))
  
  
  
cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = gr,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
)+scale_fill_manual(values = celltype_color)


expr_plot <- ExpressionPlot(
  object = adata,
  features = ge,
  assay = "RNA",
  idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
)+scale_fill_manual(values = celltype_color)


gene_plot <- AnnotationPlot(
  object = hm.integrated,
  region = gr,
)

p<-CombineTracks(
  plotlist = list(cov_plot, gene_plot),
  expression.plot = expr_plot,
  widths = c(10, 2)
)

ggsave(p, file=file.path(vizualisation.dir, paste0(ge, '.pdf')), width = 20, height = 20, units = 'cm')
}



# for (n in 1:length(ranges)){
#   ge <-ranges[n]
#   print(ge)
#   cov_plot<-CoveragePlot(
#     object = hm.integrated,
#     region = ge,
#     annotation = FALSE,
#     peaks = FALSE,
#     idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
#   )+scale_fill_manual(values = celltype_color)
#   
#   print('Annotation plot')
#   gene_plot <- AnnotationPlot(
#     object = hm.integrated,
#     region = ge,
#   )
#   
#   p<-CombineTracks(
#     plotlist = list(cov_plot, gene_plot),
#     widths = c(10, 2)
#   )
#   ggsave(p, file=file.path(vizualisation.dir, paste0(ge, '.pdf')), width = 20, height = 20, units = 'cm')
# }



####
# Expression of genes with SNPs 
filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/', 'filtering')
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))

state.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/specific_expression_SNPs/')
dir.create(file.path(state.dir, 'plots'), recursive = T)

# A) Overlap by rsID
mapSNPstoGenes<-unique(snps_all2[, c('#Uploaded_variation', 'SYMBOL')])
mapSNPstoGenes<-unique(mapSNPstoGenes)
fwrite(mapSNPstoGenes, file = file.path(state.dir, 'map_SNPs_to_genes.tsv'), sep= '\t')

# B) Make plots for all genes
Idents(adata)<-'celltypes_annotated'
for (t in unique(gn2[gn2 %in% rownames(adata@assays$RNA@meta.features) ])){

    CD44<-VlnPlot(object = adata, features = t, group.by = 'celltypes_annotated', split.by = 'patient', 
                 idents= c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs'))+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+ylim(0.00, 4)+xlab('')+ylab("Expression Level (0 truncated)")
  CD44
  ggsave(CD44, file=file.path(state.dir, 'plots', paste0(t, '.pdf')), width=40, height = 12, units = 'cm')
  
}


celltypes<-c("Oligodendrocytes", 'Immune cells', 'Astrocytes', 'OPCs')
for (ct in celltypes){
  lesion.de.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/', ct)
  deres<-fread(file.path(lesion.de.dir, paste0('DE_genes_', ct, '.tsv')), sep='\t')

  reg<-deres[gene %in% z$hgnc_symbol] 
  reg<-merge(reg, mapSNPstoGenes, by.x = 'gene', by.y = 'hgnc_symbol')
  reg<-unique(reg)
  fwrite(reg, file = file.path(state.dir, paste0(ct, '_DE_genes_associated_with_SNP.txt')), sep='\t')
}
