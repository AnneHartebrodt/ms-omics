setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell-new')
renv::repair()
set.seed(11)


require(tidyverse)
require(data.table)
require(ggplot2)
require(Seurat)
require(Signac)
require(wesanderson)


filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/', 'filtering')
adata<- readRDS(file.path(filtering.dir, "annotated.rds"))


variant.dir<-'/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/variant_calling/'
snps_all<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/data/ms-atac/snps_vep.txt', sep = '\t')
snps_all<-unique(snps_all)
colnames(snps_all)[1]<-'variant_id'
glist<-list()

samples<-c('MS1', 'MS3', 'MS6', 'MS7', 'MS9', 'C1', 'C3', 'C4', 'C1_13scattered', 'C1_15scattered', 'MS3_9scattered', 'MS3_10scattered')
for (ms in samples){
  MS1.1<-fread(file.path(variant.dir, ms, paste0(ms, '.haplotypecaller.filtered.vcf.gz')),  skip = '#CHROM')
  GWAS<-MS1.1[ID %in% snps_all$variant_id]
  GWAS$patient<-ms
  glist[[ms]]<-GWAS
}


GWAS<-rbindlist(glist)

annotated<-merge(GWAS, snps_all, by.x = 'ID', by.y = c('variant_id'))
annotated$patient<-recode(annotated$patient, 'MS3_10scattered'='MS3', 'MS3_9scattered'='MS3', 'C1_13scattered'='MS1', 'C1_15scattered'='MS1')


state.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/specific_expression_SNPs/')
mapSNPsTogenes<-fread(file = file.path(state.dir, 'map_SNPs_to_genes.tsv'), sep= '\t')

gn<-fread(file = file.path(state.dir, 'open_gene.tsv'), sep= '\t')


snp.dir<-'/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/snps'
dir.create(snp.dir)
snp.dir.open<-file.path(snp.dir, 'open_transcribed_and_snp')
dir.create(snp.dir.open)


open.transcribed.and.snp<-annotated[SYMBOL %in% gn$gene]
open.transcribed.and.snp<-open.transcribed.and.snp[SYMBOL!='-']
summary<-unique(open.transcribed.and.snp[, c('ID', 'patient', 'SYMBOL', 'REF', 'ALT', 'FORMAT', 'MS1')])

fwrite(summary, file=file.path(snp.dir.open, 'summary_atac_open_rna_transcribed_snp_called.tsv'), sep='\t')
fwrite(open.transcribed.and.snp, file=file.path(snp.dir.open, 'all_atac_open_rna_transcribed_snp_called.tsv'), sep='\t')

require(Seurat)
Idents(adata)<-'celltypes_annotated'


for(mygene in unique(summary$SYMBOL[summary$SYMBOL %in% rownames(adata@assays$RNA)])){
  pl<-VlnPlot(adata, features = mygene, split.by = 'patient')
  ggsave(pl, file=file.path(snp.dir.open, paste0(mygene, '.pdf')))
}


snp.dir.found<-file.path(snp.dir, 'snp_found')
dir.create(snp.dir.found)
summary.all<-unique(annotated[, c('ID', 'patient', 'SYMBOL', 'REF', 'ALT', 'FORMAT', 'MS1')])
summary.all<-summary.all[SYMBOL !='-']

fwrite(summary.all, file.path(snp.dir.found, 'all_snps.tsv'), sep='\t')


for(mygene in unique(summary.all$SYMBOL[summary.all$SYMBOL %in% rownames(adata@assays$RNA)])){
  pl<-VlnPlot(adata, features = mygene, split.by = 'patient')
  ggsave(pl, file=file.path(snp.dir.found, paste0(mygene, '.pdf')), width=30, height = 10, units = 'cm')
}

glist<-list()
samples<-c('MS1', 'MS3', 'MS6', 'MS7', 'MS9', 'C1', 'C3', 'C4', 'C1_13scattered', 'C1_15scattered', 'MS3_9scattered', 'MS3_10scattered')
for (ms in samples){
  MS1.1<-fread(file.path(variant.dir, ms, paste0(ms, '.haplotypecaller.filtered.vcf.gz')),  skip = '#CHROM')
  MS1.1$patient<-ms
  glist[[ms]]<-MS1.1
}


GWAS<-rbindlist(glist)
GWAS<-merge(GWAS, snps_all, by.x = 'ID', by.y = c('variant_id'))

colnames(GWAS)[2]<-'chrom'
ug<-unique(GWAS[, c('ID', 'patient', 'SYMBOL')])
ug$patient<-recode(ug$patient, 'MS3_10scattered'='MS3', 'MS3_9scattered'='MS3', 'C1_13scattered'='MS1', 'C1_15scattered'='MS1')

all.called<-ug[, .N, by=c('ID', 'SYMBOL')][N==8]



summary.all<-fread(file.path(snp.dir.found, 'all_snps.tsv'), sep='\t')

ct<-'Oligodendrocytes'

for (ct in c('Astrocytes', 'Oligodendrocytes', 'Immune cells', 'OPCs')){

de.dir<-'/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/'
de.genes<-fread(file.path(de.dir, ct, paste0('DE_genes_', ct, '.tsv')), sep='\t')
print(de.genes[gene %in% unique(summary.all$SYMBOL)])
}

summary.all[SYMBOL=='WWOX']

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

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac')
vizualisation.dir<-file.path(out.dir, 'track_viz_SNP')
dir.create(vizualisation.dir)

adata$meta<-as.ordered(paste0(adata$celltypes_annotated, '.', adata$ms))
Idents(adata)<-'meta'

gr = 'WWOX'
cov_plot<-CoveragePlot(
  object = hm.integrated,
  region = gr,
  annotation = FALSE,
  peaks = FALSE,
  idents = c('Astrocytes.MS', 'Astrocytes.C', 'Oligodendrocytes.MS','Oligodendrocytes.C', 'Immune cells.MS', 'Immune cells.C', 'OPCs.MS', 'OPCs.C')
)+scale_fill_manual(values = celltype_color)


Idents(adata)<-'meta'
expr_plot <- ExpressionPlot(
  object = adata,
  features = gr,
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

ggsave(p, file=file.path(vizualisation.dir, paste0(gr, '.pdf')), width = 20, height = 20, units = 'cm')

