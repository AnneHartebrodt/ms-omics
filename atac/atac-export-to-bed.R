setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(data.table)
require(stringr)

atac.de.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/differential_celltypes')


for (ct in c('Oligodendrocytes', 'Astrocytes', 'Immune cells', 'OPCs')){
atac<-fread(file.path(atac.de.dir, paste0(ct, '_DE_accessible_region_all.tsv')), sep = '\t')


atac[, c('chrom','start', 'stop') := tstrsplit(region, "-", fixed=TRUE)]
atac<-atac[,.(chrom, start, stop)]

fwrite(atac, file.path(atac.de.dir, paste0(ct, '_DE_accessible_region_all.bed')), sep = '\t', col.names = F)
}
