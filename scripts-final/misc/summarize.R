setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')

renv::activate('r-single-cell')
renv::repair()
set.seed(11)

require(data.table)
require(readxl)
require(tidyr)
require(dplyr)
require(stringr)
rna.results.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/')


lesions<-c('AL', 'NAWM', 'CA', 'RL')

## read and summarize bulk data
bulk_list<-list()
for (l in lesions){
  bulk<-read_excel(paste0('/work/MariaLouiseElkjær#7716/Bulk_RNA/WM_vs_', l, '.xlsx'))
  bulk<-as.data.table(bulk)
  bulk$lesion<-l
  bulk_list[[l]]<-bulk
  
}
bulk<-rbindlist(bulk_list)
colnames(bulk)<-c('ensemble', 'gene_name', 'logFC', 'FDR', 'lesion')
bulk<-bulk[FDR<0.01 & abs(logFC)>0.5]
bulk<-as.data.table(bulk)
bulk$celltype<-'bulk'

#BULK FOLD CHANGE
bulk<-bulk %>% pivot_wider(id_cols = c(ensemble, gene_name), values_from = logFC, names_from = c(celltype, lesion))


## CELLTYPE_LESION TYPE FOLD CHANGE
celltypes<-c('Astrocytes', 'Oligodendrocytes', 'Immune cells', 'OPCs')
single_list<-list()


for (ct in celltypes){
lesion.de.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/', ct)
snAstro<- read_gene_list(paste0('DE_genes_',ct,'.tsv'), lesion.de.dir)
snAstro$celltype <- ct
single_list[[ct]]<-snAstro
}


single<-rbindlist(single_list)
single<-as.data.table(single)
single<-single %>% pivot_wider(id_cols = c(gene), values_from = logFC, names_from = c(lesion, celltype))
all<-merge(single, bulk, by.x = 'gene', by.y = 'gene_name', all = T)


# FOLD CHANGE ACCESSIBILITY ATAC IN CODING GENES

atac_list<-list()
lesion.de.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-atac/differential_celltypes')


for (ct in celltypes){
  atac<- read_gene_list(paste0(ct,'_DE_accessible_region_annot.tsv'), lesion.de.dir)
  atac$celltype <- ct
  atac_list[[ct]]<-atac
  
}
atac<-rbindlist(atac_list)
atac<-as.data.table(atac)
atac<-atac[distance<5000 & p_val_adj<0.01]

atac<- atac %>%
  group_by(gene_name, celltype) %>%
  summarise(avg_log2FC := mean(avg_log2FC)) 
atac$lesion<-'MS'

atac<- atac %>% pivot_wider(id_cols = c(gene_name), values_from = avg_log2FC, names_from = c(lesion, celltype))

all<-merge(all, atac, by.x = 'gene', by.y = 'gene_name', all=T)


## AVG FOLD CHANGE CA RL in spatial clustering.
sample_name<-c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A', 'RL_1B', 'RL_1C', 'RL_1D')
de_list<-list()
for(sn in sample_name){
  name<-paste0(sn, '_differential_genes.tsv')
  data <-fread(file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/plots/',sn, name))
  data[,lesion:=stringr::str_split(sn, '_')[[1]][1]]
  de_list[[sn]]<-data
}
data<-rbindlist(de_list)
data<-data[p_val_adj<0.01]
data<- data %>%
  group_by(gene, lesion) %>%
  summarise(avg_log2FC := mean(avg_log2FC)) 
data$celltype<-'spatial'
data<- data %>% pivot_wider(id_cols = c(gene), values_from = avg_log2FC, names_from = c(lesion, celltype))

all<-merge(all, data, by.x = 'gene', by.y = 'gene', all=T)



all<-all %>% select(-ensemble)
all<-as.data.table(all)
all$evidence<-apply(all,1, function(x) sum(!is.na(x)))-1
all<-all[order(evidence, decreasing = TRUE)]

rna.results.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/summary')
dir.create(rna.results.dir)
fwrite(all, file.path(rna.results.dir, 'summary_genes.tsv'), sep='\t')
