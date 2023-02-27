setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(data.table)
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

#BiocManager::install(organism, character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
require(DOSE)
library("pathview")

require(ggvenn)

require(nVennR)

rna.differential<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR/')
enrichment.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/enrichment')
dir.create(enrichment.dir)

lesions<-c('CA', 'AL', 'NAWM', 'RL')
celltypes<-c("Astrocytes", 'Oligodendrocytes', 'Immune cells', 'OPCs')


# 1. do enrichemnt analysis for each celltype separately
for (ct in celltypes){
  go_enrichment_celltypes<-file.path(enrichment.dir, 'GO',ct)
  dir.create(go_enrichment_celltypes, recursive = T)
  lesion.de.dir<-file.path(rna.differential ,ct)
  
  
  deres<-fread(file.path(lesion.de.dir, paste0('DE_genes_', ct, '.tsv')), sep='\t')
  deres[,N:=.N, by=gene]
  #deres<-deres[N>=2]
  
  # SET THE DESIRED ORGANISM HERE
  organism = "org.Hs.eg.db"
  
  eg <- bitr(deres$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  deres<-merge(deres, eg, by.x = 'gene', by.y = 'SYMBOL')
  
  all_de_names<-deres$gene
  all_de_fc<-deres$logFC
  all_de_id<-deres$ENTREZID
  names(all_de_fc)<-all_de_id
  
  # perform GO overrepresentation analysis
  for (o in c('BP', 'MF', 'CC')){
  ego <- enrichGO(gene          = all_de_names,
                  OrgDb         = 'org.Hs.eg.db',
                  ont           = o,
                  minGSSize = 15, 
                  maxGSSize = 500, 
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.05,
                  pAdjustMethod = 'fdr',
                  readable      = TRUE,
                  keyType = "SYMBOL"
  )
  rego<-as.data.table(ego@result)
  rego<- rego[p.adjust<0.05]
  if(nrow(rego)>0){
  fwrite(rego,file.path(go_enrichment_celltypes, paste0('GO_',o,'_Terms_', ct, '.tsv')), sep='\t')
  }
  }
  
  ## COMBINE ALL LESION TYPES FOR INITAL PATHWAY ANALYSIS
  #### make KEGG pathways for all cell types
  
  candidates<-fread('/work/AnneSophiaHartebrodt#6698/ms-project/submission/kegg_candidates.tsv', header = F)
    
  kegg.dir <- file.path(enrichment.dir, 'kegg/', ct)
  dir.create(kegg.dir, recursive = T, showWarnings = F)
  setwd(kegg.dir)
  for (p in candidates$V1){
    hsa04110 <- pathview(gene.data  = all_de_fc,
                          pathway.id = p,
                          limit      = list(gene=max(abs(all_de_fc)), cpd=1),
                          kegg.dir = enrichment.dir)
  }
  
  kegg_enrichment_celltypes<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/enrichment/kegg/putative', ct)
  dir.create(kegg_enrichment_celltypes, recursive = T)

  
  kk2 <- enrichKEGG(gene = all_de_id,
              organism     = 'hsa',
              pvalueCutoff = 0.05)
  
  kk2<-as.data.table(kk2)
  fwrite(kk2, file.path(kegg_enrichment_celltypes, paste0(ct, '_putative_KEGG_pathways.tsv')), sep='\t')
  
  setwd(kegg_enrichment_celltypes)
  for (p in kk2$ID){
    hsa04110 <- pathview(gene.data  = all_de_fc,
                         pathway.id = p,
                         species    = "hsa",
                         limit      = list(gene=max(abs(all_de_fc)), cpd=1),
                         kegg.dir = enrichment.dir)
  }
  }


