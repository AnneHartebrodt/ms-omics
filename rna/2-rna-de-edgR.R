setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5filters)
require(rhdf5)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell')
renv::repair()
set.seed(11)

require(rlang)
require(Seurat)
require(data.table)
require(SingleCellExperiment)
require(ggplot2)
require(ggrepel)
require(dplyr)
require(tidyr)
require(devtools)
require(colorspace)
require(ComplexHeatmap)
library(circlize)
library(wesanderson)
require(RColorBrewer)
require(magick)
require(ggvenn)
require(nVennR)
library(stringi)
library(Cairo)
library(rsvg)
library(grImport2)
library(WGCNA)
library(hdWGCNA)
library(scran)
library(scater)
require(edgeR)
require(warrenlabSC)


res.dir <-
  file.path(
    '/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/differential-edgeR'
  )
dir.create(res.dir, recursive = T)

filtering.dir <-
  file.path(
    '/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/',
    'filtering'
  )

#1. Read data object.
adata <- readRDS(file.path(filtering.dir, "annotated.rds"))
keeplist <- list()

#2. Perform differential analysis for all celltypes separately
for (ct in c('Astrocytes', 'Oligodendrocytes', 'Immune cells', 'OPCs')) {
  print(ct)
  
  cell.dir <- file.path(res.dir, ct)
  dir.create(cell.dir)
  
  list_of_de <- list()
  for (lt in c('NAWM', 'RL', 'CA', 'AL')) {
    astro <-
      subset(adata, subset = celltypes_annotated == ct &
               lesion_type %in% c('WM', lt))
    
    Idents(astro) <- 'lesion_type'
    astro$batch <- as.factor(astro$batch)
    sse <- as.SingleCellExperiment(astro)
    cdr <- calc.cdr(astro@assays$RNA@counts)
    sse$cellular_detection_rate <- as.numeric(cdr)
    
    # Perform pseudobulk aggregation
    summed <- aggregateAcrossCells(sse,
                                   id = colData(sse)[, c("celltypes_annotated", "batch")])
    
    
    current <- summed
    #current<-sse[,label==summed$celltypes_annotated]
    # Creating up a DGEList object for use in edgeR:
    
    
    y <- DGEList(counts(current), samples = colData(current))
    y
    
    # remove groups with fewer than 10 cells
    discarded <- current$ncells < 10
    y <- y[, !discarded]
    summary(discarded)
    
    # keep genes present in at least 50% of the samples per group
    keep <- filterByExpr(y, group = current$batch, min.prop = 0.5)
    y <- y[keep, , keep.lib.sizes = FALSE]
    summary(keep)
    keeplist[[paste0(ct, lt)]] <- keep
    
    
    y <- calcNormFactors(y)
    y$samples
    
    par(mfrow = c(2, 3))
    for (i in seq_len(ncol(y))) {
      plotMD(y, column = i)
    }
    
    
    #plotMDS(cpm(y, log=TRUE),
    #        col=ifelse(y$samples$batch, "red", "blue"))
    
    
    design <- model.matrix( ~ 1 + as.factor(sex) + as.factor(ms), y$samples)
    design
    
    y <- estimateDisp(y, design)
    
    summary(y$trended.dispersion)
    
    #plotBCV(y)
    
    fit <- glmFit(y, design)
    #fit <- glmQLFit(y, design,robust = TRUE)
    summary(fit$var.prior)
    
    #summary(fit$df.prior)
    #plotQLDisp(fit)
    
    res <- glmLRT(fit, coef = ncol(design))
    #res <- glmQLFTest(fit, coef=ncol(design))
    
    s <- summary(decideTests(res))
    
    ta <- topTags(res, n = 1000)
    ta <- ta@.Data[[1]]
    na <- rownames(ta)
    ta$gene <- na
    ta$lesion <- lt
    ta$celltype <- ct
    ta <- as.data.table(ta)
    ta <- ta[FDR < 0.05]
    
    list_of_de[[lt]] <- ta
  }
  
  
  astro <- subset(adata, subset = celltypes_annotated == ct)
  
  Idents(astro) <- 'lesion_type'
  
  degenes <- rbindlist(list_of_de)
  degenes <- as.data.table(degenes)
  degenes$up <- ifelse(degenes$logFC > 0, 'upregulated', 'downregulated')
  degenes[, .N, by = c('gene', 'up')][order(N, decreasing = T)]
  fwrite(degenes, file = file.path(cell.dir, paste0('DE_genes_', ct, '.tsv')), sep =
           '\t')
  
  multiples <-
    degenes[, .N, by = c('gene', 'up')][order(N, decreasing = T)][N >= 2]
  singletons <-
    degenes[, .N, by = c('gene', 'up')][order(N, decreasing = T)][N == 1]
  
  cell.dir.pl <- file.path(cell.dir, 'plots')
  dir.create(cell.dir.pl)
  for (gg in unique(degenes$gene)) {
    vln <-
      VlnPlot(astro, features = gg, group.by = 'batch') + theme(axis.text.x = element_text(angle = 90),
                                                                legend.position = 'none') |
      VlnPlot(astro, features = gg, group.by = 'lesion_type') + theme(axis.text.x = element_text(angle = 90))
    ggsave(
      vln,
      file = file.path(cell.dir.pl, paste0(gg, '.pdf')) ,
      width = 20,
      height = 10,
      units = 'cm'
    )
  }
  
  if (length(multiples$gene) > 0) {
    exp <- FetchData(astro, multiples$gene)
    exp <- exp > 0
    exp <- as.data.table(exp)
    exp$lesion_type <- astro$lesion_type
    
    sum <- exp[, lapply(.SD, mean), by = lesion_type]
    sum <- as.data.table(sum)
    sum <- sum[, -c('lesion_type')]
    sumsum <- sum[, lapply(.SD, sum)]
    sumsum <- sumsum / nrow(sum)
    deg <- multiples[which(sumsum > 0.05)]
    
    if (length(deg[up == 'upregulated']$gene) > 0) {
      multi_up <-
        DotPlot(astro, features = deg[up == 'upregulated']$gene, scale = F) + theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) +
        scale_color_continuous_divergingx('RdBu') + ggtitle('Upregulated genes in 2 or more lesion types')
      multi_up
      ggsave(
        multi_up,
        file = file.path(
          cell.dir,
          paste0(ct, '_upregulated_2_or_more_lesion_types.pdf')
        ) ,
        width = 10 + (max(0, (
          length(deg[up == 'upregulated']$gene) - 5
        ))),
        height = 8,
        units = 'cm'
      )
    }
    
    if (length(deg[up == 'downregulated']$gene) > 0) {
      multi_down <-
        DotPlot(astro, features = deg[up == 'downregulated']$gene, scale = F) +
        theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) +
        scale_color_continuous_divergingx('RdBu') + ggtitle('Downregulated genes in 2 or more lesion types')
      multi_down
      ggsave(
        multi_down,
        file = file.path(
          cell.dir,
          paste0(ct, '_downregulated_2_or_more_lesion_types.pdf')
        ) ,
        width = 10 + (max(0, (
          length(deg[up == 'downregulated']$gene) - 5
        ))),
        height = 8,
        units = 'cm'
      )
    }
  }
  
  for (lt in c('NAWM', 'RL', 'CA', 'AL')) {
    deg <- degenes[lesion == lt & gene %in% singletons$gene]
    
    
    
    if (nrow(deg) > 0) {
      exp <- FetchData(astro, deg$gene)
      exp <- exp > 0
      exp <- as.data.table(exp)
      exp$lesion_type <- astro$lesion_type
      
      sum <- exp[, lapply(.SD, mean), by = lesion_type]
      sum <- as.data.table(sum)
      sum <- sum[, -c('lesion_type')]
      sumsum <- sum[, lapply(.SD, sum)]
      sumsum <- sumsum / nrow(sum)
      deg <- deg[which(sumsum > 0.1)]
      if (nrow(deg) > 0) {
        if (nrow(deg) < 60) {
          single_down <- DotPlot(astro, features = deg$gene, scale = F) +
            facet_wrap( ~ deg$up, scales = 'free') +
            scale_color_continuous_divergingx('RdBu') +
            ggtitle(paste0('Genes only DE in ', lt, ' vs. WM (all lesions shown)')) +
            theme(axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ))
          single_down
          ggsave(
            single_down,
            file = file.path(cell.dir, paste0(
              ct, '_up_or_down_in_', lt, '.pdf'
            )),
            width = min(40, 20 + (max(
              0, (nrow(deg) - 5)
            ))),
            height = 8,
            units = 'cm'
          )
        } else{
          single_down <- DotPlot(astro, features = deg$gene, scale = F) +
            facet_wrap( ~ deg$up, scales = 'free', nrow = 2) +
            scale_color_continuous_divergingx('RdBu') +
            ggtitle(paste0('Genes only DE in ', lt, ' vs. WM (all lesions shown)')) +
            theme(axis.text.x = element_text(
              angle = 90,
              vjust = 0.5,
              hjust = 1
            ))
          single_down
          ggsave(
            single_down,
            file = file.path(cell.dir, paste0(
              ct, '_up_or_down_in_', lt, '.pdf'
            )),
            width = min(40, 20 + (max(
              0, (nrow(deg) / 2 - 5)
            ))),
            height = 20,
            units = 'cm'
          )
          
        }
      }
    }
  }
}
