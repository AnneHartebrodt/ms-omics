setwd('/work/AnneSophiaHartebrodt#6698/ms-project/')
require(rhdf5)
require(rhdf5filters)
require(Rhdf5lib)
renv::activate('/work/AnneSophiaHartebrodt#6698/ms-project/r-single-cell-new')
renv::repair()
set.seed(11)


require(Seurat)
library(tidyverse)
require(data.table)
require(ggplot2)
require(ggrepel)
require(colorspace)
#library(circlize)
library(wesanderson)
require(RColorBrewer)
require(harmony)
#library(UCell)
library(dplyr)
library(cowplot)
library(patchwork)
require(ggplot2)
library(MuSiC)
library(CARD)
require(dplyr)
require(liana)
require(hdWGCNA)

out.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna-spatial/cellinteraction')
dir.create(out.dir, showWarnings = F, recursive = T)

de.dir<- file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/plots')

base.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-spatial/')

filtering.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/filtering')
dir.create(filtering.dir)

liana.dir<-file.path('/work/AnneSophiaHartebrodt#6698/ms-project/submission/results-final/ms-rna/liana')

adata<- readRDS(file.path(filtering.dir, "annotated.rds"))
adata$batch<-sapply(adata$batch, function(x) gsub('smpD', 'MS9_NAWM', x))

liana_objects<-readRDS(file.path(liana.dir, 'all.Rds'))

lian<-list()
for(b in unique(adata$batch)){
  liana_test_C4_WM <- liana_objects[[b]] %>%
    liana_aggregate()
  liana_test_C4_WM$batch<-b
  liana_test_C4_WM$lesion<-stringr::str_split(b, '_')[[1]][2]
  lian[[b]]<-as.data.table(liana_test_C4_WM)
  
}
lian<-rbindlist(lian)
lian_sub<-lian %>% group_by(lesion, source, target, ligand.complex, receptor.complex) %>% 
  summarise(logFC= mean(sca.LRscore), rank=min(aggregate_rank), specificity = mean(natmi.edge_specificity, ))%>% 
  group_by(source, target) %>% slice_min(rank, n=10) 

lian_sub<-as.data.table(lian_sub)

get_interactions <- function(seurat_obj ,x,y) {
  out <- tryCatch(
    {
      expr <- FetchData(object = seurat_obj, vars = c(x, y))
      so<-seurat_obj[, which(x = expr[,1] > 1 & expr[,2]>1)]
      cells<-rownames(so@meta.data)
    },
    error=function(cond) {
      # Choose a return value in case of error
      return(NA)
    }
  )    
  return(out)
}



ligand<-c("HLA-A", "C3", "APP", "APOE", "APOE")
receptor<-c("APLP2", "LRP1", "APLP2",  "RTREM2", "LRP1")
ligands<-data.table(ligand, receptor)

interactions_list<-list()

lian_unique<-unique(lian[aggregate_rank<0.01][,.(ligand.complex, receptor.complex)])

for(ds in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  seurat_obj<-readRDS(file.path(base.dir, ds, paste0(ds, '_processed_.rds')))
  for(i in 1:nrow(lian_unique)){
    inter<-get_interactions(seurat_obj, lian_unique$ligand.complex[i], lian_unique$receptor.complex[i])
    if(length(inter)>0 &&!is.na(inter)){
    interactions<-data.table(inter, lian_unique$ligand.complex[i], lian_unique$receptor.complex[i],  ds)
    interactions_list[[paste0(ds, i)]]<-interactions
    }
  }
}

spot_pred<-list()
for(ds in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  card_obj<-readRDS(file.path(base.dir, ds, paste0(ds, '_card_object.rds')))
  rown<-rownames(card_obj@Proportion_CARD)
  prop<-as.data.table(card_obj@Proportion_CARD)
  prop<-round(prop, digits = 2)
  prop$rownames<-rown
  prop$dataset<-ds
  spot_pred[[ds]]<-prop
}
prop<-rbindlist(spot_pred)

interactions<-rbindlist(interactions_list, fill = T)
interactions<-interactions[!is.na(inter)]
interactions.summary<- interactions[, .N, by = c('V2', 'V3','ds')][order(N, decreasing = T)] %>%
  pivot_wider(names_from = 'ds', values_from = 'N')
interactions.summary<-as.data.table(interactions.summary)

fwrite(interactions, file=file.path(out.dir, 'interaction_spots.tsv'), sep='\t')
fwrite(interactions.summary, file=file.path(out.dir, 'interaction_summary.tsv'), sep='\t')

#interactions<-fread(file.path(out.dir, 'interaction_spots.tsv'))
#interactions.summary<-fread(file.path(out.dir, 'interaction_summary.tsv'))

interactions<-merge(interactions, prop, by.x = c('inter', 'ds'), by.y = c('rownames', 'dataset'))

interactions.weighted<-merge(interactions, lian[aggregate_rank<0.01], by.x = c('V2', 'V3'), by.y= c('ligand.complex', 'receptor.complex'), allow.cartesian = T)
interactions.summary.weighted<- interactions[, .N, by = c('V2', 'V3')][order(N, decreasing = T)]

interactions.summary.global<- interactions[, .N, by = c('V2', 'V3')][order(N, decreasing = T)]


interactions.summary<- unique(rbind(interactions.summary.global[1:20, 1:2], interactions.summary[1:20, 1:2]))

plot_list<-list()
for(dss in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  seurat_obj<-readRDS(file.path(base.dir, dss, paste0(dss, '_processed_.rds')))
  for (i in 1:20){
    spots<-interactions[(V2 == interactions.summary$V2[i] & V3 == interactions.summary$V3[i] & ds == dss)]
    dp<-SpatialDimPlot(seurat_obj, cells.highlight = spots$inter , cols.highlight = c('yellow', 'grey'), stroke = 0)+
          guides(fill='none')
    plot_list[[dss]][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]]<-dp
  }
}


dir.create(file.path(out.dir, 'interaction_plots'))
for (i in 1:20){
    pg<-plot_grid(plot_list[['CA_2A']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]], 
     plot_list[['CA_2B']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]],
     plot_list[['CA_2C']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]],
     plot_list[['RL_1A']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]],
     plot_list[['RL_1B']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]],
     plot_list[['RL_1C']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]],
     plot_list[['RL_1D']][[paste0(interactions.summary$V2[i], interactions.summary$V3[i])]],
     nrow = 1)
    
    title <- ggdraw() + 
      draw_label(
       paste('Interactions between', interactions.summary$V2[i], 'and', interactions.summary$V3[i]),
        fontface = 'bold',
        x = 0,
        hjust = 0
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 0.1)
      )
    pg<-plot_grid(
      title, pg,
      ncol = 1,
      # rel_heights values control vertical title margins
      rel_heights = c(0.1, 1)
    )
    ggsave(pg, file = file.path(out.dir, 'interaction_plots', paste0(interactions.summary$V2[i], '_', interactions.summary$V3[i], '.pdf')), height = 5, width = 30, units = 'cm')
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i1<-c('TF')
i2<-c('LRP2')

plot_list<-list()
for(dss in c('CA_2A', 'CA_2B', 'CA_2C', 'RL_1A','RL_1C', 'RL_1B', 'RL_1D')){
  seurat_obj<-readRDS(file.path(base.dir, dss, paste0(dss, '_processed_.rds')))
  for (i in 1:length(i1)){
    spots<-interactions[(V2 == i1[i] & V3 == i2[i] & ds == dss)]
    dp<-SpatialDimPlot(seurat_obj, cells.highlight = spots$inter , cols.highlight = c('yellow', 'grey'), stroke = 0)+
      guides(fill='none')
    plot_list[[dss]][[paste0(i1[i], i2[i])]]<-dp
  }
}


dir.create(file.path(out.dir, 'interaction_plots'))
for (i in 1:length(i1)){
  pg<-plot_grid(plot_list[['CA_2A']][[paste0(i1[i], i2[i])]], 
                plot_list[['CA_2B']][[paste0(i1[i], i2[i])]],
                plot_list[['CA_2C']][[paste0(i1[i], i2[i])]],
                plot_list[['RL_1A']][[paste0(i1[i], i2[i])]],
                plot_list[['RL_1B']][[paste0(i1[i], i2[i])]],
                plot_list[['RL_1C']][[paste0(i1[i], i2[i])]],
                plot_list[['RL_1D']][[paste0(i1[i], i2[i])]],
                nrow = 1)
  
  title <- ggdraw() + 
    draw_label(
      paste('Interactions between', i1[i], 'and', i2[i]),
      fontface = 'bold',
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 0.1)
    )
  pg<-plot_grid(
    title, pg,
    ncol = 1,
    # rel_heights values control vertical title margins
    rel_heights = c(0.1, 1)
  )
  ggsave(pg, file = file.path(out.dir, 'interaction_plots', paste0(i1[i], '_', i2[i], '.pdf')), height = 5, width = 30, units = 'cm')
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lian$source<-recode(lian$source, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPC', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')
lian$target<-recode(lian$target, 'Oligodendrocytes'='OL', 'Astrocytes'='AS', 'Immune cells'='IM', 'OPCs'='OPC', 'B-cells'='BC', 'Neurons'='NC','Pericytes'='PC', 'Endothelial cells'='EC')

gp<-ggplot(lian[aggregate_rank<0.01  &(ligand.complex %in% interactions.summary$V2[1:20] & receptor.complex %in% interactions.summary$V3[1:20])], 
           aes(x = lesion, y = paste0(ligand.complex, '->', receptor.complex), size=natmi.edge_specificity, color=sca.LRscore))+geom_point()+
  scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        title = element_text(size = 18, vjust = -0.5),
        axis.title.y.right = element_text())+
  guides(color=guide_colorbar('Log(FC)'))+
  scale_size_continuous(range = c(3,8)) +
  ggh4x::facet_nested("Source" + source ~ "Target" + target, scales = 'free',  space="free")+
  ylab("")+xlab("")
  
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('top20_in_spatial.pdf'))), width = 20, height = 10, units = 'in')


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# CELLTYPE INTERACTIONS PLOTS
lian_sub.1<-lian %>% group_by(lesion, source, target, ligand.complex, receptor.complex) %>% 
  summarise(logFC= mean(sca.LRscore), rank=min(aggregate_rank), specificity = mean(natmi.edge_specificity, ))%>% 
  group_by(source, target)

lian_sub.1<-as.data.table(lian_sub.1)
lian_sub.1<-lian_sub.1[specificity>0.4]

sum<-interactions[, .N,  by = c('V2', 'V3')][order(N, decreasing = T)]

dir.create(file.path(out.dir, 'liana_filtered'))
cct<-'OL'
inter<-sum[N>50]
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0(cct, ' >> other celltypes'))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(paste0(cct, '_to_others.pdf')))), width = 10.25, height = 5.5, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'vertical')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 9.25, units = 'in')





cct<-'AS'
gp<-ggplot(lian_sub.1[source==cct & specificity >0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('Astrocytes >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 10.5, height = 10, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 12, units = 'in')



cct<-'NC'
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('Neurons >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 12.5, height = 12, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 7, units = 'in')



cct<-'EC'
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('Endothelial cells >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 10, height = 8, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 10, units = 'in')

cct<-'PC'
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('Pericytes >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 12.5, height = 22, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 12, units = 'in')


cct<-'BC'
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('B-cells >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 12.5, height = 10, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 8, units = 'in')


cct<-'IM'
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('Immune cells >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 12.5, height = 12, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 7, units = 'in')



cct<-'OPC'
gp<-ggplot(lian_sub.1[source==cct & specificity>0.4 & (ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3)], aes(x = target, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c('ligand.complex', 'lesion'),scales = 'free',  space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y = element_blank(),
        title = element_text(size = 18, vjust = -0.5))+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle('OPCs >> other celltypes')
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0(cct,'_to_others.pdf'))), width = 12.5, height = 8, units = 'in')

gp<-ggplot(lian_sub.1[target==cct & specificity>0.4 &(ligand.complex %in% inter$V2 | receptor.complex %in% inter$V3) ], 
           aes(x = source, y = paste0(ligand.complex, '->', receptor.complex), size=specificity, color=logFC))+geom_point()+
  facet_grid(c( 'receptor.complex', 'lesion'),scales = 'free', space="free")+scale_color_viridis_b()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 18), 
        axis.text.y = element_text(size = 18),
        strip.background = element_blank(), 
        strip.text = element_text(size = 18),
        strip.text.y =  element_blank(),
        title = element_text(size = 18),
        legend.box = 'horizontal')+
  guides(color=guide_colorbar('Log(FC)'))+
  ylab("")+xlab('')+
  scale_size_continuous(range = c(3,8)) +ylab("")+xlab('')+
  ggtitle(paste0('Other celltypes >>', cct))
gp
ggsave(gp, file = file.path(file.path(out.dir, 'liana_filtered', paste0('others_to_',cct,'.pdf'))), width = 12.5, height = 10, units = 'in')





