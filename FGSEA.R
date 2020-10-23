library(tidyverse)
library(fgsea)
library(msigdbr)

# Create a directory to save figures and tables

mainDir <- "../Figs"
subDir <- "FGSEA"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "FGSEA"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")



# MOdeling separately.
# Run FGSEA analysis 
trt <- c("HRZE","NTZ")[1]

pathway_list <- list()
for(trt in c("HRZE","NTZ")){
  # Limma results for genes
  res_trt <- read.csv(paste0("../Tables/Limma_rna_drug/",trt,"_Post_Vs_Pre_50.csv")
                      ,row.names = 1)
  res_TRT <- res_trt
  res_TRT$genes <- rownames(res_TRT)
  res_TRT$SYMBOL <- gsub(".*_","",rownames(res_trt))  
  
  library(tidyverse)
  res <- res_TRT
  
  rownames(res) <- gsub("\\..*","",rownames(res))
  res$row = rownames(res)
  
  res2 <- res %>% 
    dplyr::select(SYMBOL, t) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(t))
  res2
  
  ranks <- deframe(res2)
  head(ranks, 20)
  
  db_sets <-  c("H","C2","C5","C6","C7")
  
  res_sum_trt_list <- list()
  
  p_set <- db_sets[1]
  # Pathway set 
  m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == p_set)
  m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene,gene_symbol) %>% as.data.frame()
  
  # Select kegg,
  #m_t2g <- m_t2g[grep("KEGG",m_t2g$gs_name),]
  
  m_p_set <- m_t2g[,c("entrez_gene","gs_name")]
  
  #m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
  pathway_set <- unstack(m_p_set)
  pathways.hallmark <- unstack(m_t2g[,c("gene_symbol","gs_name")])
  
  
  
  # Look at them all if you want (uncomment)
  pathways.hallmark
  
  # Show the first few pathways, and within those, show only the first few genes. 
  pathways.hallmark %>% 
    head() %>% 
    lapply(head)
  
  set.seed(3000)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=20000,
                    minSize = 5, maxSize = 500)
  
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy$Group <- trt
  
  
  pathway_list[[trt]] <- fgseaResTidy
}  
# Draw an arrow shape instead of points

sig_path <- do.call("rbind",pathway_list)
# Import Pathways annotation
p_ann <-  read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv",
                   stringsAsFactors = F)
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
p_ann <- p_ann[p_ann$Hallmark.Name %in% as.character(sig_path$pathway),]
p_ann <- p_ann[match(sig_path$pathway,p_ann$Hallmark.Name),]
sig_path$Category <-  p_ann$Process.Category

# Save csv file for Fig 4 A
fig_dt <-  sig_path

fig_dt$leadingEdge <- unlist(lapply(fig_dt$leadingEdge,function(x) { paste(x,collapse = ",") }))


write.csv(fig_dt,paste0(tab_folder,"/Hallmark_pathways_FGSEA_HRZE_NTZ.csv"))



sig_path <- sig_path[sig_path$padj< 0.01,]
sig_path$Regulation <-  ifelse(sig_path$NES>0,"Up","Down")
sig_path$pathway <- gsub("HALLMARK_","",sig_path$pathway)

sig_path <- sig_path[order(sig_path$Category),]
sig_path$pathway <- factor(sig_path$pathway)

sig_path$leadingEdge <- unlist(lapply(sig_path$leadingEdge,function(x) { paste(x,collapse = ",") }))


write.csv(sig_path,paste0(tab_folder,"/Fig_4A.csv"))




library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
#col_dir <-  c("red","white","blue")
shape_vec <-  c("Up" = 24,"Down" = 25)
col_dir <- colfunc(7) 
p_hall <- ggplot(sig_path,aes(x = Group, y = pathway,fill = NES)) +
  geom_point(aes(shape = Regulation,size = -log10(padj)), color = "black", stroke = 1)+
  scale_shape_manual(values = shape_vec)+
  scale_fill_gradientn(limit = c(-4,4),colors = col_dir) +
  facet_wrap(~Category, strip.position = "top", scales = "free_y",ncol = 1)+
  theme_bw()+
  theme(axis.text.y=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        strip.placement = "outside")+
  ggtitle("Sig Hallmark Pathways across groups (padj < 0.01) ")+
  guides(size = guide_legend(override.aes = list(shape=24)))+
  ylab("Pathway")
#p_hall  

pdf(paste0(fig_folder,"/Hallmark_pathways_FGSEA_HRZE_NTZ.pdf"),width=8,height=14)
print(p_hall)
dev.off()



# Long Obs Analysis


# MOdeling separately.
# Run FGSEA analysis 
#trt <- c("TimeD14","TimeD56","TimeD56_TimeD14")[1]

limma_long_obs_dir <- "../Tables/Long_Obs_analysis/Limma_rnaseq/"

pathway_list <- list()

trt <-  c("Time_D14_Vs_D0","Time_D56_Vs_D0","TimeD56_Vs_D14")[3]

for(trt in c("Time_D14_Vs_D0","Time_D56_Vs_D0","TimeD56_Vs_D14")){
  # Limma results for genes
  res_trt <- read.csv(paste0(limma_long_obs_dir,trt,"_100.csv"),
                      row.names = 1)
  res_TRT <- res_trt
  res_TRT$genes <- rownames(res_TRT)
  res_TRT$SYMBOL <- gsub(".*_","",rownames(res_trt))  
  
  library(tidyverse)
  res <- res_TRT
  
  rownames(res) <- gsub("\\..*","",rownames(res))
  res$row = rownames(res)
  
  res2 <- res %>% 
    dplyr::select(SYMBOL, t) %>% 
    na.omit() %>% 
    distinct() %>% 
    group_by(SYMBOL) %>% 
    summarize(stat=mean(t))
  res2
  
  library(tidyverse)
  library(fgsea)
  ranks <- deframe(res2)
  head(ranks, 20)
  library(msigdbr)
  db_sets <-  c("H","C2","C5","C6","C7")
  
  res_sum_trt_list <- list()
  
  p_set <- db_sets[1]
  # Pathway set 
  m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == p_set)
  m_t2g = m_df %>% dplyr::select(gs_name, entrez_gene,gene_symbol) %>% as.data.frame()
  
  # Select kegg,
  #m_t2g <- m_t2g[grep("KEGG",m_t2g$gs_name),]
  
  m_p_set <- m_t2g[,c("entrez_gene","gs_name")]
  
  #m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
  pathway_set <- unstack(m_p_set)
  pathways.hallmark <- unstack(m_t2g[,c("gene_symbol","gs_name")])
  
  
  
  # Look at them all if you want (uncomment)
  pathways.hallmark
  
  # Show the first few pathways, and within those, show only the first few genes. 
  pathways.hallmark %>% 
    head() %>% 
    lapply(head)
  
  set.seed(3000)
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=20000,
                    minSize = 5, maxSize = 500)
  
  
  fgseaResTidy <- fgseaRes %>%
    as_tibble() %>%
    arrange(desc(NES))
  
  fgseaResTidy$Group <- trt
  pathway_list[[trt]] <- fgseaResTidy
  
}  
# Draw an arrow shape instead of points

sig_path <- do.call("rbind",pathway_list)
# Import Pathways annotation
p_ann <-  read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv",
                   stringsAsFactors = F)
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
p_ann <- p_ann[p_ann$Hallmark.Name %in% as.character(sig_path$pathway),]
p_ann <- p_ann[match(sig_path$pathway,p_ann$Hallmark.Name),]
sig_path$Category <-  p_ann$Process.Category

# Save csv file for Fig 5 D
fig_dt <-  sig_path

fig_dt$leadingEdge <- unlist(lapply(fig_dt$leadingEdge,function(x) { paste(x,collapse = ",") }))


write.csv(fig_dt,paste0(tab_folder,"/Hallmark_pathways_FGSEA_Long_Obs.csv"))



sig_path <- sig_path[sig_path$padj< 0.01,]
sig_path$Regulation <-  ifelse(sig_path$NES>0,"Up","Down")
sig_path$pathway <- gsub("HALLMARK_","",sig_path$pathway)

sig_path <- sig_path[order(sig_path$Category),]
sig_path$pathway <- factor(sig_path$pathway)

sig_path$Group <-  factor(sig_path$Group,levels = c("Time_D14_Vs_D0","TimeD56_Vs_D14","Time_D56_Vs_D0"))


sig_path$leadingEdge <- unlist(lapply(sig_path$leadingEdge,function(x) { paste(x,collapse = ",") }))


write.csv(sig_path,paste0(tab_folder,"/Fig_5D.csv"))



library(RColorBrewer)
colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))
#col_dir <-  c("red","white","blue")
shape_vec <-  c("Up" = 24,"Down" = 25)
col_dir <- colfunc(7) 
p_hall <- ggplot(sig_path,aes(x = Group, y = pathway,fill = NES)) +
  geom_point(aes(shape = Regulation,size = -log10(padj)), color = "black", stroke = 1)+
  scale_shape_manual(values = shape_vec)+
  scale_fill_gradientn(limit = c(-4,4),colors = col_dir) +
  facet_wrap(~Category, strip.position = "top", scales = "free_y",ncol = 1)+
  theme_bw()+
  theme(axis.text.y=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(size = 10,face = "bold"))+
  theme(axis.text.x=element_text(angle = -90, hjust = 0))+
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank(),
        strip.placement = "outside")+
  ggtitle("Sig Hallmark Pathways across groups (padj < 0.01) ")+
  guides(size = guide_legend(override.aes = list(shape=24)))
#p_hall  


pdf(paste0(fig_folder,"/Hallmark_pathways_FGSEA_Long_Obs.pdf"),width=8,height=14)
print(p_hall)
dev.off()
