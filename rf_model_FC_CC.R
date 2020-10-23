#Load libraries
library(DESeq2)
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(ggthemes)
library(msigdbr)
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(GSVA)




# Create a directory to save figures and tables

mainDir <- "../Figs"
subDir <- "RF_FC_CC"
#subDir <- "Figs_Clostridia_kegg"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")


# Import microbiome data

phy_fc_cc <- readRDS("../data/Control/phy_mic_fc_cc.rds")


#RNASeq deseq normalized
phy_vst_rna <-  readRDS("../data/Control/phy_rna_vst_fc_cc.rds")


phy_gsva <- phy_vst_rna

counts <- otu_table(phy_gsva) %>% as.data.frame %>% as.matrix() #or as.matrix() #TPM data here with gene rownames and sample colnames

#metadata from RNAseq phyloseq object
metadata <- data.frame(sample_data(phy_gsva))

#gene lists
msigdbr_show_species()
m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == "H" )

m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_t2g <- m_t2g[,c("gene_symbol","gs_name")]
#m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
hall_set <- unstack(m_t2g)

rownames(counts) <- gsub(".*_","",rownames(counts))

set.seed(1057)
gsvaRes_ssgsea <- gsva(counts, hall_set,
                       min.sz=5, max.sz=500,
                       method = "gsva",
                       kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1)

met  <- metadata
mat  <- gsvaRes_ssgsea 

met <- met[match(colnames(gsvaRes_ssgsea),met$rna_seq_ID),]

library(ComplexHeatmap)
library(RColorBrewer)
status_col <-  c("family_contact" = "red","community_control" = "blue")

library(gtools)
ha_column = HeatmapAnnotation(Status =  met$TB_status,
                              Batch = met$batchID,
                              Sex = met$sex,
                              col=list(Status = status_col))

split_cols<-  met$TB_status
split_cols <- factor(split_cols, levels= c("community_control","family_contact","Pre HRZE","Post HRZE","Pre NTZ","Post NTZ"))

colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))
cols <- colfunc(8)
rownames(mat) <-  gsub("HALLMARK_","",rownames(mat))


ht1 = Heatmap(mat, name = "ES", column_title = NA, 
              top_annotation = ha_column,
              #clustering_distance_rows = "euclidean",
              col = cols,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 10,fontface = "bold"), 
              row_title_rot = 0,
              column_split = split_cols,
              width=2, cluster_columns = T, 
              row_names_gp = gpar(fontsize = 9),
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              show_column_names = F, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
# pdf("GSVA_deseq_norm_vst.pdf",width = 18,height = 8)
# draw(ht1)
# dev.off()


# Microbiome:
phy_sel_mic <-  prune_taxa(taxa_sums(phy_fc_cc)>0,phy_fc_cc)

phy <- phy_sel_mic
colnames(tax_table(phy))[3] <- "Kingdom"


# Prevalence 10 % 
perc_samples <-  0.1
# # 
phy_fil = filter_taxa(phy, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)
phy_fil <- prune_taxa(taxa_sums(phy_fil)> 0, phy_fil)

#phy <- yingtools2::phy.collapse(phy,c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))
############# ############# #############

#view the rownames as Phylum, Species, rather than OTU_number
t <- get.tax(phy_fil) %>% mutate(PhySpec=paste(Species,otu))
taxa_names(phy_fil) <- t$PhySpec
taxa_names(phy_fil) <- gsub("\\[","",taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub("\\]","",taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub(';.*','\\1', taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub("ASV","(ASV",taxa_names(phy_fil))
taxa_names(phy_fil) <- gsub("$",")",taxa_names(phy_fil)) #edit this between old/new methods
taxa_names(phy_fil) <- gsub(" ","_",taxa_names(phy_fil)) #edit this between old/new methods
taxa_names(phy_fil)

# # Some manual changes in taxa table
# 
tax_table(phy_fil) <-  tax_table(phy_fil)[,3:9]


phy_sel_mic <- phy_fil


# # Deseq normalized counts data
dig <- phyloseq_to_deseq2(phy_sel_mic, ~  type) #replace this with any sample variable(s)

dig <- estimateSizeFactors(dig,"poscounts")
# Takes a while to run
dig <- estimateDispersions(dig)
dig <- DESeq(dig,fitType= "local")

plotDispEsts(dig)
# VST phyloseq for Heatmap later

phy_vst_mic <- phy_sel_mic
vst_dt <- getVarianceStabilizedData(dig)
otu_table(phy_vst_mic) <-  otu_table(vst_dt, taxa_are_rows = T)

# Lets regression each pathway against the microbes
# For RNA data
pathways_dt <-   data.frame(t(gsvaRes_ssgsea)) 

dt_SV <-  data.frame(t(otu_table(phy_vst_mic)))
dt_SV <-  dt_SV[match(rownames(pathways_dt),rownames(dt_SV)),]

rownames(dt_SV)
rownames(pathways_dt)



###################### Microbiome and Pathway#########################
pth_dt <- pathways_dt

# Store all the RF plots inside the folder below

mainDir <- "../Figs/RF_FC_CC"
subDir <- "H_Pathway"

fig_folder <-  paste0(mainDir,"/",subDir)
dir.create(fig_folder,showWarnings = TRUE,recursive = TRUE)


mic_dt <-  dt_SV
match(rownames(mic_dt),rownames(pth_dt))
names(pth_dt) 
names(mic_dt)



# Loop across the pathways
comp <-  names(pth_dt)[20]

imp_dt_list <- list()
for(comp in names(pth_dt)[1:2]){
  
  t_data <- cbind(mic_dt,Y_var = pth_dt[,comp])
  
  
  library(randomForest)
  # Overall model using all data
  set.seed(1057)
  rf_model <- randomForest(Y_var~.,data=t_data,
                           #importance = "impurity_corrected",
                           #importance = "permutation",
                           importance = T,
                           #probability = T,
                           ntree = 5000)
  
  library(vita)
  
  #t_data <- dt_orig
  train_data <-  t_data[,!colnames(t_data) %in% "Y_var"]
  pimp <- PIMP(train_data, t_data$Y_var,rf_model,parallel=TRUE, ncores=10 ,seed = 340)
  p_test <- PimpTest(pimp)
  #pimp_all <- data.frame( summary(p_test,pless = 0.05))
  pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
  #pimp_all$padj <- p.adjust(pimp_all$p.value,method = "BH")
  pimp_all$padj <- pimp_all$p.value
  imp_dt <- pimp_all[pimp_all$VarImp > 0 & pimp_all$padj <= 0.05,]     
  
  
  
  imp_dt$variable <- rownames(imp_dt)
  names(imp_dt)
  imp_dt <- imp_dt[order(imp_dt$VarImp,decreasing = T),]
  
  if(nrow(imp_dt)>0){
    
    sel_imp_dt <- imp_dt
    
    library(ggplot2)
    p_var <- ggplot(sel_imp_dt, aes(x=reorder(variable,VarImp), y=VarImp))+ 
      geom_bar(stat="identity", position="dodge")+ coord_flip()+
      ylab("Variable Importance")+
      xlab("")+theme_bw()
    
    pdf(paste0(fig_folder,"/Variable_imp_RF_",comp,".pdf"), width = 8, height = 10)
    print(p_var)
    dev.off()
    
    
    # ALE Plots 
    # ALE analysis
    library(ALEPlot)
    library(parallel)
    yhat <- function(X.model, newdata) as.numeric(randomForest:::predict.randomForest(X.model, newdata))
    
    ALE_summary_list <- list()
    imp_bugs <-  as.character(sel_imp_dt$variable)
    
    train_data <-  t_data[,!colnames(t_data) %in% "Y_var"]
    
    ale_func <-  function(rep){
      
      if(imp_bugs[len] %in% colnames(train_data)){
        idx <- which(colnames(train_data) == imp_bugs[len]) 
        ale_dt  <- data.frame(ALEPlot(X = train_data,X.model =  rf_model, 
                                      pred.fun = yhat, J=idx, K = 40, NA.plot =F))
        
        ale_dt$pred <-  imp_bugs[len]
        ale_dt$nrep <- rep
        res <-  ale_dt
        # rm(ale_dt)
      }
    }
    
    pp_bug_list <- list()
    for(len in seq_along(imp_bugs)){
      print(len)
      ale_dt <- ale_func(1)
      pp_bug_list[[len]] <- ale_dt
    }
    
    pp_dt_final <- do.call("rbind",pp_bug_list)
    
    pp_dt_final$f.values <-  as.numeric(pp_dt_final$f.values)
    pp_dt_final$x.values <-  as.numeric(pp_dt_final$x.values)
    pp_dt_final <-  na.omit(pp_dt_final)
    
    
    pp_dt_final$pred <- factor(pp_dt_final$pred,levels = imp_bugs)
    
    library(ggthemes)
    p_p <- ggplot(data = pp_dt_final ,aes(x = x.values, y = f.values)) +
      geom_line(aes(group = nrep), alpha = 0.6)+
      #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = trial))+
      #geom_smooth()+
      #geom_line(data = pp_dt_final_median ,aes(x = x, y = y_med, color = "red"))+
      theme_base()+
      facet_wrap(~pred,scales = "free")+
      xlab("Microbes")+
      theme(legend.position="none")
    #  print(p_p)
    pdf(paste0(fig_folder,'/RF_ALE_',comp,'.pdf'),height = 15, width = 15)
    print(p_p)
    dev.off()
    
    library(tidyr)
    library(dplyr)
    slope_dt <- pp_dt_final %>%
      group_by(pred) %>%
      do({
        mod <-  data.frame(summary(lm( f.values ~ x.values , data = .))$coefficients)[,c(1,4)]
        names(mod) <- c("val","p_val")
        mod$coeff <- rownames(mod)
        
        data.frame(Intercept = mod$val[1],
                   Slope = mod$val[2],
                   p_val_int =  mod$p_val[1],
                   p_val_slope =  mod$p_val[2] )
      }) %>% as.data.frame()
    
    head(slope_dt)
    
    
    
    imp_slope <-  sel_imp_dt
    imp_slope$Slope <-  slope_dt$Slope
    imp_slope$Relation <- imp_slope$Slope
    imp_slope$Relation <-  ifelse(imp_slope$Relation > 0, "Positive","Negative")
    imp_slope$Intercept <- slope_dt$Intercept
    imp_slope$p_val_int <- slope_dt$p_val_int
    imp_slope$p_val_slope <- slope_dt$p_val_slope
    
    imp_slope$Comp <- comp
    
    imp_dt_list[[comp]] <- imp_slope 
    
    #Predictors to visualize
    pred_vec <-  as.character(sel_imp_dt$variable)
    dt_sel <- t_data[,c(pred_vec,"Y_var")]
    dt_sel$sample_id <- rownames(dt_sel)
    head(dt_sel)
    #rownames(X_sel)
    
    library(reshape2)
    melt_dt <-  reshape2::melt(dt_sel, id.vars  = c("Y_var","sample_id"))
    
    
    head(melt_dt)
    
    library(tidyverse)
    
    xy_p <- ggplot(melt_dt, aes(y = Y_var, x = value)) +
      geom_point()+
      facet_wrap(~ variable, scales = "free") +
      theme_bw() +
      ylab("Abundance") +
      xlab("Induction")
    pdf(paste0(fig_folder,"/XY_plot_",comp,".pdf"), width = 10, height = 10)
    print(xy_p)
    dev.off()
  }
  
}

# Based on the importance , lets draw the heatmap which only demonstrates the
# importance of the species towards the Compounds


# Save the run
saveRDS(imp_dt_list,"../Figs/RF_FC_CC/RF_imp_FC_CC.rds")
imp_dt_list <- readRDS("../Figs/RF_FC_CC/RF_imp_FC_CC.rds")
comp_imp <- do.call("rbind",imp_dt_list)

#Select only top 20
comp_imp <- comp_imp %>%
  dplyr::group_by(Comp) %>%
  dplyr::arrange(desc(VarImp)) %>%
  dplyr::slice_max(n= 10,VarImp)%>% data.frame()



comp_imp$value <- 1
comp_imp$value[comp_imp$Slope < 0 ] <- -1 


comp_imp <- comp_imp[,c("variable","Comp","value")]
library(tidyr)

comp_imp_w <- comp_imp %>%
  pivot_wider(names_from = Comp, values_from = value,values_fill = list(value = 0))%>%
  data.frame()

rownames(comp_imp_w) <-  comp_imp_w$variable
comp_imp_w$variable <- NULL

mat <- as.matrix(comp_imp_w) 

#col_binary <- colorRampPalette(rev(brewer.pal(n = 3, name = "RdYlBu")))

col_binary <-  c("#91BFDB","white","#FC8D59")

# Import R2 information
rsq_dt <- read.csv("../Random_Forest/Control/Predictions_R2_Mic_H_Pathways.csv")
#rsq_dt <-  rsq_dt[rsq_dt$Grp %in% vec_phy[loop],]
rsq_dt <-  rsq_dt[rsq_dt$Comp %in% colnames(mat),]
rsq_dt <-  rsq_dt[rsq_dt$R2 >= 0.1,]

# Order by R2 
rsq_dt <- rsq_dt[order(rsq_dt$R2,decreasing = T),]


mat <-  mat[,colnames(mat) %in% rsq_dt$Comp]
mat <- mat[apply(mat, 1, function(x) !all(x==0)),]
mat <-  mat[,match(rsq_dt$Comp,colnames(mat))]


# Import Pathways annotation
p_ann <-  read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv",
                   stringsAsFactors = F)
p_ann$Hallmark.Name <-  paste0("HALLMARK_",p_ann$Hallmark.Name)
p_ann <- p_ann[p_ann$Hallmark.Name %in% as.character(rsq_dt$Comp),]
p_ann <- p_ann[match(rsq_dt$Comp,p_ann$Hallmark.Name),]
# # Match it with the colnames of the matrix
# rsq_dt <-  rsq_dt[match(colnames(mat),rsq_dt$Comp),]
# 
sp_cols <- factor(p_ann$Process.Category)

library(RColorBrewer)
library(circlize)
#sel_cols = colorRampPalette(c("white", "darkblue"))(2)
col_bar = colorRamp2(c(0,1), c("white", "darkblue"))

# col_ha = HeatmapAnnotation(Rsq =  anno_simple(rsq_dt$R2 
#                                        ) ,
#                            col = list(Rsq = col_bar),
#                            border = TRUE)
ha1 = HeatmapAnnotation(show_annotation_name = F,
                        Rsq = rsq_dt$R2, 
                        col = list(Rsq = col_bar),
                        Rq = anno_barplot(
                          rsq_dt$R2, 
                          bar_width = 1, 
                          gp = gpar(col = "white", fill = "black"), 
                          border = T,
                          ylim = c(0,1),
                          axis_param = list(at = c(0, 0.2, 0.4,0.6,0.8,1.0),
                                            labels = as.character(c(0, 0.2, 0.4,0.6,0.8,1.0))),
                          height = unit(1.5, "cm")
                        ),
                        border = T
)


# Add tax dt info
phy_ht <- phy_sel_mic
taxa_names(phy_ht) <- make.names(taxa_names(phy_ht))

phy_ht <-  prune_taxa(rownames(mat),phy_ht)
phy_melt <-  get.otu.melt(phy_ht)
tax_dt <-  unique(phy_melt[,c("otu","Order","Phylum")])
topN <- length(unique(tax_dt$Order))
match_tax_dt <- data.frame(tax_table(phy_ht))
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]
rownames(match_tax_dt)


phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Tenericutes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Cyanobacteria" = "#808080",
                "Euryarchaeota" = "#8c8c8c",
                "Spirochaetes" = "grey50",
                "Candidatus Melainabacteria" = "#999999")

dt_col_phylum <- data.frame(phylum_col)
dt_col_phylum$Phylum <-  rownames(dt_col_phylum)
shades_num <- match_tax_dt %>%
  group_by(Phylum) %>%
  mutate(unique_types = n_distinct(Order))%>% select(Phylum,unique_types) %>% unique %>% as.data.frame()

shades_num <- merge(shades_num, dt_col_phylum, by = "Phylum")
shades_num$phylum_col <-  as.character(shades_num$phylum_col )
tax_sel <-  unique(match_tax_dt[,c("Phylum","Order")])

#shades_num$phy_col <- as.vector(phylum_col[order(names(phylum_col))])
order_col <- mapply(function(x, y){ shades(color = x, ncolor = y, variation = 1)}, x= shades_num$phylum_col ,y = shades_num$unique_types,
                    SIMPLIFY = F)

order_col <- as.vector(unlist(order_col))
names(order_col) <- tax_sel$Order

mycol <- order_col

match_tax_dt <- match_tax_dt[match(rownames(mat),rownames(match_tax_dt)),]




ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))


# Replace dots with space
rownames(mat)  <- gsub("\\."," ",rownames(mat))
rownames(mat)  <- gsub("X","",rownames(mat))

colnames(mat) <- gsub("HALLMARK_","",colnames(mat))
colnames(mat) <- gsub("_"," ",colnames(mat))

#col_ha = rowAnnotation(Rsq = rsq_dt$R2)

#mat <- t(mat)
cn = colnames(mat)
bottom_annotation = HeatmapAnnotation(
  text = anno_text(cn, rot = 90, offset = unit(1, "npc"), just = "right"),
  annotation_height = max_text_width(cn)
)


rownames(mat) <- gsub("_"," ",rownames(mat))  
ht1  = Heatmap(mat,name = "Relation",
               column_split = sp_cols,
               row_split =  split_rows, row_gap = unit(2, "mm"),border = TRUE,
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 90,
               row_names_gp = gpar(fontsize = 12,fontface = "italic"),
               left_annotation = ha1,
               show_parent_dend_line = F,
               #top_annotation = ha1,
               cluster_columns = F,
               #top_annotation = col_ha,
               rect_gp = gpar(col = "gainsboro"),
               row_names_side = "left",
               row_dend_side = c("left"),
               col = col_binary,
               
               row_names_max_width = max_text_width(rownames(mat),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat),
                                                        gp = gpar(fontsize = 12)),
               column_names_rot =45,
               
               heatmap_legend_param = list(
                 title = "Relations", at = c(-1, 1), 
                 labels = c("Neg", "Pos")
               ) )


pdf(paste0("../Figs/RF_FC_CC/Heatmap_imp_Mic_H_Pathway.pdf"),height = 14, width = 10)
draw(ht1)
dev.off()



# Save the data from the heatmap

# Matrix
mat

# Row Annotation
rowann_dt <- cbind(match_tax_dt)

# Column Annotation
colann_dt <-as.data.frame(p_ann)

# Full df
f_df <-  cbind(rowann_dt,mat)

write.csv(f_df,"../Figs/RF_FC_CC/Fig_7_B_data.csv")

write.csv(colann_dt,"../Figs/RF_FC_CC/Fig_7_B_col_ann.csv")



# Fisher Exact test

# Check if the Clostridia are overrepresented predicting pathways


# Total number of clostridia in the model

tax_dt <- tax_table(phy_vst_mic)@.Data
tax_dt <- data.frame(tax_dt)
tab_freq  <- data.frame(table(tax_dt$Order))

# Clostridiales
num_clos <-  tab_freq$Freq[tab_freq$Var1 == "Clostridiales"]
num_total <-  nrow(tax_dt)

tax_after_dt <- data.frame(tax_table(phy_ht)@.Data)
tab_freq_after  <- data.frame(table(tax_after_dt$Order))

num_clos_after <- tab_freq_after$Freq[tab_freq_after$Var1 == "Clostridiales"]
num_total_after <- nrow(tax_after_dt)


num_clos_after/num_total_after
num_clos/num_total

num_clos_after <- 30

chisq.test(matrix(c(num_clos, num_total-num_clos,
                    num_clos_after, num_total_after-num_clos_after), ncol=2))







