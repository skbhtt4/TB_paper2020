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



# Create a directory to save figures and tables

mainDir <- ".."
subDir <- "GSVA"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")



#RNASeq 
rna_phy <-  readRDS("../data/genes/phy_rna_all.rds")

rna_phy <-  prune_taxa(taxa_sums(rna_phy)>0, rna_phy)

phy_rna_long <-  readRDS("../data/Long_Obs/phy_rna_tb.rds")
sm_dt_long <-  data.frame(sample_data(phy_rna_long))

sm_dt_long <- sm_dt_long[,c("Patient.ID","SEX","GROUP")]
names(sm_dt_long) <- c("Patient.ID","sex","TB_status")
sample_data(phy_rna_long) <-  sample_data(sm_dt_long)

ps_comb <-  merge_phyloseq(rna_phy,phy_rna_long)

tmp <- data.frame(sample_data(ps_comb))

phy_gene_sel <-  ps_comb


dds <- phyloseq_to_deseq2(phy_gene_sel , ~  TB_status) #replace this with any sample variable(s)

dds$TB_status
dds <- estimateSizeFactors(dds,"poscounts")


# Takes a while to run
dds <- estimateDispersions(dds)
dds <- DESeq(dds,fitType= "local")

# VST phyloseq for Heatmap later
phy_vst_rna <- phy_gene_sel
vst_dt <- getVarianceStabilizedData(dds)
otu_table(phy_vst_rna) <-  otu_table(vst_dt, taxa_are_rows = T)
saveRDS(phy_vst_rna,paste0(fig_folder,"/phy_rna_vst_tb_status_all.rds"))


phy_vst_rna <-  readRDS("../data/GSVA/phy_rna_vst_tb_status_all.rds")

#phy_gsva <- rna_phy
phy_gsva <- phy_vst_rna







#I think this should be TPM--or some other normalized count 
counts <- otu_table(phy_gsva) %>% as.data.frame %>% as.matrix() #or as.matrix() #TPM data here with gene rownames and sample colnames

#metadata from RNAseq phyloseq object
metadata <- data.frame(sample_data(phy_gsva))
#%>% column_to_rownames(var = "sample") #data frame with sample metadata with rownames = colnames of counts df


#gene lists
msigdbr_show_species()
m_df = msigdbr("Homo sapiens") %>% dplyr::filter(gs_cat == "H" )
m_t2g = m_df %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
m_t2g <- m_t2g[,c("gene_symbol","gs_name")]
#m_t2g$entrez_gene <- as.character(m_t2g$entrez_gene )
hall_set <- unstack(m_t2g)

rownames(counts) <- gsub(".*_","",rownames(counts))

library(GSVA)
set.seed(1057)
gsvaRes_ssgsea <- gsva(counts, hall_set,
                       min.sz=5, max.sz=500,
                       method = "gsva",
                       kcdf="Gaussian", mx.diff=TRUE, verbose=FALSE, parallel.sz=1
)

write.csv(gsvaRes_ssgsea,paste0(fig_folder,"/GSVA_all.csv"))

gsvaRes_ssgsea <-  read.csv("../GSVA/GSVA_all.csv")
rownames(gsvaRes_ssgsea) <- gsvaRes_ssgsea$X
gsvaRes_ssgsea$X <-  NULL
colnames(gsvaRes_ssgsea) <- gsub("^X","",colnames(gsvaRes_ssgsea))



phy_gsva <- phyloseq(otu_table(gsvaRes_ssgsea,taxa_are_rows = T),sample_data(metadata))
saveRDS(phy_gsva,paste0(fig_folder,"/phy_GSVA_all_hallmark.rds"))

met  <- metadata
mat  <- gsvaRes_ssgsea 

met$SampleID  <- rownames(metadata)
met <- met[match(colnames(gsvaRes_ssgsea),met$SampleID),]

#make phyloseq object out of gsva data
met_samp <- met
met_samp$sample  <- rownames(met_samp)
rownames(met_samp) <- NULL
gvsa_phy <- phyloseq(otu_table(mat,taxa_are_rows = T),set.samp(met_samp))
gvsa_phy

status_col <-  brewer.pal(length(unique(met$TB_status)), "Set1")
names(status_col) <- unique(unique(met$TB_status))

# status_col <- c("pretreatment" = "#537c4a","NTZ" = "#e46983", "HRZE" = "#984ea3",
#                 "family_contact" =	"#7F007F","community_control" = "#00C000",
#                 "EBA_0" = "#537c4a", "EBA_14" = "#984ea3", "EBA_56" = "#F781BF")
#ord_id <- order(met$type,met$Patient.ID)
library(gtools)
mixedsort(met$TB_status)
ha_column = HeatmapAnnotation(Status =  met$TB_status,
                              col=list(Status = status_col))

split_cols<-  met$TB_status
split_cols <- factor(split_cols, 
                     levels= c("pretreatment","EBA_0","NTZ","HRZE","EBA_14","EBA_56",
                               "community_control","family_contact"))


levels(split_cols) <- c("Clinical_Trial_pre","Long_obs_pre","NTZ_14","HRZE_14","Long_obs_14","Long_obs_56","Comm_control","Family_contact")






colfunc <- colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))
cols <- colfunc(8)
rownames(mat) <-  gsub("HALLMARK_","",rownames(mat))

library(tibble)
Hallmark_pathway_info <- read.csv("../data/Hallmark/Hallmark_pathway_annotation.csv")
rownames(Hallmark_pathway_info) <- Hallmark_pathway_info$Hallmark.Name

setdiff(rownames(mat) ,rownames(Hallmark_pathway_info))

Hallmark_pathway_info <-  Hallmark_pathway_info[match(rownames(mat),
                                                      rownames(Hallmark_pathway_info)),]
split_rows <- Hallmark_pathway_info$Process.Category

rownames(mat) <- gsub("_"," ",rownames(mat))
rownames(mat) <- str_to_sentence(rownames(mat))


ht1 = Heatmap(mat, name = "NES", column_title = NA, 
              top_annotation = ha_column,
              col = cols,
              cluster_column_slices = F,
              row_split = split_rows,
              row_names_side = "left",
              row_title_gp = gpar(fontsize = 10,fontface = "bold"), 
              row_title_rot = 0,
              column_split = split_cols,
              show_parent_dend_line = F,
              width=2, cluster_columns = T, 
              row_names_gp = gpar(fontsize = 9),
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              show_column_names = F, show_row_names = T,
              column_names_side = "bottom",na_col="white",
              show_row_dend = F)
#pdf("GSVA_raw.pdf",width = 14,height = 8)
pdf(paste0(fig_folder,"/GSVA_deseq_norm_vst.pdf"),width = 17,height = 9)
draw(ht1)
dev.off()

# Save the data from the heatmap

# Matrix
mat

# Row Annotation
rowann_dt <- cbind(Hallmark_pathway_info)

# Column Annotation
colann_dt <-as.data.frame(met)
colann_dt$Cat <- split_cols

# Full df
f_df <-  cbind(rowann_dt,mat)

write.csv(f_df,paste(fig_folder,paste("Fig_7_A_data.csv",sep=""),sep = "/"))

write.csv(colann_dt,paste(fig_folder,paste("Fig_7_A_col_ann.csv",sep=""),sep = "/"))



#  compare distances in the matrix across groups


phy_gsva

sm_dt <-  data.frame(sample_data(phy_gsva))

sm_dt$Group <- as.character(sm_dt$TB_status)

sm_dt$Group[sm_dt$Group %in% c("community_control" ,"family_contact")] <- "Control"
sm_dt$Group[sm_dt$Group %in% c("pretreatment")] <- "Clinical_trial_pre"
sm_dt$Group[sm_dt$Group %in% c("HRZE")] <- "HRZE_14"
sm_dt$Group[sm_dt$Group %in% c("NTZ")] <- "NTZ_14"
sm_dt$Group[sm_dt$Group %in% c("EBA_0")] <- "Long_Obs_0"
sm_dt$Group[sm_dt$Group %in% c("EBA_14")] <- "Long_Obs_14"
sm_dt$Group[sm_dt$Group %in% c("EBA_56")] <- "Long_Obs_56"


sm_dt$sam_ID <-  rownames(sm_dt)

sample_data(phy_gsva) <- sample_data(sm_dt)

# every group 

grps <-  unique(sm_dt$Group)[!unique(sm_dt$Group) %in% "Control"]


phy_gsva 

p = phy_gsva
m = "euclidean"
s = "sam_ID"
d = "Group"

# calc distances
wu = phyloseq::distance(p, m)
wu.m = reshape2::melt(as.matrix(wu))

# remove self-comparisons
wu.m = wu.m %>%
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor,as.character)

# get sample data (S4 error OK and expected)
sd = data.frame(sample_data(p)) %>%
  select(s, d) %>%
  mutate_if(is.factor,as.character)

# combined distances with sample data
colnames(sd) = c("Var1", "Type1")
wu.sd = left_join(wu.m, sd, by = "Var1")

colnames(sd) = c("Var2", "Type2")
wu.sd = left_join(wu.sd, sd, by = "Var2")

# Only select the distance between control and rest of the samples
wu.sd <-  wu.sd[wu.sd$Type1 == "Control",]

wu.sd <-  wu.sd[wu.sd$Type2 != "Control",]


wu.sd$Type2 <- factor(wu.sd$Type2,levels = c("Control","HRZE_14","Long_Obs_14","Long_Obs_56","Clinical_trial_pre", "NTZ_14","Long_Obs_0" ))

# plot
p_box <- ggplot(wu.sd, aes(x = Type2, y = value)) +
  theme_bw() +
  #geom_point() +
  #  geom_violin()+
  geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
  scale_color_identity() +
  facet_wrap(~ Type1, scales = "free_x") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle(paste0("Distance Metric = ", m))+
  xlab("Groups Vs Control")+
  ylab("Distance")

p_box
# pdf(paste0(fig_folder,"/Boxplot_euc_dist_control.pdf"),width = 6,height = 4,useDingbats = F)
# p_box
# dev.off()


# Lets just look at clinical trial

clinical_dt <-  wu.sd[wu.sd$Type2 %in% c("Clinical_trial_pre","HRZE_14","NTZ_14"),]


# plot
p_box <- ggplot(clinical_dt, aes(x = Type2, y = value)) +
  theme_bw() +
  #geom_point() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  #geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
  scale_color_identity() +
  facet_wrap(~ Type1, scales = "free_x") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle(paste0("Distance Metric = ", m))+
  xlab("Groups Vs Control")+
  ylab("Distance")

p_box
pdf(paste0(fig_folder,"/Boxplot_euc_dist_clinical_trial.pdf"),width = 4,height = 4,useDingbats = F)
p_box
dev.off()









# Get the patient ID

clinical_sm_dt <- sm_dt[,c("sam_ID","Patient.ID")]

clincal_dt_mer <- merge(clinical_dt,clinical_sm_dt,by.x = "Var2", by.y = "sam_ID")

library(nlme)
lme_dt <-  clincal_dt_mer
lme_dt$Type2 <- relevel(lme_dt$Type2,ref = "Clinical_trial_pre")
lme_dt$ID <- factor(lme_dt$Patient.ID)
lme_dt$Treatment <- "HRZE"
hrze_pat <- unique(lme_dt$Patient.ID[lme_dt$Type2 %in% "HRZE_14"] )
ntz_pat <- unique(lme_dt$Patient.ID[lme_dt$Type2 %in% "NTZ_14"] )
lme_dt$Treatment[lme_dt$Patient.ID %in% ntz_pat] <- "NTZ"
lme_dt$Treatment <- factor(lme_dt$Treatment,levels =  c("NTZ","HRZE"))


lme_dt$Day <- "Pre"
lme_dt$Day[lme_dt$Type2 %in% c("HRZE_14","NTZ_14")] <- "Post"

lme_dt$Day <- factor(lme_dt$Day, levels = c("Pre","Post"))


#lme_dt$Av_TTP <- log(lme_dt$Av_TTP)
names(lme_dt)

unique(lme_dt$Treatment)

table(lme_dt$Treatment)
table(lme_dt$Day,lme_dt$Treatment)
# fit_ttp <-lme( Av_TTP ~ sex + Age + pool_batch + Treatment * Day ,
#                random = ~ 1|ID,
#                data = lme_dt)

fit_ttp <-lme( value ~ Type2 ,
               random = ~ 1|ID,
               data = lme_dt)


sum_fit <- summary(fit_ttp)
dt_sum <-  data.frame(sum_fit$tTable)
# 
# dt_sum
# Value  Std.Error   DF    t.value       p.value
# (Intercept)   3.0471183 0.09511607 2450  32.035787 2.152348e-188
# Type2HRZE_14 -0.4874600 0.04303907 2450 -11.325989  5.071016e-29
# Type2NTZ_14   0.2372406 0.03258294 2450   7.281129  4.438862e-13


# Save the LME results
write.csv(dt_sum,paste(fig_folder,"LME_Clinical.csv",sep = "/"))





# Lets use the long data
long_dt <-  wu.sd[wu.sd$Type2 %in% c("Long_Obs_0" ,"Long_Obs_14" ,"Long_Obs_56"),]

# plot
p_box <- ggplot(long_dt, aes(x = Type2, y = value)) +
  theme_bw() +
  #geom_point() +
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75))+
  #geom_jitter(height = 0, width = 0.1)+
  #geom_boxplot(aes(color = ifelse(Type1 == Type2, "red", "black"))) +
  scale_color_identity() +
  facet_wrap(~ Type1, scales = "free_x") +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) + 
  ggtitle(paste0("Distance Metric = ", m))+
  xlab("Groups Vs Control")+
  ylab("Distance")

p_box
pdf(paste0(fig_folder,"/Boxplot_euc_dist_long_obs.pdf"),width = 4,height = 4,useDingbats = F)
p_box
dev.off()
















# Get the patient ID

clinical_sm_dt <- sm_dt[,c("sam_ID","Patient.ID")]

clincal_dt_mer <- merge(long_dt,clinical_sm_dt,by.x = "Var2", by.y = "sam_ID")

library(nlme)
lme_dt <-  clincal_dt_mer
lme_dt$Type2 <- relevel(lme_dt$Type2,ref = "Long_Obs_0")
lme_dt$ID <- factor(lme_dt$Patient.ID)



#lme_dt$Av_TTP <- log(lme_dt$Av_TTP)
names(lme_dt)

table(lme_dt$Type2)
# fit_ttp <-lme( Av_TTP ~ sex + Age + pool_batch + Treatment * Day ,
#                random = ~ 1|ID,
#                data = lme_dt)

fit_ttp <-lme( value ~ Type2 ,
               random = ~ 1|ID,
               data = lme_dt)


sum_fit <- summary(fit_ttp)
dt_sum <-  data.frame(sum_fit$tTable)

# LME
# dt_sum
# Value  Std.Error   DF   t.value       p.value
# (Intercept)       2.84766175 0.08172139 2784 34.845979 4.128264e-221
# Type2Long_Obs_14 -0.17800675 0.02500896 2784 -7.117719  1.390260e-12
# Type2Long_Obs_56 -0.05686512 0.02898653 2784 -1.961777  4.988794e-02



# Save the LME results
write.csv(dt_sum,paste(fig_folder,"LME_Long_Obs.csv",sep = "/"))











