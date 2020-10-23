#Load libraries
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(gtools)
library(ComplexHeatmap)
library(circlize)
library(variancePartition)
library(edgeR)
library(BiocParallel)
library(ggthemes)

# Create a directory to save figures and tables

mainDir <- "../Figs"
subDir <- "Long_Obs_analysis"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Long_Obs_analysis"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


########## Read Microbiome data #################################################

phy_pair <- readRDS("../data/Long_Obs/phy_mic_long_obs.rds")

# Now we are ready to analyze the data

############## Fig 5  #########################


## TTP analysis#########

ttp_data <- data.frame(sample_data(phy_pair))


ttp_data$Average.TTP <-  as.numeric(as.character(ttp_data$Average.TTP))

names(ttp_data)

ttp_data <- ttp_data[, c("Patient.ID","Time","Average.TTP","TTP1" ,"TTP2")]


# Similar to NTZ HRZE
ttp_data <- ttp_data[!is.na(ttp_data$Average.TTP),]

ttp_data$Day <-  as.character(ttp_data$Time)
ttp_data$Day[ttp_data$Time == "Day0"] <- 0
ttp_data$Day[ttp_data$Time == "Day14"] <- 14
ttp_data$Day[ttp_data$Time == "OneMonth"] <- 30
ttp_data$Day[ttp_data$Time == "TwoMonths"] <- 56

ttp_data$Day <-  as.numeric(ttp_data$Day)

#ttp_data <- ttp_data[ttp_data$Day %in% c(0, 14),]
ttp_data$Day <-  factor((ttp_data$Day))

p <- ggplot()+
  geom_point(data=ttp_data, aes(x = Day, y=Average.TTP, group = Patient.ID),
             color= "#984ea3", fill = "#984ea3", 
             size=2, shape = 16, alpha =0.7)+
  geom_line(data=ttp_data, aes(x = Day, y=Average.TTP, group = Patient.ID),
            color= "#984ea3", 
            size=0.5, alpha = 0.8)+
  theme_classic()+
  xlab("Day")+
  ylab("Time to positivity (Hours)")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        #axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15)) 
#+
#        scale_x_continuous(breaks=c(0,14,30,56),
#                     labels=c("0", "14", "30", "56"))


pdf(paste(fig_folder,"Long_Obs_TTP_Fig_5_A.pdf",sep = "/"), height = 4,width = 4,useDingbats = F)
print(p)
dev.off()

# Save csv file for each figure:
write.csv(ttp_data,paste(tab_folder,"Fig_5_A.csv",sep = "/"))


# LME Model 

# Model TTP ~ Treatment + Day + Treatment:Day
# Day -  0, 14
# Treatment - HRZE NTZ 
library(nlme)
lme_dt <-  ttp_data
#lme_dt$Day <- factor(lme_dt$Day, levels = c(0,14))
lme_dt$ID <- factor(lme_dt$Patient.ID)
lme_dt$Av_TTP <- log(lme_dt$Average.TTP)
names(lme_dt)

lme_dt <- lme_dt[!is.na(lme_dt$Av_TTP),]

fit_ttp <-lme( Av_TTP ~  Day ,
               random = ~ 1|ID,
               data = lme_dt)

sum_fit <- summary(fit_ttp)
dt_sum <-  data.frame(sum_fit$tTable)

#rownames(dt_sum)[4:5] <- c("0 NTZ - 0 HRZE","14 HRZE - 0 HRZE")

# The only comparisons missing is "TreatmentNTZ_Day_14_Vs_TreatmentNTZ_Day_0"

# It can be done by flipping the reference level in lme or using emmeans
# lme_dt$Treatment <- relevel(lme_dt$Treatment ,ref = "NTZ")
# fit_ttp_ntz <-lme( Av_TTP ~ sex + Age + Treatment * Day ,
#                random = ~ 1| ID,
#                data = lme_dt)
# 
# sum_fit_ntz <- summary(fit_ttp_ntz)
# dt_sum_ntz <-  data.frame(sum_fit_ntz$tTable)

# Using emmeans
library(emmeans)
em_sum <- emmeans(fit_ttp,specs = pairwise ~ Day,adjust="bonferroni")
contrast_dt  <-  data.frame(em_sum$contrasts)
rownames(contrast_dt) <- paste0("Contrast ",contrast_dt$contrast)
contrast_dt$contrast <- NULL

names(contrast_dt) <- names(dt_sum)

dt_sum <- rbind(dt_sum,contrast_dt)


# Save the LME results
write.csv(dt_sum,paste(tab_folder,"LME_TTP_5_A.csv",sep = "/"))


###########Diversity boxplot ###############
raw_phy <-  phy_pair
alpha <- estimate_richness(raw_phy) #%>% mutate(sample=row.names(.))
data <- get.samp(raw_phy,stats = T) %>% as.data.frame()




# Line plot similar to  TTP
div_shan <- data
div_shan$Day <-  as.character(div_shan$Time)

div_shan$Day[div_shan$Time == "Day0"] <- 0
div_shan$Day[div_shan$Time == "Day7"] <- 7
div_shan$Day[div_shan$Time == "Day14"] <- 14
div_shan$Day[div_shan$Time == "OneMonth"] <- 30
div_shan$Day[div_shan$Time == "TwoMonths"] <- 56
div_shan$Day[div_shan$Time == "SixMonthFollowup"] <- 180

div_shan$Day <-  as.numeric(div_shan$Day)
#div_shan <- div_shan[div_shan$Day != 180,] 
div_shan$Day <-  factor((div_shan$Day))


p_div <- ggplot()+
  
  geom_point(data=div_shan, aes(x = Day, y=InvSimpson, group = Patient.ID),
             color= "#984ea3", fill = "#984ea3", 
             size=2, shape = 16, alpha =0.8)+
  geom_line(data=div_shan, aes(x = Day, y=InvSimpson, group = Patient.ID),
            color= "#984ea3",
            size=0.5, alpha = 0.8)+
  theme_classic()+
  xlab("Day")+
  ylab("InvSimpson")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        #axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15)) 
#+
# scale_x_continuous(breaks=c(0,7,14,30,56,180),
#                   labels=c("0", "7","14", "30", "56","180"))  

#dev.off()
pdf(paste(fig_folder,"LinePlot_Inv_Simpson_diversity_Long_Obs_Fig_5_B.pdf",sep = "/"),height = 4,width = 4) #paired sample analysis 
print(p_div)
dev.off()

# Save csv file for each figure:
write.csv(div_shan,paste(tab_folder,"Fig_5_A_B.csv",sep = "/"))




# Inv Simpson 
lme_div_dt <- div_shan
#lme_div_dt$Day <- factor(lme_div_dt$Day,levels = c("0", "14"))
lme_div_dt$ID <- factor(lme_div_dt$Patient.ID)

lme_div_dt

# Missing sex for repeated patient ID
# Need to fix that
sex_dt <-  lme_div_dt[lme_div_dt$Sex != "", ]
sex_dt <- sex_dt[,c("Patient.ID","Sex")]
lme_div_dt$Sex <- NULL


lm_dt <-  merge(lme_div_dt,sex_dt,by = "Patient.ID")

#lm_dt <- lm_dt[,c("Patient.ID")]

#lm_dt <-  lm_dt[is.na(lm_dt$Sex.y]
fit_div <-lme( InvSimpson ~ Sex +  Day ,
               random = ~ 1| ID,
               data = lm_dt,na.action = )
sum_fit <- summary(fit_div)
dt_sum <-  data.frame(sum_fit$tTable)
dt_sum



# Using emmeans
library(emmeans)
em_sum <- emmeans(sum_fit,specs = pairwise ~ Day,adjust="bonferroni")
contrast_dt  <-  data.frame(em_sum$contrasts)
rownames(contrast_dt) <- paste0("Contrast ",contrast_dt$contrast)
contrast_dt$contrast <- NULL

names(contrast_dt) <- names(dt_sum)

dt_sum <- rbind(dt_sum,contrast_dt)


# Save the LME results
write.csv(dt_sum,paste(tab_folder,"LME_InvSimpson_div_Fig_5_B.csv",sep = "/"))



################# Differential Analysis##############################

################ Microbiome #####################################
phy_trt <-  phy_pair

sum_counts <- c(10,seq(100,500,100))

cnt_filt <- 100

# Loop over the number of filtered counts

#for(cnt_filt in sum_counts){


phy_sel <- phy_trt
phy_sel <-  prune_taxa(taxa_sums(phy_sel)>0,phy_sel)

# Differential analysis via Limma
cnt_dt <- as.matrix(data.frame(otu_table(phy_sel)))
sm_dt <-  data.frame(sample_data(phy_sel))
sm_dt$ind <- factor(sm_dt$ind)
sm_dt$Time 
sm_dt$Sex
cnt_dt <-  cnt_dt[,match(rownames(sm_dt),colnames(cnt_dt))]

dim(cnt_dt)
# filter genes by number of counts
sum(rowSums(cnt_dt) >= cnt_filt)
isexpr = rowSums(cnt_dt) >= cnt_filt
table(rowSums(cnt_dt) >=cnt_filt)
# Standard usage of limma/voom
gExpr = DGEList( cnt_dt[isexpr,] )
gExpr = calcNormFactors( gExpr )
dim(gExpr$counts)
# filter genes by number of counts
#table(rowSums(cpm(cnt_dt)>0.1)==8)
#isexpr = rowSums(cpm(cnt_dt)>0.1) >= 4

sm_dt$Time
sm_dt$ind
sm_dt$Sex


#Visualize
library(RColorBrewer)
lcpm <- cpm(gExpr$counts, log=TRUE)
# group <- paste0(sm_dt$Time,"_",sm_dt$ind,"_",sm_dt$SEX,"_",sm_dt$EXTRACTION)
group <- paste0(sm_dt$Time,"_",sm_dt$ind)
plotMDS(lcpm, labels=group)


# Limma 
# apply duplicateCorrelation is two rounds
design = model.matrix( ~  Time  , sm_dt)
# apply duplicateCorrelation is two rounds
vobj_tmp = voom(gExpr, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=sm_dt$ind)

# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
# Otherwise, use the results from the first voom run
pdf(paste0(tab_folder,"/Mean_Dispersion_plot_EBA.pdf"),width = 5,height = 5)
vobj = voom( gExpr, design, plot=T, block=sm_dt$ind, correlation=dupcor$consensus)
dev.off()

# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj, design, block=sm_dt$ind)

# But this step uses only the genome-wide average for the random effect
fitDupCor <- lmFit(vobj, design, block=sm_dt$ind, correlation=dupcor$consensus)

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )

# get names of available coefficients and contrasts for testing
colnames(fitDupCor)


topTable( fitDupCor, coef=c("TimeDay7"),n=10)
topTable( fitDupCor, coef=c("TimeDay14"),n=30)

p_cut_off <-  0.05 

colnames(fitDupCor)
# Create a contrast matrix for the required comparisons
cont.mat <- cbind(Day_7_Vs_Day_0 = c(0,1,0,0,0,0) ,
                  Day_14_Vs_Day_0 = c(0,0,1,0,0,0) ,
                  One_Month_Vs_Day_0 = c(0,0,0,1,0,0),
                  Two_Month_Vs_Day_0 = c(0,0,0,0,1,0),
                  Six_Month_Vs_Day_0 = c(0,0,0,0,0,1),
                  All_Vs_Day_0 = c(0,1,1,1,1,0))

final_fit <-  contrasts.fit(fitDupCor, cont.mat)
final_fit <- eBayes(final_fit)

sum_eba <- summary( decideTests(final_fit,p.value = p_cut_off,
                                method = "global",adjust.method = "BH"))
sum_dt <-  as.data.frame.matrix(sum_eba) 
write.csv(sum_dt ,paste(tab_folder,paste("Summary_mic_Long_Obs.csv",sep=""),sep = "/"))

# Global pval adjustment
nom_pval <-  final_fit$p.value
rnames <- rownames(nom_pval)

vec_pval <- c(nom_pval)
global_pval <- p.adjust(nom_pval,method = "BH")
global_mat <- data.frame(matrix(global_pval,ncol = ncol(nom_pval)))
names(global_mat) <- colnames(nom_pval)
rownames(global_mat) <-  rnames


colnames(final_fit)
# DE Species (DES) D0 Vs D7 
res_7_Vs_0 <- topTable( final_fit, coef=c("Day_7_Vs_Day_0"),n=Inf)
p_cut_off <-  0.05
res_7_Vs_0 <- data.frame(res_7_Vs_0)

# Match the order with global_mat
res_7_Vs_0 <- res_7_Vs_0[match(rownames(global_mat),rownames(res_7_Vs_0)),]

res_7_Vs_0$global_padj <- global_mat$Day_7_Vs_Day_0

res_7_Vs_0 <- res_7_Vs_0[!is.na(res_7_Vs_0$adj.P.Val),]
res_7_Vs_0 <- res_7_Vs_0[order(res_7_Vs_0$global_padj,decreasing = F),]
res_7_Vs_0_sig <- res_7_Vs_0[res_7_Vs_0$global_padj <= p_cut_off,]
write.csv(res_7_Vs_0,paste(tab_folder,paste("Day_7_Vs_Day_0","_",cnt_filt,".csv",sep=""),sep = "/"))


# DES D0 Vs D14 
res_14_Vs_0 <- topTable( final_fit, coef=c("Day_14_Vs_Day_0"),n=Inf)
p_cut_off <-  0.05
res_14_Vs_0 <- data.frame(res_14_Vs_0)

# Match the order with global_mat
res_14_Vs_0 <- res_14_Vs_0[match(rownames(global_mat),rownames(res_14_Vs_0)),]

res_14_Vs_0$global_padj <- global_mat$Day_14_Vs_Day_0
res_14_Vs_0 <- res_14_Vs_0[!is.na(res_14_Vs_0$adj.P.Val),]
res_14_Vs_0 <- res_14_Vs_0[order(res_14_Vs_0$global_padj,decreasing = F),]
res_14_Vs_0_sig <- res_14_Vs_0[res_14_Vs_0$global_padj <= p_cut_off,]
write.csv(res_14_Vs_0,paste(tab_folder,paste("Day_14_Vs_Day_0","_",cnt_filt,".csv",sep=""),sep = "/"))


# DEG D0 Vs OneMonth 
res_one_month_Vs_0 <- topTable( final_fit, coef=c("One_Month_Vs_Day_0"),n=Inf)
p_cut_off <-  0.05
res_one_month_Vs_0 <- data.frame(res_one_month_Vs_0)
# Match the order with global_mat
res_one_month_Vs_0 <- res_one_month_Vs_0[match(rownames(global_mat),rownames(res_one_month_Vs_0)),]

res_one_month_Vs_0$global_padj <- global_mat$One_Month_Vs_Day_0

res_one_month_Vs_0 <- res_one_month_Vs_0[!is.na(res_one_month_Vs_0$adj.P.Val),]
res_one_month_Vs_0 <- res_one_month_Vs_0[order(res_one_month_Vs_0$global_padj,decreasing = F),]
res_one_month_Vs_0_sig <- res_one_month_Vs_0[res_one_month_Vs_0$global_padj <= p_cut_off,]
write.csv(res_one_month_Vs_0,paste(tab_folder,paste("One_Month_Vs_Day_0_",cnt_filt,".csv",sep=""),sep = "/"))


# DEG D0 Vs TwoMonth 
res_two_month_Vs_0 <- topTable( final_fit, coef=c("Two_Month_Vs_Day_0"),n=Inf)
p_cut_off <-  0.05
res_two_month_Vs_0 <- data.frame(res_two_month_Vs_0)
# Match the order with global_mat
res_two_month_Vs_0 <- res_two_month_Vs_0[match(rownames(global_mat),rownames(res_two_month_Vs_0)),]

res_two_month_Vs_0$global_padj <- global_mat$Two_Month_Vs_Day_0


res_two_month_Vs_0 <- res_two_month_Vs_0[!is.na(res_two_month_Vs_0$adj.P.Val),]
res_two_month_Vs_0 <- res_two_month_Vs_0[order(res_two_month_Vs_0$global_padj,decreasing = F),]
res_two_month_Vs_0_sig <- res_two_month_Vs_0[res_two_month_Vs_0$global_padj <= p_cut_off,]
write.csv(res_two_month_Vs_0,paste(tab_folder,paste("Two_Month_Vs_Day_0_",cnt_filt,".csv",sep=""),sep = "/"))


#TimeSixMonthFollowup
# DEG D0 Vs SixMonth 
res_six_month_Vs_0 <- topTable( final_fit, coef=c("Six_Month_Vs_Day_0"),n=Inf)
p_cut_off <-  0.05
res_six_month_Vs_0 <- data.frame(res_six_month_Vs_0)
# Match the order with global_mat
res_six_month_Vs_0 <- res_six_month_Vs_0[match(rownames(global_mat),rownames(res_six_month_Vs_0)),]

res_six_month_Vs_0$global_padj <- global_mat$Six_Month_Vs_Day_0


res_six_month_Vs_0 <- res_six_month_Vs_0[!is.na(res_six_month_Vs_0$adj.P.Val),]
res_six_month_Vs_0 <- res_six_month_Vs_0[order(res_six_month_Vs_0$global_padj,decreasing = F),]
res_six_month_Vs_0_sig <- res_six_month_Vs_0[res_six_month_Vs_0$global_padj <= p_cut_off,]
write.csv(res_six_month_Vs_0,paste(tab_folder,paste("Time_SixMonth_Vs_0_",cnt_filt,".csv",sep=""),sep = "/"))

# ASVs that responds to treatment at any time except SixMonth
# DES (D7, D14, One Month ,Two month)  Vs D0

res_all_Vs_0 <- topTable( final_fit,coef = c("All_Vs_Day_0"), n=Inf)
p_cut_off <-  0.05
res_all_Vs_0 <- data.frame(res_all_Vs_0)
# Match the order with global_mat
res_all_Vs_0 <- res_all_Vs_0[match(rownames(global_mat),rownames(res_all_Vs_0)),]

res_all_Vs_0$global_padj <- global_mat$All_Vs_Day_0


res_all_Vs_0 <- res_all_Vs_0[!is.na(res_all_Vs_0$adj.P.Val),]
res_all_Vs_0 <- res_all_Vs_0[order(res_all_Vs_0$global_padj,decreasing = F),]
res_all_Vs_0_sig <- res_all_Vs_0[res_all_Vs_0$global_padj <= p_cut_off,]

write.csv(res_all_Vs_0,paste(tab_folder,paste("Time_All_Vs_0_",cnt_filt,".csv",sep=""),sep = "/"))



res_7_Vs_0
res_14_Vs_0
res_one_month_Vs_0
res_two_month_Vs_0
res_six_month_Vs_0


mat <- cpm(cnt_dt, log=TRUE)
phy_ht <- phy_pair
otu_table(phy_ht) <-  otu_table(mat, taxa_are_rows = T)
#common_otus <- intersect(rownames(res_hrze_post),rownames(res_ntz_post))
all_sp <- Reduce(union,list(rownames(res_7_Vs_0_sig),
                            rownames(res_14_Vs_0_sig),
                            rownames(res_one_month_Vs_0_sig),
                            rownames(res_two_month_Vs_0_sig),
                            rownames(res_six_month_Vs_0_sig)))
phy_ht <- prune_taxa(all_sp,phy_ht)

phy_melt <-  get.otu.melt(phy_ht)
tax_dt <-  unique(phy_melt[,c("otu","Order","Phylum")])

#tax_dt$Phylum[tax_dt$Phylum == "xxxx"] <- "Tenericutes" 
topN <- length(unique(tax_dt$Order))
match_tax_dt <- unique(tax_dt [,c("Order","Phylum")])
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]

phylum_col <- c("Firmicutes"= "#9C854E",
                "Bacteroidetes" = "#51AB9B",
                "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red",
                "Verrucomicrobia" = "#33445e",
                "Tenericutes" = "#e2eab7",
                "Fusobacteria"= "#653131",
                "Cyanobacteria" = "#808080",
                "Euryarchaeota" = "#8c8c8c",
                "Candidatus Melainabacteria" = "#999999",
                "Elusimicrobia" = "#777777",
                "Lentisphaerae" = "#555555",
                "Spirochaetes" = "#595959",
                "Streptophyta" = "#4c4c4c",
                "Synergistetes" = "#a6a6a6"
)

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


# Volcano plots for Day 7 Vs Day 0
# HRZE
fc_dt <-  res_7_Vs_0
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze

mer_hrze$Order <-  factor(mer_hrze$Order,levels = tax_sel$Order)

p_cut_off <- 0.05

# data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
# p1 <-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
#              aes(y=-log10(global_padj),x=logFC),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data= data_pvalue_less,
#              aes(y=-log10(global_padj),x=logFC,color=Order),
#              size=2,alpha=1)+
#   scale_color_manual(name="Order",values = mycol) +
#   geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
#   geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
#   geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
#   theme_base()+
#   guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
#   #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
#   xlab("Log2 Fold Change")+
#   ylab("-Log10 P-Value")+
#   xlim(-8,8)+
#   ylim(0,10)
# 
# print(p1)

p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
data_pvalue_less$otu_mod <- gsub('.*\\(',"",data_pvalue_less$otu)
data_pvalue_less$otu_mod <- gsub('\\)',"",data_pvalue_less$otu_mod)
data_pvalue_less$otu_mod <- gsub('\\_'," ",data_pvalue_less$otu_mod)
p1 <-ggplot() + 
  
  ggrepel::geom_text_repel(data= data_pvalue_less[abs(data_pvalue_less$logFC)>2.5,],
                           aes(y=-log10(global_padj),x=logFC,label=otu_mod))+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,fill="black",color='black',alpha=0.1,shape=16)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,fill=Order,color=Order),
             size=2,alpha=.9,shape=16)+
  scale_fill_manual(name="Order",values = mycol) +
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  guides(fill=FALSE,color = F)+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-8,8)+
  ylim(0,8)

print(p1)


# Volcano plots for Day 14 Vs Day 0
# HRZE
fc_dt <-  res_14_Vs_0
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze

mer_hrze$Order <-  factor(mer_hrze$Order,levels = tax_sel$Order)

p_cut_off <- 0.05

# data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
# p2 <-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
#              aes(y=-log10(global_padj),x=logFC),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data= data_pvalue_less,
#              aes(y=-log10(global_padj),x=logFC,color=Order),
#              size=2,alpha=1)+
#   scale_color_manual(name="Order",values = mycol) +
#   geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
#   geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
#   geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
#   theme_base()+
#   guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
#   #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
#   xlab("Log2 Fold Change")+
#   ylab("-Log10 P-Value")+
#   xlim(-8,8)+
#   ylim(0,10)
# 
# print(p2)


p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
data_pvalue_less$otu_mod <- gsub('.*\\(',"",data_pvalue_less$otu)
data_pvalue_less$otu_mod <- gsub('\\)',"",data_pvalue_less$otu_mod)
data_pvalue_less$otu_mod <- gsub('\\_'," ",data_pvalue_less$otu_mod)
p2 <-ggplot() + 
  ggrepel::geom_text_repel(data= data_pvalue_less[abs(data_pvalue_less$logFC)>2.5,],
                           aes(y=-log10(global_padj),x=logFC,label=otu_mod))+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,fill="black",color='black',alpha=0.1,shape=16)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,fill=Order,color=Order),
             size=2,alpha=.9,shape=16)+
  scale_fill_manual(name="Order",values = mycol) +
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  guides(fill=FALSE,color = F)+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-8,8)+
  ylim(0,8)

print(p2)


# Volcano plots for One Month Vs Day 0
# HRZE
fc_dt <-  res_one_month_Vs_0
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze
mer_hrze$Order <-  factor(mer_hrze$Order,levels = tax_sel$Order)


# p_cut_off <- 0.05
# 
# data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
# p3 <-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
#              aes(y=-log10(global_padj),x=logFC),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data= data_pvalue_less,
#              aes(y=-log10(global_padj),x=logFC,color=Order),
#              size=2,alpha=1)+
#   scale_color_manual(name="Order",values = mycol) +
#   geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
#   geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
#   geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
#   theme_base()+
#   guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
#   #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
#   xlab("Log2 Fold Change")+
#   ylab("-Log10 P-Value")+
#   xlim(-8,8)+
#   ylim(0,10)
# 
# print(p3)
p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
data_pvalue_less$otu_mod <- gsub('.*\\(',"",data_pvalue_less$otu)
data_pvalue_less$otu_mod <- gsub('\\)',"",data_pvalue_less$otu_mod)
data_pvalue_less$otu_mod <- gsub('\\_'," ",data_pvalue_less$otu_mod)
p3 <-ggplot() + 
  ggrepel::geom_text_repel(data= data_pvalue_less[abs(data_pvalue_less$logFC)>2.5,],
                           aes(y=-log10(global_padj),x=logFC,label=otu_mod))+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,fill="black",color='black',alpha=0.1,shape=16)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,fill=Order,color=Order),
             size=2,alpha=.9,shape=16)+
  scale_fill_manual(name="Order",values = mycol) +
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  guides(fill=FALSE,color = F)+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-8,8)+
  ylim(0,8)


print(p3)




# Volcano plots for Two Month Vs Day 0
# HRZE
fc_dt <-  res_two_month_Vs_0
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze
mer_hrze$Order <-  factor(mer_hrze$Order,levels = tax_sel$Order)


# p_cut_off <- 0.05
# 
# data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
# p4 <-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
#              aes(y=-log10(global_padj),x=logFC),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data= data_pvalue_less,
#              aes(y=-log10(global_padj),x=logFC,color=Order),
#              size=2,alpha=1)+
#   scale_color_manual(name="Order",values = mycol) +
#   geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
#   geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
#   geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
#   theme_base()+
#   guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
#   #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
#   xlab("Log2 Fold Change")+
#   ylab("-Log10 P-Value")+
#   xlim(-8,8)+
#   ylim(0,10)
# 
# print(p4)
p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
data_pvalue_less$otu_mod <- gsub('.*\\(',"",data_pvalue_less$otu)
data_pvalue_less$otu_mod <- gsub('\\)',"",data_pvalue_less$otu_mod)
data_pvalue_less$otu_mod <- gsub('\\_'," ",data_pvalue_less$otu_mod)
p4 <-ggplot() + 
  ggrepel::geom_text_repel(data= data_pvalue_less[abs(data_pvalue_less$logFC)>2.5,],
                           aes(y=-log10(global_padj),x=logFC,label=otu_mod))+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,fill="black",color='black',alpha=0.1,shape=16)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,fill=Order,color=Order),
             size=2,alpha=.9,shape=16)+
  scale_fill_manual(name="Order",values = mycol) +
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  guides(fill=FALSE,color = F)+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-8,8)+
  ylim(0,8)

print(p4)



# Volcano plots for SixMonth Vs Day 0
# HRZE
fc_dt <-  res_six_month_Vs_0
fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze
mer_hrze$Order <-  factor(mer_hrze$Order,levels = tax_sel$Order)


# p_cut_off <- 0.05
# 
# data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
# p5 <-ggplot() + 
#   geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
#              aes(y=-log10(global_padj),x=logFC),
#              size=2,color="black",alpha=0.1)+
#   geom_point(data= data_pvalue_less,
#              aes(y=-log10(global_padj),x=logFC,color=Order),
#              size=2,alpha=1)+
#   scale_color_manual(name="Order",values = mycol) +
#   geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
#   geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
#   geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
#   theme_base()+
#   guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
#   #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
#   xlab("Log2 Fold Change")+
#   ylab("-Log10 P-Value")+
#   xlim(-8,8)+
#   ylim(0,10)
# 
# print(p5)

p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
data_pvalue_less$otu_mod <- gsub('.*\\(',"",data_pvalue_less$otu)
data_pvalue_less$otu_mod <- gsub('\\)',"",data_pvalue_less$otu_mod)
data_pvalue_less$otu_mod <- gsub('\\_'," ",data_pvalue_less$otu_mod)

p5 <-ggplot() + 
  ggrepel::geom_text_repel(data= data_pvalue_less[abs(data_pvalue_less$logFC)>2.5,],
                           aes(y=-log10(global_padj),x=logFC,label=otu_mod))+
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,fill="black",color='black',alpha=0.1,shape=16)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,fill=Order,color=Order),
             size=2,alpha=.9,shape=16)+
  scale_fill_manual(name="Order",values = mycol) +
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  guides(fill=FALSE,color = F)+
  theme_base()+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 P-Value")+
  xlim(-8,8)+
  ylim(0,8)

print(p5)


library(ggpubr)
combined<-ggarrange(p1,p2,p3,p4,p5)
combined
pdf(paste(fig_folder,"volcano_plot_combined_5C.pdf",sep = "/"),height = 10, width = 20, useDingbats = F)
print(combined)
dev.off()



# Now make a table of up and down ASVs in different comparisons
list_dt_sig <- list(res_7_Vs_0_sig,res_14_Vs_0_sig,res_one_month_Vs_0_sig,
                    res_two_month_Vs_0_sig,res_six_month_Vs_0_sig)


comb_dt_sig  <- bind_rows(list_dt_sig, .id = "column_label")

comb_dt_sig$column_label <- factor(comb_dt_sig$column_label)
levels(comb_dt_sig$column_label) <- c("Day_7_Vs_0","Day_14_Vs_0","one_month_Vs_0",
                                      "two_month_Vs_0","six_month_Vs_0")











# Combine all the plots

list_dt <- list(res_7_Vs_0,res_14_Vs_0,res_one_month_Vs_0,
                res_two_month_Vs_0,res_six_month_Vs_0)


comb_dt  <- bind_rows(list_dt, .id = "column_label")
comb_dt$otu <- c(rownames(res_7_Vs_0),rownames(res_14_Vs_0),rownames(res_one_month_Vs_0),
                 rownames(res_two_month_Vs_0),rownames(res_two_month_Vs_0))
# HRZE
fc_dt <-  comb_dt
#fc_dt$otu <- rownames(fc_dt)
mer_hrze <- merge(tax_dt,fc_dt,by = "otu")
mer_hrze$Type <- "HRZE"
mer_hrze
mer_hrze$Order <-  factor(mer_hrze$Order,levels = tax_sel$Order)

mer_hrze$column_label[mer_hrze$column_label == 1 ] <- "Day 7 Vs Day 0"
mer_hrze$column_label[mer_hrze$column_label == 2 ] <- "Day 14 Vs Day 0"
mer_hrze$column_label[mer_hrze$column_label == 3 ] <- "Day 30 Vs Day 0"
mer_hrze$column_label[mer_hrze$column_label == 4 ] <- "Day 56 Vs Day 0"
mer_hrze$column_label[mer_hrze$column_label == 5 ] <- "Day 180 Vs Day 0"

mer_hrze$column_label <- factor(mer_hrze$column_label,levels = c("Day 7 Vs Day 0",
                                                                 "Day 14 Vs Day 0",
                                                                 "Day 30 Vs Day 0",
                                                                 "Day 56 Vs Day 0",
                                                                 "Day 180 Vs Day 0"))
p_cut_off <- 0.05

data_pvalue_less <- mer_hrze[which(mer_hrze$global_padj< p_cut_off),]
p5 <-ggplot() + 
  geom_point(data=mer_hrze[which(mer_hrze$global_padj>=p_cut_off),],
             aes(y=-log10(global_padj),x=logFC),
             size=2,color="black",alpha=0.1)+
  geom_point(data= data_pvalue_less,
             aes(y=-log10(global_padj),x=logFC,color=Order),
             size=2,alpha=1)+
  
  scale_color_manual(name="Order",values = mycol) +
  geom_hline(yintercept = -log10(p_cut_off),size = .1,linetype = "dashed")+
  geom_vline(xintercept = -1,size = .1,linetype = "dashed")+
  geom_vline(xintercept = 1,size = .1,linetype = "dashed")+
  theme_base()+
  facet_wrap(~column_label)+
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 7)))+
  #scale_fill_manual(name="",values = mycol[[1]],breaks=mycol[[2]]) +
  xlab("Log2 Fold Change")+
  ylab("-Log10 FDR")+
  xlim(-8,8)+
  ylim(0,8)

print(p5)

# Using the cowplot package
legend <- cowplot::get_legend(p5)

pdf(paste(fig_folder,"Fig_5C_legend.pdf",sep = "/"),height = 6, width = 6,
    useDingbats = F)
grid.newpage()
grid.draw(legend)
dev.off()

pdf(paste(fig_folder,"volcano_plot_combined_mic_Long_obs_nolabel.pdf",sep = "/"),height = 6, width = 14,
    useDingbats = F)
print(p5)
dev.off()


write.csv(mer_hrze, paste(tab_folder,paste("mic_volcano_plot_Fig_5_C",".csv",sep=""),sep = "/"))


# Visualize the species in in the heatmap
# Lets have a heatmap of the species across different time points
# Heatmaps for the differential species


phy_sel <- phy_ht                                            

mat <-data.frame(otu_table(phy_sel))


match_tax_dt <- data.frame(tax_table(phy_sel))
match_tax_dt <- match_tax_dt[order(match_tax_dt$Phylum,match_tax_dt$Order),]
rownames(match_tax_dt)


phylum_col <- c("Firmicutes"= "#9C854E","Bacteroidetes" = "#51AB9B",   "Actinobacteria" = "#A77097",
                "Proteobacteria" = "red","Tenericutes" = "brown","Lentisphaerae" = "purple",
                "Euryarchaeota" = "green","Verrucomicrobia" = "orange",
                "Fusobacteria"= "#653131",
                "Candidatus Melainabacteria" = "#999999",
                "Spirochaetes" = "#9790ab"
)
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

order_col

match_tax_dt <- match_tax_dt[match(rownames(mat),rownames(match_tax_dt)),]


match_tax_dt$Day_7_Vs_0 <- "NS"
match_tax_dt$Day_7_Vs_0[rownames(match_tax_dt) %in%
                          rownames(res_7_Vs_0_sig[res_7_Vs_0_sig$logFC > 0,])] <- "Pos"
match_tax_dt$Day_7_Vs_0[rownames(match_tax_dt) %in%
                          rownames(res_7_Vs_0_sig[res_7_Vs_0_sig$logFC < 0,])] <- "Neg"

match_tax_dt$Day_14_Vs_0 <- "NS"
match_tax_dt$Day_14_Vs_0[rownames(match_tax_dt) %in% 
                           rownames(res_14_Vs_0_sig[res_14_Vs_0_sig$logFC > 0,])] <- "Pos"
match_tax_dt$Day_14_Vs_0[rownames(match_tax_dt) %in% 
                           rownames(res_14_Vs_0_sig[res_14_Vs_0_sig$logFC < 0,])] <- "Neg"



match_tax_dt$One_month_Vs_0 <- "NS"
match_tax_dt$One_month_Vs_0[rownames(match_tax_dt) %in% 
                              rownames(res_one_month_Vs_0_sig[res_one_month_Vs_0_sig$logFC > 0,])] <- "Pos"
match_tax_dt$One_month_Vs_0[rownames(match_tax_dt) %in% 
                              rownames(res_one_month_Vs_0_sig[res_one_month_Vs_0_sig$logFC < 0,])] <- "Neg"


match_tax_dt$Two_month_Vs_0 <- "NS"
match_tax_dt$Two_month_Vs_0[rownames(match_tax_dt) %in% 
                              rownames(res_two_month_Vs_0_sig[res_two_month_Vs_0_sig$logFC > 0,])] <- "Pos"
match_tax_dt$Two_month_Vs_0[rownames(match_tax_dt) %in% 
                              rownames(res_two_month_Vs_0_sig[res_two_month_Vs_0_sig$logFC < 0,])] <- "Neg"

match_tax_dt$Six_month_Vs_0 <- "NS"
match_tax_dt$Six_month_Vs_0[rownames(match_tax_dt) %in% 
                              rownames(res_six_month_Vs_0_sig[res_six_month_Vs_0_sig$logFC > 0,])] <- "Pos"
match_tax_dt$Six_month_Vs_0[rownames(match_tax_dt) %in% 
                              rownames(res_six_month_Vs_0_sig[res_six_month_Vs_0_sig$logFC < 0,])] <- "Neg"






#match_tax_dt$Six_month_Vs_0 <- "Sig"
col_binary <- c("NS" = "grey","Pos"="blue","Neg" = "red")
ha2 = HeatmapAnnotation( Day_7_Vs_0 = match_tax_dt$Day_7_Vs_0,
                         Day_14_Vs_0 = match_tax_dt$Day_14_Vs_0,
                         One_month_Vs_0 = match_tax_dt$One_month_Vs_0,
                         Two_month_Vs_0 = match_tax_dt$Two_month_Vs_0,
                         Six_month_Vs_0 = match_tax_dt$Six_month_Vs_0,
                         col = list(Day_7_Vs_0 = col_binary,
                                    Day_14_Vs_0 = col_binary,
                                    One_month_Vs_0 = col_binary,
                                    Two_month_Vs_0 =  col_binary,
                                    Six_month_Vs_0 = col_binary ),which = "row")


library(ComplexHeatmap)
ha1 = HeatmapAnnotation(Order = match_tax_dt$Order,
                        col = list(Order = order_col), which = "row")
df <- as.data.frame(sample_data(phy_sel)[,c("Patient.ID","Time")]) #"days_ON_HRZE","TB_status"

#d5b8da  #984ea3 #6a3672 #3c1f41 #1e0f20

ha_column = HeatmapAnnotation(Time =  df$Time,
                              col=list(Time = c("Day0" = "#537c4a","Day7"="#d5b8da",
                                                "Day14" = "#984ea3","OneMonth" = "#6a3672",
                                                "TwoMonths" = "#3c1f41","SixMonthFollowup" = "#3c1f41")))



split_rows <-  match_tax_dt$Phylum
split_rows <- factor(split_rows, levels= unique(split_rows))
library(circlize)
library(RColorBrewer) 

colfunc <- colorRampPalette(rev(brewer.pal(n = 7, name =
                                             "RdYlBu")))

ht_col <- colfunc(8)


# Split the column
split_cols<-  df$Time
split_cols <- factor(split_cols, levels= c(  "Day0","Day7","Day14","OneMonth", "TwoMonths","SixMonthFollowup"))


# Remove underscores 
rownames(mat) <- gsub("_"," ",rownames(mat))

ht1 = Heatmap(as.matrix(mat), name = "vst", column_title = NA, 
              col = ht_col,
              #clustering_distance_rows = "euclidean",
              split = split_rows,
              column_split = split_cols ,
              cluster_column_slices = F,
              cluster_row_slices = F,
              top_annotation = ha_column,
              row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              show_parent_dend_line =  F,
              left_annotation = ha2,
              right_annotation = ha1,
              row_names_side = "right", km=1, color_space = "LAB",
              #row_dend_side="right",
              #col=cols,
              #col= colorRamp2(c(-2, -1,0, 1, 2),teddy_cols),
              #heatmap_legend_param = list(),
              #clustering_method_columns = "ward.D",
              width=2, show_column_names= F,
              row_names_gp = gpar(fontsize = 9),
              #cluster_columns = T,cluster_rows = T,
              row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 12)),
              #column_names_max_height = max_text_width(colnames(mat_logic), gp = gpar(fontsize = 12)),
              #column_names_side = "top",
              na_col="white")




pdf(paste(fig_folder,"Heatmap_Long_obs_Species.pdf",sep = "/"),height = 45, width = 20,useDingbats = F)
draw(ht1)
dev.off()



# Save the data from the heatmap

# Matrix
mat

# Row Annotation
rowann_dt <- match_tax_dt

# Column Annotation
colann_dt <-as.data.frame(df)

# Full df
f_df <-  cbind(rowann_dt,mat)

write.csv(f_df,paste(tab_folder,paste("Heatmap_Long_obs_Species.csv",sep=""),sep = "/"))

write.csv(colann_dt,paste(tab_folder,paste("Heatmap_Long_obs_Species_col_ann.csv",sep=""),sep = "/"))



