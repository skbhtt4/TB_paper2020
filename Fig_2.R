#Load libraries
library(DESeq2)
library(phyloseq)
library(dplyr)
library(yingtools2)
library(data.table)
library(ifultools)
library(stringr)
library(ggplot2)
library(gtools)



# Create a directory to save figures and tables

mainDir <- "../Figs"
subDir <- "Figs_2"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
fig_folder <- paste(mainDir,subDir,sep="/")

mainDir <- "../Tables"
subDir <- "Tabs_2"
dir.create(file.path(mainDir, subDir), showWarnings = TRUE,recursive = TRUE)
tab_folder <- paste(mainDir,subDir,sep="/")


########## Read Microbiome data #################################################

phy_pair <- readRDS("../data/mic/phy_mic_fil_paired.rds")

# Sample data
sm_dt <- data.frame(sample_data(phy_pair))

sm_dt$Treatment <- "HRZE"
sm_dt$Treatment[grep("NTZ",sm_dt$type)] <- "NTZ"


# Now we are ready to analyze the data

############## Fig 2  #########################

pal_Prt_HRZE_NTZ <- c("pretreatment" = "#537c4a",
                      "NTZ" = "#356983",
                      "HRZE" = "#984ea3") 
sm_dt$Treatment <- factor(sm_dt$Treatment,levels = c("HRZE","NTZ"))

ttp_dt <-  sm_dt[,c("Patient.ID","Day","Av_TTP","Treatment","TTP1", "TTP2","TTP_min","TTP_max", "Av_TTP" )]

ttp_dt <- ttp_dt[!is.na(ttp_dt$Av_TTP),]

set.seed(100)
p <- ggplot()+
  geom_errorbar(data=sm_dt, aes(x = Day, ymin = TTP_min, ymax= TTP_max, color= 
                                  Treatment, group = Patient.ID), width=.4, 
                position=position_dodge(.5),
                alpha = 0.9)+
  geom_point(data=ttp_dt, aes(x = Day, y=Av_TTP, color= Treatment,
                              group = Patient.ID, fill = Treatment), 
             position=position_dodge(.5), size=2, shape = 19, alpha =0.9)+
  geom_line(data=sm_dt, aes(x = Day, y=Av_TTP, color= 
                              Treatment, group = Patient.ID), 
            position=position_dodge(.5), 
            size=0.5, alpha = 0.9)+
  theme_classic()+
  scale_color_manual(values = pal_Prt_HRZE_NTZ)+
  scale_fill_manual(values = pal_Prt_HRZE_NTZ)+
  facet_wrap(~Treatment)+
  xlab("Day")+
  ylab("Time to positivity (Hours)")+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        #axis.text.x = element_text(size=15,angle=90, hjust=1),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15))+
  scale_x_continuous(breaks=c(0,14),
                     labels=c("0","14"))  


pdf(paste(fig_folder,"Fig_2_B.pdf",sep = "/"), height = 6,width = 8,useDingbats = F)
print(p)
dev.off()

# Save csv file for each figure:
write.csv(ttp_dt,paste(tab_folder,"Fig_2_B.csv",sep = "/"))


# LME Model 

# Model TTP ~ Treatment + Day + Treatment:Day
# Day -  0, 14
# Treatment - HRZE NTZ 
library(nlme)
lme_dt <-  sm_dt
lme_dt$Day <- factor(lme_dt$Day, levels = c(0,14))
lme_dt$ID <- factor(lme_dt$Patient.ID)
lme_dt$Treatment <- relevel(lme_dt$Treatment ,ref = "NTZ")
lme_dt$Av_TTP <- log(lme_dt$Av_TTP)
names(lme_dt)

lme_dt <- lme_dt[!is.na(lme_dt$Av_TTP),]

# fit_ttp <-lme( Av_TTP ~ sex + Age + pool_batch + Treatment * Day ,
#                random = ~ 1|ID,
#                data = lme_dt)

fit_ttp <-lme( Av_TTP ~ sex + Age  + Treatment * Day ,
               random = ~ 1|ID,
               data = lme_dt)


sum_fit <- summary(fit_ttp)
dt_sum <-  data.frame(sum_fit$tTable)

#rownames(dt_sum)[7:8] <- c("TreatmentNTZ_Day_0_Vs_TreatmentHRZE_Day_0","TreatmentHRZE_Day_14_Vs_TreatmentHRZE_Day_0")

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
em_sum <- emmeans(fit_ttp,specs = pairwise ~ Day:Treatment,adjust="bonferroni")
contrast_dt  <-  data.frame(em_sum$contrasts)
rownames(contrast_dt) <- paste0("Contrast ",contrast_dt$contrast)
contrast_dt$contrast <- NULL

names(contrast_dt) <- names(dt_sum)

dt_sum <- rbind(dt_sum,contrast_dt)


# Save the LME results
write.csv(dt_sum,paste(tab_folder,"LME_TTP_2_B.csv",sep = "/"))



##### Principle component analysis #####
library(vegan)
phy_trt <-  phy_pair
sm_dt <- data.frame(sample_data(phy_trt))
ps <- phy_trt
# Take a log transform 
pslog <- transform_sample_counts(ps, function(x) log(1 + x))


#Well first perform a PCoA using Bray-Curtis dissimilarity.
set.seed(1057)
out.pcoa.log <- ordinate(pslog,  method = "PCoA", distance = "bray")
evals <- out.pcoa.log$values[,1]
p <- plot_ordination(pslog, out.pcoa.log,axes = c(1,2)) 

data_pca <- p$data

data_pca$Group <-  gsub("Pre |Post ","",data_pca$type)

pal_sam_type <-  c("Pre NTZ" = "#537c4a","Pre HRZE" = "#537c4a",
                   "Post NTZ" = "#356983", "Post HRZE" = "#984ea3") 

data_pca[1,]
table(data_pca$type,data_pca$pool_batch)


# Fig 2 D
library(ggrepel)
p_pca <- ggplot(data_pca, aes(x =Axis.1, y = Axis.2,label = Patient.ID ))+ 
  geom_line(aes(group=Patient.ID),color="grey",size=0.5) + 
  geom_point(aes(color = type, group=Patient.ID, shape = drug),
             alpha = .8,size=5) + 
  stat_ellipse(data=data_pca,aes(color=type)) + #, linetype = 2
  # geom_text_repel()+
  theme_classic()+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15)) +
  #scale_color_manual(name = "",values=pal_Prt_HRZE_NTZ) +
  scale_color_manual(name = "",values=pal_sam_type) +
  scale_shape_manual(name = "",values = c("pretreatment" = 15,"HRZE" = 16,"NTZ" =  17))+
  xlab(p$labels$x)+
  ylab (p$labels$y)
pdf(paste(fig_folder,"PCoA_Fig_2_C.pdf",sep = "/"),width = 12,height = 10,useDingbats = F)
print(p_pca)
dev.off()


# Save the PCA
write.csv(data_pca,paste(tab_folder,"PCA_2_C.csv",sep = "/"))



library(vegan)
metadata <- as(sample_data(pslog), "data.frame")

p_manova  <- adonis(phyloseq::distance(pslog, method="bray") ~ type,
                    data = metadata)

write.csv(p_manova$aov.tab,paste(tab_folder,"permanova_Fig_2_C.csv",sep = "/"))




###########Diversity boxplot ###############
raw_phy <-  phy_pair
alpha <- estimate_richness(raw_phy) #%>% mutate(sample=row.names(.))
data <- get.samp(raw_phy,stats = T) %>% as.data.frame()

data$type <-  factor(data$type , levels = c("Pre HRZE","Post HRZE","Pre NTZ","Post NTZ"))

data$Treatment <- "HRZE"
data$Treatment[grep("NTZ",data$type)] <- "NTZ"
# LME Model diversity

# Inv Simpson 
lme_div_dt <- data
lme_div_dt$Day <- factor(lme_div_dt$Day,levels = c("0", "14"))
lme_div_dt$Treatment <- factor(lme_div_dt$Treatment,levels = c("HRZE", "NTZ"))
lme_div_dt$ID <- factor(lme_div_dt$Patient.ID)
lme_div_dt$Age

lme_sub_div_dt <- lme_div_dt[!is.na(lme_div_dt$Age),]

fit_div <-lme( InvSimpson ~ sex + Age + pool_batch +  Treatment * Day ,
               random = ~ 1| ID,
               data = lme_sub_div_dt)
sum_fit <- summary(fit_div)
dt_sum <-  data.frame(sum_fit$tTable)
dt_sum


#rownames(dt_sum)[7:8] <- c("TreatmentNTZ_Day_0_Vs_TreatmentHRZE_Day_0","TreatmentHRZE_Day_14_Vs_TreatmentHRZE_Day_0")


# Using emmeans
library(emmeans)
em_sum <- emmeans(sum_fit,specs = pairwise ~ Day:Treatment,adjust="bonferroni")
contrast_dt  <-  data.frame(em_sum$contrasts)
rownames(contrast_dt) <- paste0("Contrast ",contrast_dt$contrast)
contrast_dt$contrast <- NULL

names(contrast_dt) <- names(dt_sum)

dt_sum <- rbind(dt_sum,contrast_dt)


# Save the LME results
write.csv(dt_sum,paste(tab_folder,"LME_InvSimpson_div_1_D.csv",sep = "/"))




# Line plot InvSimpson
div_dt_inv <- lme_div_dt
div_dt_inv$Day_num <- as.numeric(as.character(div_dt_inv$Day))
# Line plot InvSimpson
p_div <- ggplot() +
  geom_point(data=div_dt_inv, aes(x = Day_num, y=InvSimpson, color= 
                                    Treatment, group = ID, fill = Treatment), size=3, shape = 19, alpha =0.7)+
  geom_line(data=div_dt_inv, aes(x = Day_num, y=InvSimpson, color= 
                                   Treatment, group = ID), 
            size=0.5, alpha = 0.7)+
  #stat_summary(data=div_dt_inv, aes(x = Day_num, y=InvSimpson, group = Day_num), fun = median, fun.min = median, fun.max = median,
  #             geom = "crossbar", width = 4)+
  theme_classic()+
  scale_color_manual(values = pal_Prt_HRZE_NTZ)+
  scale_fill_manual(values = pal_Prt_HRZE_NTZ)+
  facet_wrap(~Treatment)+
  xlab("Day")+
  ylab("InvSimpson")+
  #scale_x_continuous(breaks = c(0,5,10,15),labels = c(0,"","","14"))+
  theme(legend.text = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=15, hjust=1),
        axis.text.y = element_text(size=15),
        axis.ticks.x=element_blank())
pdf(paste(fig_folder,"LinePlot_InvSimpson_diversity_Fig_2_D.pdf",sep = "/"),height = 6,width = 8, useDingbats = F) #paired sample analysis 
print(p_div)
dev.off()


# Save csv file for each figure:
write.csv(div_dt_inv,paste(tab_folder,"Fig_2_D.csv",sep = "/"))

############## Fig 2 complete #########################
