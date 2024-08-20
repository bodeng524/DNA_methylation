library(knitr)
library(limma)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(RColorBrewer)
library(missMethyl)
library(minfiData)
library(Gviz)
library(DMRcate)
library(stringr)
library(dplyr)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
library(gridExtra)

install.packages("BiocManager")
BiocManager::install("minfi")
BiocManager::install("IlluminaHumanMethylationEPICmanifest")
BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
BiocManager::install("IlluminaHumanMethylationEPIC")
install.packages("RColorBrewer")
install.packages("FactoMineR")
BiocManager::install("methylationArrayAnalysis")
install.packages("DMRcate")
install.packages("Gviz")
install.packages("plotly")
install.packages("pheatmap")
install.packages("remotes")
remotes::install_version("estimability", version="1.4.1")
install.packages("estimability")
install.packages("FactoMineR")
install.packages("factoextra")
install.packages("data.table")
library(data.table)
install.packages("devtools")
install.packages("mgmtstp27_0.6-4.tar.gz",repos=NULL)
BiocManager::install("lumi")
install.packages("ade4")
install.packages("ggalluvial")
library(mgmtstp27)
BiocManager::install("ChAMP")
library(ChAMP)
#######做个一个samplesheet大杂烩（包含MGMT预测值，平均甲基化值，可能性最高的亚型，原复发等等）
mgmt_meanmethy_subgroup_samplesheet <- 
  merge(Global_Mvalue_allsample_66,subgroupswitch_sheet[,c("sample","subtype","score")],by.x="sample",by.y="sample")
############制作samplesheet
csv_files <- list.files(path=basedir,pattern=".csv",recursive = T,full.names = T)
all_data <- data.frame()
for (i in 1:length(csv_files)) {
  temp_data <- read.csv(csv_files[i],fill=TRUE, check.names=FALSE)
  all_data <- rbind(all_data, temp_data)
}

targets2 <- read.metharray.sheet(base="/home/r102673/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/EPIC2023-387-plate2/STS",pattern=".csv")
targets1 <- read.metharray.sheet(base="/home/r102673/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/EPIC2023-387-plate1/STS",pattern=".csv")
targets3 <- read.metharray.sheet(base="/home/r102673/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/MET2017-126-014/MET2017-126-014/STS",pattern=".csv")
targets4 <- read.metharray.sheet(base="/home/r102673/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays/MET2022-350-014",pattern=".csv")
targets <- bind_rows(targets1, targets2, targets3, targets4)
targets$Array_Slide <- paste(targets$Slide,targets$Array,sep = "_")
targets_unique <- targets %>%
  distinct(Array_Slide, .keep_all = TRUE)
targets_GSAM <- targets_unique %>%
  filter(grepl("GSAM", Sample_Name))
targets_GSAM$Group <- select(targets_GSAM$Group,select =-Sample_name)
targets_GSAM$Group <- targets_GSAM$Sample_Name %>% substr(6,9)
targets_GSAM$Sample_Name <- NULL
targets_GSAM <- targets_GSAM%>%rename(Sample_name=Group)
targets_GSAM$Sample_name <- targets_GSAM$Group$Sample_name 


MGMT_folder <- list.dirs(path="~/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131",recursive = FALSE,full.names = T)
MGMT_sheet <- data.frame()
for (i in 1:length(MGMT_folder)) {aaa <- MGMT_folder[i]%>% substr(142,171)
                                  MGMT_sheet <- rbind(MGMT_sheet,aaa)
}
MGMT_sheet<- MGMT_sheet%>%rename(name=X.206467010123_R01C01__GSAM_EBF2.)
MGMT_sheet$Sample_name <- MGMT_sheet$name%>%substr(27,31)

sample_names_MGMT <- as.vector(MGMT_sheet$Sample_name)
sample_names_targets <- as.vector(targets_GSAM$Sample_name)
# 查找在 MGMT_sheet 中但不在 targets_GSAM 中的值
diff_MGMT_to_targets <- setdiff(sample_names_MGMT, sample_names_targets)
# 查找在 targets_GSAM 中但不在 MGMT_sheet 中的值
diff_targets_to_MGMT <- setdiff(sample_names_targets, sample_names_MGMT)
targets_GSAM <- targets_GSAM %>% mutate("Group"=if_else(grepl("1",Sample_name),"primary","recurrence"))

targets_GSAM$Array_position <- paste(targets_GSAM$Slide,targets_GSAM$Array,sep="_")
targets_GSAM <- targets_GSAM[,-4]


#############Remove columns from rgSet that do not have a corresponding SampleName(eg.207331540058_R05C01)
basedir <- "/home/r102673/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays"

rgSet <- read.metharray.exp(base=basedir,recursive=T,force=T)
dim(rgSet)#1051539     103

keep_columns <- targets_GSAM$Array_position
filteredRgSet <- rgSet[, colnames(rgSet) %in% keep_columns]

all_idatfiles <- list.files(basedir, pattern = "\\.idat$", full.names = TRUE, recursive = TRUE)
filteredsample_names <- c("204073510041_R01C01", "204073510041_R02C01", "204073510041_R03C01", "204073510041_R04C01", "204073510041_R05C01",
                  "204073510041_R06C01", "204073510041_R07C01", "204073510041_R08C01", "204073520105_R01C01", "204073520105_R02C01",
                  "204073520105_R03C01", "204073520105_R04C01", "204073520105_R05C01", "204073520105_R06C01", "204073520105_R07C01",
                  "204073520105_R08C01")

# Create a regex pattern to match any of the sample names
pattern <- paste(filteredsample_names, collapse = "|")
# Filter all_idatfiles to find files that include any of the sample names
filtered_files <- grep(pattern, all_idatfiles, value = TRUE)

###############calculate P value
sessionInfo("minfi")
data("IlluminaHumanMethylationEPICmanifest")
detP <- detectionP(filteredRgSet)
head(detP)
dim(detP)#865859    103 #865859     87

pal <- brewer.pal(8,"Dark2")
par(mfrow=c(1,2))
barplot(colMeans(detP), col=pal[factor(targets_GSAM$Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets_GSAM$Group)), fill=pal, 
       bg="white")


qcReport(rgSet, sampNames=targets_GSAM$Array_position, sampGroups=targets_GSAM$Group, 
         pdf="qcReport.pdf")


keep <- colMeans(detP) < 0.01
table(keep)
filteredRgSet <- filteredRgSet[,keep]


# keep <- rowSums(detP< 0.01) == ncol(detP)
# table(keep)#FALSE   TRUE 13774 852085
# rgSet1 <- rgSet1[keep,]
# dim(rgSet1)#[1] 1035498     103


mSetSqNoob <- preprocessNoob(filteredRgSet)
dim(mSetSqNoob)# 850544    103 # 865859 103（第二次）#865859     87

keep <- rowSums(detP < 0.01) == ncol(mSetSqNoob) 
table(keep)#FALSE   TRUE 13774 852085
mSetSqFltNoob <- mSetSqNoob[keep,]
dim(mSetSqFltNoob)# 852085    103  #853868     87

EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
keep <- !(featureNames(mSetSqFltNoob) %in% EPICanno$Name[EPICanno$chr %in% 
                                                           c("chrX","chrY")])
table(keep)# FALSE   TRUE 18665 833420 
mSetSqFltNoob <- mSetSqFltNoob[keep,]
dim(mSetSqFltNoob)# 833420    103  #835134     87


mSetSqFltNoob <- mapToGenome(mSetSqFltNoob)
mSetSqFltNoob <- dropLociWithSnps(mSetSqFltNoob)
dim(mSetSqFltNoob)  ##807258     87


# exclude cross reactive probes 
dataDirectory <- ("/home/r102673/R/x86_64-redhat-linux-gnu-library/4.2/methylationArrayAnalysis/extdata")
xReactiveProbes <- read.csv(file=paste(dataDirectory,
                                       "48639-non-specific-probes-Illumina450k.csv",
                                       sep="/"), stringsAsFactors=FALSE)
keep1 <- !(featureNames(mSetSqFltNoob) %in% xReactiveProbes$TargetID)
table(keep1)#FALSE   TRUE  24752 790515 
mSetSqFltNoob <- mSetSqFltNoob[keep1,] 
dim(mSetSqFltNoob)#781174    103  ##782665     87
# aaa <- match(featureNames(mSetSqFltNoob),xReactiveProbes$TargetID)
# mSetSqFltNoob222 <- mSetSqFltNoob[is.na(aaa),]#这里可以用is.na函数把整数类型的aaa转换为逻辑值，然后应用到mSetSqFltNoob中。
# identical(mSetSqFltNoob222, mSetSqFltNoob)

mVals <- getM(mSetSqFltNoob)
head(mVals)
bVals <- getBeta(mSetSqFltNoob)
head(bVals)
dim(mVals)#781174    103  #782665     87

par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets_GSAM$Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets_GSAM$Group)),
       text.col=brewer.pal(8,"Dark2"))

densityPlot(mVals, sampGroups=targets_GSAM$Group,
            main="Normalized(Noob)", legend=FALSE)
legend("top", legend = levels(factor(targets_GSAM$Group)),
       text.col=brewer.pal(8,"Dark2"))



par(mfrow=c(1,2))
plotMDS(mVals, top=1000, gene.selection="common",
        col=pal[factor(targets_GSAM$Group)])
legend("top", legend=levels(factor(targets_GSAM$Group)), text.col=pal,
       bg="white", cex=0.7)


library(FactoMineR)
library(factoextra)
dat <- t(mVals)
#这个函数只能从行里面读取分组信息，所以一般情况下需要对bVals矩阵进行行列对调，然后经PCA函数处理后才可以用于这里函数中.
dat.pca <- PCA(dat, graph = F) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = targets_GSAM$Group, 
             addEllipses = TRUE, 
             legend.title = "Groups")

# 热图
library(pheatmap)
tmp=as.data.frame(mVals)
#这里的“1”指定要应用函数的维度,1 表示行，2 表示列.
#var 是 R 的内置函数，用于计算方差。
tmp$var=apply(tmp,1,var)
#order() 函数返回的是一个整数向量，表示各行的新位置。例如，如果 tmp$var 是 c(30, 10, 20)，order(tmp$var) 将返回 c(2, 3, 1)。
tmp=tmp[rev(order(tmp$var)), ]
tmp$var = NULL
tmpaaa <- tmp[1:1000,]

heatmap(as.matrix(tmp))

pheatmap::pheatmap(tmp[1:1000,],show_rownames = F, show_colnames = F)


designMatrix <- model.matrix(~ 0 + factor(targets_GSAM$Group))
colnames(designMatrix) <- levels(factor(c("primary", "recurrence")))
contrastMatrix <- makeContrasts(primaryVsrecurrence = primary - recurrence, levels = designMatrix)
fit <- lmFit(mVals, designMatrix)
fit2 <- contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit2)
logFC <- fit2$coefficients
head(logFC)
resultDataFrame <- topTable(fit2, coef="primaryVsrecurrence", number=Inf, sort.by="p", p.value=1)

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
cpg_info <- EPICanno[EPICanno$Name == "cg03732056",]



########################用原始数据计算MGMT甲基化的预测值
EPICanno <- as.data.table(EPICanno)
mgmt_probes <- data[UCSC_RefGene_Name == "MGMT", .(IlmnID = Name, geneSymbol = UCSC_RefGene_Name, chromosome = chr, mapInfo = pos)]
print(mgmt_probes)
cpg_ids <- mgmt_probes$IlmnID
target_ids <- c("cg12434587", "cg12981137")
id_presence <- target_ids %in% rownames(mSetSqFltNoob)

###The gene encoding MGMT is located at chromosome 10q26.3, where a CpG island of > 700 bp spans the promoter region and into the first exon.

##The CpG island constitutes 97 CpG dinucleotides, and two differentially methylated regions (DMR1 and DMR2) are defined as part of these CpG sites.

##For both DMRs, the methylation status correlates inversely with MGMT gene expression. DMR1 is upstream from the transcriptional start site (nucleotides -250 to -101), and DMR2 is situated within exon 1 (nucleotides +97 to +196). If any of the CpG dinucleotides within the DMR2 are changed, an effect is seen on MGMT expression, and if DMR2 is hypermethylated, DMR1 is too.

#这是计算两个cpg岛的M值的自定义函数
cpg_of_sample <- function(parameter1) 
{data.frame(
  t(
    getM(
      subset(mSetSqFltNoob, rownames(mSetSqFltNoob) %in% target_ids, colnames(mSetSqFltNoob) %in% parameter1))))
}

#######单独计算MGMT两个cpg岛的B值
#计算87个样本cg12434587beta值的
Bvalue_cg12434587 <- function(parameter1) 
{data.frame(
  t(
    getBeta(
      subset(mSetSqFltNoob, rownames(mSetSqFltNoob) %in% "cg12434587", colnames(mSetSqFltNoob) %in% parameter1))))
}
B_value_cg12434587 <- list()

for (i in 1:nrow(targets_GSAM)){aaa <- Bvalue_cg12434587(targets_GSAM$Array_position[i])
B_value_cg12434587 <- rbind(B_value_cg12434587 ,aaa)}
B_value_cg12434587$sentrixposition <- rownames(B_value_cg12434587)
rownames(B_value_cg12434587) <- NULL
#计算87个样本cg12981137beta值的
Bvalue_cg12981137 <- function(parameter1) 
{data.frame(
  t(
    getBeta(
      subset(mSetSqFltNoob, rownames(mSetSqFltNoob) %in% "cg12981137", colnames(mSetSqFltNoob) %in% parameter1))))
}
B_value_cg12981137 <- list()

for (i in 1:nrow(targets_GSAM)){aaa <- Bvalue_cg12981137(targets_GSAM$Array_position[i])
B_value_cg12981137 <- rbind(B_value_cg12981137 ,aaa)}
######修剪两个cpg岛的表格，添加原发复发等等信息列
B_value_cg12434587 <- merge(B_value_cg12434587, predict_MGMT_allsample[, c("Group", "Sample_name","sample_group","sample")], 
                                                   by.x ="sentrixposition" , by.y = "sample", all.x = TRUE)
B_value_cg12434587 <- rename("BETA_value","cg12434587")

recurrence_pred_MGMT_cg12434587 <- B_value_cg12434587 %>%
  filter(Group == "recurrence") %>%
  dplyr::select(cg12434587)

primary_pred_MGMT_cg12434587 <- B_value_cg12434587 %>%
  filter(Group == "primary") %>%
  dplyr::select(cg12434587)
###计算cg12434587的T值并画线图可视化
t_cg12434587 <- t.test(recurrence_pred_MGMT_cg12434587, primary_pred_MGMT_cg12434587)

ggplot(B_value_cg12434587, aes(x = Group, y = cg12434587, group = sample_group, label=Sample_name)) +
  geom_point(position = position_jitter(width = 0.02), aes(color = sample_group)) +
  geom_line(aes(color = sample_group), position = position_dodge(width = 0.02))+
  labs(title = "BETA Values by Sample and Recurrence Status(cg12434587)",
       x = "Recurrence Status",
       y = "Estimated Value") +
  geom_hline(yintercept=0.3582, linetype="dashed", 
             color = "red", size=0.25) +
  #geom_text(size=2.4) +
  geom_text(aes(label = Sample_name), vjust = -0.5, hjust = 0.5, size = 2.4, position = position_jitter(width = 0.02)) +
  theme_minimal()
###计算cg12981137的T值并画线图可视化(未完成)




###这里是计算两个岛一起的值
mgmt_data_allsample <- list()

for (i in 1:nrow(targets_GSAM)){aaa <- cpg_of_sample(targets_GSAM$Array_position[i])
                                 mgmt_data_allsample<- rbind(mgmt_data_allsample,aaa)}

install.packages("lumi")#‘lumi’, ‘ade4’, ‘methylumi’
BiocManager::install("methylumi")
install.packages("~/R/packages.bo/mgmtstp27_0.6-4.tar.gz",repos=NULL)
library(mgmtstp27)

predict_MGMT_allsample <- data.frame()
for(i in 1:nrow(mgmt_data_allsample)){
  single_row_df <- mgmt_data_allsample[i, , drop = FALSE]
  predict_MGMT_allsample <- rbind(MGMTpredict(single_row_df),predict_MGMT_allsample)}

predict_MGMT_allsample <- merge(predict_MGMT_allsample, targets_GSAM[, c("Group", "Sample_name","Array_position")], 
                                by.x = "sample", by.y = "Array_position", all.x = TRUE)

predict_MGMT_allsample <- predict_MGMT_allsample %>%
  mutate(sample_group = substr(Sample_name, 1, 3))
####计算原发和复发的MGMT预计值有没有统计学差异
# 将数据划分为recurrence组和primary组
recurrence_pred_MGMT <- predict_MGMT_allsample %>%
  filter(Group == "recurrence") %>%
  dplyr::select(pred)

primary_pred_MGMT <- predict_MGMT_allsample %>%
  filter(Group == "primary") %>%
  dplyr::select(pred)
# 进行t检验
t_test_MGMT_between_TandR <- t.test(recurrence_pred_MGMT, primary_pred_MGMT)


######读取已经存在的MGMT预测值
mgmt_files <- list.files(path="~/mnt/neuro-genomic-1-ro/gsam/DNA Methylation - EPIC arrays - MNP CNS classifier/brain_classifier_v12.8_sample_report__v1.1__131",pattern="mgmt.csv",recursive = T,full.names = T)

mgmt_exist<- data.frame()

for (i in 1:length(mgmt_files)) {
  temp_data <- read.csv(mgmt_files[i],fill=TRUE, check.names=FALSE)
  mgmt_exist <- rbind(mgmt_exist, temp_data)
}

mgmt_exist <- mgmt_exist[,-1]
mgmt_exist <- mgmt_exist %>% mutate(sample_name=substr(mgmt_files,142,171))
mgmt_exist <- mgmt_exist %>% mutate(Sample_name=substr(mgmt_exist$sample_name,27,30))
mgmt_exist <- mgmt_exist %>% mutate(Group=if_else(grepl("1",mgmt_exist$Sample_name),"primary","recurrence"))
mgmt_exist <- mgmt_exist %>% mutate(sample=substr(mgmt_exist$sample_name,1,19))
mgmt_exist <- mgmt_exist %>% mutate(sample_group=substr(mgmt_exist$Sample_name,1,3))
###计算从表格提取出来的MGMT预计值的原发与复发之间的统计学差异
recurrence_pred_MGMT_exist <- mgmt_exist %>%
  filter(Group == "recurrence") %>%
  dplyr::select(Estimated)

primary_pred_MGMT_exist <- mgmt_exist %>%
  filter(Group == "primary") %>%
  dplyr::select(Estimated)
# 进行t检验
t_test_MGMT_exist_between_TandR <- t.test(recurrence_pred_MGMT_exist, primary_pred_MGMT_exist)


########
plot1 <-ggplot(predict_MGMT_allsample, aes(x = Group, y = pred, group = sample_group, label=Sample_name)) +
  geom_point(position = position_jitter(width = 0.02), aes(color = sample_group)) +
  geom_line(aes(color = sample_group), position = position_dodge(width = 0.02))+
  labs(title = "Estimated Values by Sample and Recurrence Status",
       x = "Recurrence Status",
       y = "Estimated Value") +
  geom_hline(yintercept=0.3582, linetype="dashed", 
             color = "red", size=0.25) +
  #geom_text(size=2.4) +
  geom_text(aes(label = Sample_name), vjust = -0.5, hjust = 0.5, size = 2.4, position = position_jitter(width = 0.02)) +
  theme_minimal()

plot2 <-ggplot(mgmt_exist, aes(x = Group, y = Estimated, group = sample_group, label=Sample_name)) +
  geom_point(position = position_jitter(width = 0.02), aes(color = sample_group)) +
  geom_line(aes(color = sample_group), position = position_dodge(width = 0.02))+
  labs(title = "Estimated Values by Sample and Recurrence Status(heidelberg)",
       x = "Recurrence Status",
       y = "Estimated Value") +
  geom_hline(yintercept=0.3582, linetype="dashed", 
             color = "red", size=0.25) +
  #geom_text(size=2.4) +
  geom_text(aes(label = sample_group), vjust = -0.5, hjust = 0.5, size = 2.4, position = position_jitter(width = 0.02)) +
  theme_minimal()
grid.arrange(plot1, plot2, ncol = 2)

###绘制MGMT预测值的条形图（那个文件里面的那样）
mgmt_exist1 <- mgmt_exist[1,]
plot <- ggplot(mgmt_exist1, aes(x = 0.5, y = Estimated)) + # x = 0.5 to center the point in the middle of the plot
  geom_rect(aes(xmin = 0.40, xmax = 0.6, ymin = -0.1, ymax = 1.1), fill = NA, color = "black", size = 1.5) +
  geom_point(size = 3) +  # Black point at the pred value
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.025) +  # Horizontal line for the CI
  geom_segment(aes(x = 0.40, xend = 0.6, y = 0.3582, yend = 0.3582), color = "red", linetype = "solid", size = 1) +  # Red dashed line for the cutoff
  scale_x_continuous(limits = c(0, 1), name = "Score (red line cutoff=0.3582)") +
  scale_y_continuous(name = "MGMT promoter status prediction") +
  theme_minimal() +
  theme(
    axis.title.y = element_text(angle = 0, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

ggplot(mgmt_exist1, aes(x = Estimated, y = 0.5)) +  # 交换x和y的值，以适应横向布局
  geom_rect(aes(ymin = 0.40, ymax = 0.6, xmin = -0.1, xmax = 1.1), fill = NA, color = "black", size = 1.5) +
  geom_point(size = 8) +  # Black point at the pred value
  geom_errorbar(aes(xmin = CI_Lower, xmax = CI_Upper), width = 0.025,size = 1.25) +  # Vertical line for the CI
  geom_segment(aes(y = 0.40, yend = 0.6, x = 0.3582, xend = 0.3582), color = "red", linetype = "solid", size = 1) +  # Red solid line for the cutoff
  scale_y_continuous(limits = c(0, 1), name = "Score (red line cutoff=0.3582)") +
  scale_x_continuous(name = "MGMT promoter status prediction", breaks = seq(0, 1, by = 0.25)) +
  theme_minimal() +
  ggtitle("MGMT promoter status prediction")+
  labs(caption = "Score (red line cutoff=0.3582)")+
  theme(
    axis.title.x = element_text(angle = 0, vjust = 0.5),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold",size = 30),
    plot.caption = element_text(hjust = 0.5, face = "bold", size = 15, vjust = 1)
  )



####################################################计算全局平均的M或者B值
###计算并制作samplesheet
Global_Mvalue <- function(parameter1){mean(
  getM(subset(mSetSqFltNoob,select=colnames(mSetSqFltNoob) %in% parameter1)))
}

Global_Mvalue_allsample <- list()
parameter_name <- colnames(mSetSqFltNoob)

for(i in 1:ncol(mSetSqFltNoob)){Global_Mvalue_allsample[[i]] <- Global_Mvalue(parameter_name[i])
names(Global_Mvalue_allsample)[i] <- colnames(mSetSqFltNoob)[i]}

Global_Mvalue_allsample <- do.call(cbind, Global_Mvalue_allsample)
Global_Mvalue_allsample <- data.frame(Global_Mvalue_allsample)
Global_Mvalue_allsample <- t(Global_Mvalue_allsample)
Global_Mvalue_allsample <- data.frame(Global_Mvalue_allsample)
Global_Mvalue_allsample$sample <- substr(rownames(Global_Mvalue_allsample),2,20)
Global_Mvalue_allsample <- merge(Global_Mvalue_allsample, mgmt_exist[, c("sample", "Sample_name")], by = "sample", all.x = TRUE)

Global_Mvalue_allsample <- Global_Mvalue_allsample %>%
  mutate(`primary_or_recurrence` = if_else(grepl("2", Sample_name), "recurrence", "primary"))
Global_Mvalue_allsample$sample_group <- substr(Global_Mvalue_allsample$Sample_name,1,3)
Global_Mvalue_allsample <- Global_Mvalue_allsample %>%
  rename(mean_methy = Global_Mvalue_allsample)
#########只保留66个有原发和复发对照的样本
Global_Mvalue_allsample_66 <- merge(Global_Mvalue_allsample,mgmt_exist_filtered[,c("Status","Estimated","sample")],by.x="sample",by.y="sample")

####计算原复发全局平均甲基化值及其T检验
primary_methy <- Global_Mvalue_allsample_66$mean_methy[Global_Mvalue_allsample_66$primary_or_recurrence == "primary"]
recurrent_methy <- Global_Mvalue_allsample_66$mean_methy[Global_Mvalue_allsample_66$primary_or_recurrence == "recurrence"]
mean(primary_methy)
mean(recurrent_methy)
test_result <- t.test(primary_methy, recurrent_methy)
p_value <- test_result$p.value

##绘制全局甲基化的对比图
Global_Mvalue_allsample_66_plot <- ggplot(Global_Mvalue_allsample_66, aes(x = primary_or_recurrence, y = mean_methy, group = sample_group)) +
  geom_point(position = position_jitter(width = 0.02), aes(color = sample_group)) +
  geom_line(aes(color = sample_group), position = position_dodge(width = 0.02))+
  geom_boxplot(aes(group = primary_or_recurrence),alpha = 0.5, width = 0.5, outlier.shape = NA)+
  labs(title = "Mean methylation between promary and recurrent tumor",
       x = "Recurrence Status",
       y = "Estimated Value") +
  geom_hline(yintercept=0.3582, linetype="dashed", 
             color = "red", size=0.25) +
  #geom_text(aes(label = sample), vjust = -0.5, hjust = 0.5, size = 2.4, position = position_jitter(width = 0.02)) +
  annotate("text", x = 1.5, y = 0.8, label = sprintf("p = %.3f", p_value), size = 4)+
  geom_segment(aes(x = 1.0, y = 0.75, xend = 2.0, yend = 0.75), color = "black",size = 0.4) +
  geom_segment(aes(x = 1.0, y = 0.75, xend = 1.0, yend = 0.72), arrow = arrow(type = "closed", length = unit(0.015, "npc")), color = "black")+
  geom_segment(aes(x = 2.0, y = 0.75, xend = 2.0, yend = 0.72), arrow = arrow(type = "closed", length = unit(0.015, "npc")), color = "black")+
  theme_minimal()

###计算亚型的归属概率
testfunc <- function(fn) {
  df <- read.csv(fn) %>% 
    dplyr::rename(subtype = 1,
                  score=2) %>% 
  dplyr::arrange(-score) %>% 
  #dplyr::slice_head(n=1) %>% 
    dplyr::filter(subtype == "GBM_RTK1") %>% 
  dplyr::mutate(sample = gsub("^.+/([0-9]+_R[0-9]{2}C[0-9]{2}).+$","\\1",fn) )

  return(df)
}
subgroupsheet <- list.files(pattern = "cal\\.csv$", full.names = TRUE,recursive=T)
tt = pbapply::pblapply(subgroupsheet, testfunc)
ty = dplyr::bind_rows(tt)

predict_MGMT_allsample

RTK1_score <- merge(ty, targets_GSAM[, c("Group", "Sample_name","Array_position")], 
                                by.x = "sample", by.y = "Array_position", all.x = TRUE)

RTK1_score <-RTK1_score %>% 
  mutate(type=substr(RTK1_score$Sample_name,1,3))
  mutate%>% mutate(RTK1=if_else(RTK1_score$score>0.9,"RTK1","non-RTK1"))

###subgroup switch figure(Sankey Diagram)
  library(ggplot2)
  library(ggalluvial)
  
type_counts <- RTK1_score %>%
    group_by(type) %>%
    summarise(count = n(), .groups = 'drop')
# 筛选出那些出现至少两次的 type
RTK1_score_filtered <- RTK1_score %>%
  filter(type %in% type_counts$type[type_counts$count > 1])
# RTK1_score_filtered<-RTK1_score_filtered%>%mutate(RTK1 = paste(RTK1, Sample_name,sep = "_"))
RTK1_score_filtered$RTK1 <-substr(RTK1_score_filtered$RTK1,1,8)
  
primary_RTK1_values <- RTK1_score_filtered %>%
    filter(Group == "primary") %>%
     dplyr::select(RTK1)  %>% rename(primary_RTK1_values=RTK1)
  
recurrence_RTK1_values <- RTK1_score_filtered %>%
  filter(Group == "recurrence") %>%
  dplyr::select(RTK1) %>% rename(recurrence_RTK1_values=RTK1)
Freq <- c(rep(1, 29), 3, 1)
RTK1data <- data.frame(primary_RTK1_values,recurrence_RTK1_values,Freq = c(rep(1, 29), 3, 1))
RTK1data$primary_RTK1_values <- paste("primary",RTK1data$primary_RTK1_values,sep="_")
RTK1data$recurrence_RTK1_values <- paste("recurrence",RTK1data$recurrence_RTK1_values,sep="_")

RTK1data1 <- data.frame(
  primary_RTK1_values = c( "RTK I", "RTK I","non-RTK I","non-RTK I"),
  recurrence_RTK1_values = c("RTK I", "non-RTK I","RTK I","non-RTK I"),
  freq=(c(1,3,0,29))
)


ggplot(data = RTK1data1,
       aes(axis1 = primary_RTK1_values, axis2 = recurrence_RTK1_values, y = freq)) +
  geom_alluvium(aes(fill = recurrence_RTK1_values,width = 1/32)) + # 调整宽度
  geom_stratum(width = 1/32) +  # 为每个状态绘制条带
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("primary_RTK1_values", "recurrence_RTK1_values"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = c("non-RTK I" = "lightgreen", "RTK I" = "lightblue")) +
  theme_void() +  # 使用简洁主题
  theme(axis.text.x = element_blank(),  # 移除x轴的文字
        axis.ticks.x = element_blank())  # 移除x轴的刻度














barplot(avg_mVals_mgmt_JAB, main="MGMT Methylation Levels", beside=T, col=c("blue", "red"), names.arg=c("JAB2", "JAB1"))
legend("topright", legend=c("JAB2", "JAB1"), fill=c("blue", "red"))



data <- data.frame(
  Group = c("JAB2", "JAB1"),
  mVals_mgmt_JAB = avg_mVals_mgmt_JAB
)
ggplot(data, aes(x = Group, y = mVals_mgmt_JAB)) +
  geom_col() +
  ggsci::scale_color_jama()+
  theme_minimal() +
  labs(title = "MGMT Methylation Levels", x = "Group", y = "M Values")

data_JAB12 <- as.data.frame(mVals_mgmt_JAB1_2)
data_JAB12 <- rownames_to_column(data_JAB12, var = "CpG")
data_long <- pivot_longer(data_JAB12, 
                          cols = c("207331540058_R05C01", "207331540058_R06C01"),
                          names_to = "Group",
                          values_to = "M_Values")
data_long$Group <- recode(data_long$Group,
                          `207331540058_R05C01` = "JAB2",
                          `207331540058_R06C01` = "JAB1")



ggplot(data_long, aes(x = Group, y = M_Values)) +
  geom_boxplot() +
  ggsci::scale_color_jama()+
  theme_minimal() +
  labs(title = "MGMT Methylation Levels", x = "Group", y = "M Values")


####创建一个显示 mgmt_data_JAB2 数据集中甲基化和未甲基化探针占比的柱状图
library(ggplot2)
library(data.table)

total_meth <- sum(getMeth(mgmt_data_JAB1))
total_unmeth <- sum(getUnmeth(mgmt_data_JAB1))

total <- total_meth + total_unmeth
meth_percent_JAB1 <- total_meth / total * 100
unmeth_percent_JAB1 <- total_unmeth / total * 100
total_meth_unmeth <- meth_percent+unmeth_percent

plot_data <- data.frame(
  Type = c("Methylated", "Unmethylated"),
  Percentage = c(meth_percent, unmeth_percent)
)

ggplot(plot_data, aes(x = "", y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack",width = 0.5)  +
  labs(title = "Methylation and Unmethylation Percentages(JAB1)",
       x = "",
       y = "Percentage (%)") +
  scale_fill_manual(values = c("Methylated" = "#1976D2", "Unmethylated" = "#BBDEFB")) +
  theme_minimal()





total_meth <- sum(getMeth(mgmt_data_JAB2))
total_unmeth <- sum(getUnmeth(mgmt_data_JAB2))

total <- total_meth + total_unmeth
meth_percent_JAB2 <- total_meth / total * 100
unmeth_percent_JAB2 <- total_unmeth / total * 100
total_meth_unmeth <- meth_percent+unmeth_percent

plot_data <- data.frame(
  Type = c("Methylated", "Unmethylated"),
  Percentage = c(meth_percent, unmeth_percent)
)

ggplot(plot_data, aes(x = "", y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack",width = 0.5)  +
  labs(title = "Methylation and Unmethylation Percentages(JAB2)",
       x = "",
       y = "Percentage (%)") +
  scale_fill_manual(values = c("Methylated" = "#1976D2", "Unmethylated" = "#BBDEFB")) +
  theme_minimal()

####合并代码
data_JAB1 <- data.frame(
  Sample = "JAB1",
  Type = c("Methylated", "Unmethylated"),
  Percentage = c(meth_percent_JAB1, unmeth_percent_JAB1)  # 这些值需要你根据实际情况填写
)

data_JAB2 <- data.frame(
  Sample = "JAB2",
  Type = c("Methylated", "Unmethylated"),
  Percentage = c(meth_percent_JAB2, unmeth_percent_JAB2)  # 这些值需要你根据实际情况填写
)

combined_data <- rbind(data_JAB1, data_JAB2)

data <- data.frame(
  Sample = rep(c("JAB1", "JAB2"), each = 2),
  Type = rep(c("Methylated", "Unmethylated"), 2),
  Percentage = c(meth_percent_JAB1, unmeth_percent_JAB1, meth_percent_JAB2, unmeth_percent_JAB2)  # 使用你的实际百分比
)

data <- data %>%
  group_by(Sample) %>%
  mutate(CumSum = cumsum(Percentage) - 0.5 * Percentage)

data$Type <- factor(data$Type, levels = c("Unmethylated", "Methylated"))

ggplot(data, aes(x = Sample, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  geom_point(aes(y = CumSum, color = Type), size = 4, position = position_dodge(width = 0.5)) +
  geom_line(aes(y = CumSum, group = Type, color = Type), position = position_dodge(width = 0.5), size = 1) +
  scale_fill_manual(values = c("Methylated" = "#1976D2", "Unmethylated" = "#BBDEFB")) +
  scale_color_manual(values = c("Methylated" = "#1976D2", "Unmethylated" = "#BBDEFB")) +
  labs(title = "Methylation and Unmethylation Percentages (JAB1 vs JAB2)",
       x = "Sample",
       y = "Percentage (%)") +
  theme_minimal() +
  theme(legend.position = "top") 

