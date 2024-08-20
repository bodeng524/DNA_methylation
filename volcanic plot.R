library(limma)
mgmt_exist_filtered_64 <- mgmt_exist_filtered[-c(1,2),]

designMatrix <- model.matrix(~Group+sample_group,data=mgmt_exist_filtered_64)
####重设行名只需要把行名设为空值，它会自动重排的
rownames(mgmt_exist_filtered_64) <- NULL

mSetSqFltNoob_64 <- mSetSqFltNoob[,c(mgmt_exist_filtered_64$sample)]
mVals_64 <- getM(mSetSqFltNoob_64)
#colnames(design) = levels(factor(group))
#rownames(design) = colnames(exprSet)
con <- "Grouprecurrence-Groupprimary"
#这里生成的对比是 recurrence-primary，表示计算的是 recurrence 组相对于 primary 组的差异。

contrastMatrix <- makeContrasts(contrasts = c(con), levels = designMatrix)
#makeContrasts这一步和contrasts.fit这一步是为了更加复杂的对比指定，通常不需要
fit <- lmFit(mVals_64, designMatrix)
#fit2 <- contrasts.fit(fit, contrastMatrix)
fit2 <- eBayes(fit)

resultDataFrame <- topTable(fit2, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')


table(data$change)
####找出的差异cpg点位的名称
discrepancy_up_cpg <- data %>% filter(change == "Up") %>% select(gene)
discrepancy_down_cpg <- data %>% filter(change == "Down") %>% select(gene)
######找出这些cpg位点的染色体位置
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
cpg_info <- EPICanno[EPICanno$Name %in% c(as.vector(discrepancy_up_cpg$gene)),]
cpg_info <- EPICanno[EPICanno$Name %in% c(as.vector(discrepancy_down_cpg$gene)),]
########把差异值从大到小从小到大进行排列
data[order(-data$logFC), ]
data[order(data$logFC), ]
#######计算最大差异以及最小差异的cpg岛之间的平均甲基化值有没有差异
Global_Mvalue <- function(parameter1){mean(
  getM(mSetSqFltNoob[parameter1, ]))
}

MeanMethy_down_cpg <- Global_Mvalue(discrepancy_down_cpg$gene)
MeanMethy_up_cpg <- Global_Mvalue(discrepancy_up_cpg$gene)
########计算每个差异化的每个cpg位点的M值，并分析两种差异化cpg位点的值有没有统计学的差异
Global_Mvalue_allsample <- list()
parameter_name <- colnames(mSetSqFltNoob)

for(i in 1:nrow(discrepancy_up_cpg)){Global_Mvalue_allsample[[i]] <- mean(getM(mSetSqFltNoob[discrepancy_up_cpg[i,]]))
}
meanmethy_everydown_cpg <- as.data.frame(Global_Mvalue_allsample)
meanmethy_everyup_cpg <- as.data.frame(Global_Mvalue_allsample)

test_result <- t.test(meanmethy_everydown_cpg, meanmethy_everyup_cpg)
####差异基因的火山图，
volcanicplot_for_ <- ggplot(data,aes(logFC, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 2, 
             alpha = 0.5,
             shape = 19)+
  theme_bw(base_size = 24)+ 
  ggsci::scale_color_jama()+ 
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  # 添加标签：含义是：在ggplot2图形中，仅为那些表达变化超过2倍且具有极高统计显著性（调整后P值的负对数大于5）的数据点添加避免重叠的文本标签。这样的设置通常用于突出显示特别重要或显著的结果。
  geom_text_repel(data = filter(data, abs(logFC) > 1 & adj.P.Val < 0.05),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)") #xlab和ylab函数分别用来设置图形的X轴和Y轴标签
###########热图，根本画不了，火山图已经把差异基因筛选出来了，要热图干屁用
#heatmap_66discrepancy_cpg <- pheatmap(mVals_66, show_colnames = F, show_rownames = F,
#scale = "row",
#cluster_cols = F,
#annotation_col = annotation_col,
#breaks = seq(-3, 3, length.out = 100)) 

########用PCA给66给基因分个层
dat <- t(mVals_64)
#这个函数只能从行里面读取分组信息，所以一般情况下需要对bVals矩阵进行行列对调，然后经PCA函数处理后才可以用于这里函数中.
dat.pca <- PCA(dat, graph = F) 
fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = mgmt_exist_filtered_64$Group, 
             addEllipses = TRUE, 
             legend.title = "Groups")
######用预后再来分一个组画火山图呢（前面的是通过原发和复发配对的，这里通过预后配对，但是样本选哪些呢，是44个原发还是33个？？）

###########寻找EGFR基因区域的cpg位点
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
EGFR_cpgs <- EPICanno[grepl("EGFR", EPICanno$UCSC_RefGene_Name),]
EGFR_cpgs_name <- EGFR_cpgs[,"Name"]
length(EGFR_cpgs_name)

EGFR_mVals <- mVals_64[rownames(mVals) %in% rownames(EGFR_cpgs), ]
nrow(EGFR_mVals)
# Create a data frame for plotting
cpgsite <- c(rownames(EGFR_mVals))
EGFR_mVals_long <- data.frame(cpgsite=cpgsite)
mean_Mvalue <- 
EGFR_mVals_long <- EGFR_mVals_long%>%mutate(mean_Mvalue=rowMeans(EGFR_mVals))
EGFR_mVals_long_ordered <- arrange(EGFR_mVals_long,-mean_Mvalue)

##plot the data
ggplot(EGFR_mVals_long, aes(x = cpgsite, y = mean_Mvalue)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  xlab("CpG Sites in EGFR Gene Region") +
  ylab("M-value")
#############用85个数据算出来的M值再画一次带cpg点位的图
EGFR_mVals_85 <- mVals_85[rownames(mVals_85) %in% rownames(EGFR_cpgs), ]
EGFR_mVals_long <- EGFR_mVals_long%>%mutate(mean_Mvalue=rowMeans(EGFR_mVals_85))
EGFR_mVals_long_ordered <- arrange(EGFR_mVals_long,-mean_Mvalue)
######对EGFR中的这些cpg位点做差异分析
designMatrix <- model.matrix(~Group+sample_group,data=mgmt_exist_filtered_64)
fit <- lmFit(EGFR_mVals, designMatrix)
fit2 <- eBayes(fit)
resultDataFrame <- topTable(fit2, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')

table(data$change)
####volcanic map
ggplot(data,aes(logFC, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 2, 
             alpha = 0.5,
             shape = 19)+
  theme_bw(base_size = 24)+ 
  ggsci::scale_color_jama()+ 
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  # 添加标签：含义是：在ggplot2图形中，仅为那些表达变化超过2倍且具有极高统计显著性（调整后P值的负对数大于5）的数据点添加避免重叠的文本标签。这样的设置通常用于突出显示特别重要或显著的结果。
  geom_text_repel(data = filter(data, abs(logFC) > 1 & adj.P.Val < 0.05),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                  color = change),
                  size = 2) +
  xlab("logFC")+
  ylab("-Log10(adj.P.Val)")
  ###xlim(-2,2)
#########画带cpg位点的图
ggplot(data, aes(x = gene, y = t)) +
  geom_point(alpha = 0.5) +
  theme_bw(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold")
  ) +
  xlab("CpG Sites in EGFR Gene Region") +
  ylab("t-score")
#######用85个样本的数据再来筛选一次差异基因
#######制作样本表
remaining_samples <- mgmt_exist[!mgmt_exist$sample_group %in% mgmt_exist_filtered_64$sample_group,]
rownames(remaining_samples) <- NULL
remaining_samples$sample_group <- "remaining"
remaining_samples <- remaining_samples[-(2:3),]
MGMT_85samples <- rbind(remaining_samples,mgmt_exist_filtered_64)
rownames(MGMT_85samples) <- NULL
##拟合
designMatrix <- model.matrix(~Group+sample_group,data=MGMT_85samples)
mSetSqFltNoob_85 <- mSetSqFltNoob[,c(MGMT_85samples$sample)]
mVals_85 <- getM(mSetSqFltNoob_85)
fit <- lmFit(mVals_85, designMatrix)
#fit2 <- contrasts.fit(fit, contrastMatrix)
fit <- eBayes(fit)

resultDataFrame <- topTable(fit, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')


table(data$change)
######使用85个样本的数据再次进行EGFR的差异基因筛选
designMatrix <- model.matrix(~Group+sample_group,data=MGMT_85samples)
EGFR_mVals_85 <- mVals_85[rownames(mVals_85) %in% rownames(EGFR_cpgs), ]
fit <- lmFit(EGFR_mVals_85, designMatrix)
fit2 <- eBayes(fit)
resultDataFrame <- topTable(fit2, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')

table(data$change)
####volcanic map
ggplot(data,aes(t, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  #geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 2, 
             alpha = 0.5,
             shape = 19)+
  theme_bw(base_size = 24)+ 
  ggsci::scale_color_jama()+ 
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  # 添加标签：含义是：在ggplot2图形中，仅为那些表达变化超过2倍且具有极高统计显著性（调整后P值的负对数大于5）的数据点添加避免重叠的文本标签。这样的设置通常用于突出显示特别重要或显著的结果。
  #geom_text_repel(data = filter(data, abs(logFC) > 0 & adj.P.Val < 0.05),
  #max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
  # aes(label = gene, 
  #    color = change),
  #size = 2) +
  xlab("t")+
  ylab("-Log10(adj.P.Val)")+
  xlim(-5,5)

########找一找64个样本中的DMP
group_list <- mgmt_exist_filtered_64$Group
dmp <- dmpFinder(mVals_64, pheno=group_list, type="categorical")
BiocManager::install("DMRcate")
library(DMRcate)
myAnnotation <- cpg.annotate(object = mVals_64, datatype = "array", what = "M", 
                             analysis.type = "differential", design = designMatrix, 
                             contrasts = F, 
                            arraytype = "450K")
DMRs <- dmrcate(myAnnotation, lambda=1000, C=2)
results.ranges <- extractRanges(DMRs)
results.ranges

contrastMatrix <- makeContrasts(contrasts = c(con), levels = designMatrix)
designMatrix <- model.matrix(~Group+sample_group,data=mgmt_exist_filtered_64)

library(TCGAbiolinks)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(limma)
library(missMethyl)
library(DMRcate)
library(Gviz)
library(ggplot2)
library(RColorBrewer)
library(edgeR)
BiocManager::install("missMethyl")

########################################bumphunter计算DMRs
BiocManager::install("bumphunter")
library(bumphunter)
designMatrix <- model.matrix(~Group+sample_group,data=mgmt_exist_filtered_64)
annotation_mSetSqFltNoob_64 <- getAnnotation(mSetSqFltNoob_64)
chr <- as.character(annotation_mSetSqFltNoob_64$chr)
pos <- as.integer(annotation_mSetSqFltNoob_64$pos)

myDMR <- bumphunter(mVals_64,
                    design = designMatrix,
                    coef = 2, # design的第2列
                    cutoff = 0.2,
                    nullMethod = "bootstrap",
                    type = "M",
                    chr = chr,
                    pos = pos
)

head(myDMR$table)
######使用DMRcate函数计算
####安装这个包出现了failed to lock directory错误，解决方法看codes文档2024-07-11所记录
BiocManager::install("DMRcate", dependencies = TRUE, INSTALL_opts = '--no-lock')







library(GenomicRanges)
library(rtracklayer)
BiocManager::install("annotatr")
library(annotatr)
# 从 bumphunter 结果中提取信息
dmr_table <- myDMR$table
# 创建 GRanges 对象
dmr_granges <- GRanges(
  seqnames = dmr_table$chr,
  ranges = IRanges(start = dmr_table$start, end = dmr_table$end),
  strand = "*",
  value = dmr_table$value
)
# 查看 GRanges 对象
BiocManager::install("genomation")
library(genomation)
annotations <- build_annotations(genome = 'hg19', annotations = 'hg19_basicgenes')

dmr_annotated <- annotate_regions(
  regions = dmr_granges,
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE
)
dmr_genes <- as.data.frame(dmr_annotated)
head(dmr_genes)
unique_genes <- unique(dmr_genes$annot.symbol)
#####Circos plot
install.packages("circlize")
library(circlize)
DMR_64 <- myDMR$table
DMR_64$status <- ifelse(DMR_64$value > 0, "hypermethylated", "hypomethylated")
DMR_64 <- DMR_64[order(DMR_64$chr, DMR_64$start), ]

circos.initializeWithIdeogram(species = "hg19")
circos.genomicTrackPlotRegion(dmr_data, 
                              ylim = c(-1.5, 1.5),
                              panel.fun = function(region, value, ...) {
                                # 遍历每个区域并单独绘制
                                for (i in 1:nrow(region)) {
                                  circos.genomicLines(region[i, , drop = FALSE], 
                                                      value[i, , drop = FALSE],
                                                      col = ifelse(value$status[i] == "hypomethylated", "blue", "red"),
                                                      lwd = 2,
                                                      ...)
                                }
                              })


circos.genomicTrackPlotRegion(dmr_data,
                              ylim = c(-1.5, 1.5),
                              panel.fun = function(region, value, ...) {
                                circos.genomicLines(region, value,
                                                   col = ifelse(value$status == "hypomethylated", "blue", "red"),
                                                   lwd = 2,
                                                   ...)
                              })

# 添加图例
legend("topright", legend = c("Hypomethylated", "Hypermethylated"), 
       fill = c("blue", "red"), border = NA)


dmr_data <- DMR_64[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)




















