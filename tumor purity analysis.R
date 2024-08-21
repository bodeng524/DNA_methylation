BiocManager::install("InfiniumPurify")
library(InfiniumPurify)
test <- MGMT_85samples$Sample_name
basedir <- "~/DNA Methylation - EPIC arrays"
test1<- list.files("~/DNA Methylation - EPIC arrays/7147027")
test1 <- grep(test,test1)
test1 <- lapply(test, function(x) grep(x, test1, value = TRUE))
######Youri计算的基于CNA的tumor只有79个，还是使用infiniumpurity再次计算85个的
bVals_85 <- getBeta(mSetSqFltNoob_85)
tumor_purity <- getPurity(bVals_85, tumor.type = "GBM")
tumor_purity <- as.data.frame(tumor_purity)
tumor_purity$sample <- rownames(tumor_purity)
#######
high_purity <- subset(tumor_purity, tumor_purity > 0.7)
#######
high_purity <-data.frame(high_purity)
high_purity$sample <- rownames(high_purity)
########
test <- MGMT_85samples %>% filter(sample %in% high_purity$sample)
test <- cbind(test,high_purity$high_purity)
######PCA先分层tumor purity的值
library(FactoMineR)
#这个函数只能从行里面读取分组信息，所以一般情况下需要对bVals矩阵进行行列对调，然后经PCA函数处理后才可以用于这里函数中.

tumor_purity <- tumor_purity %>%
  mutate(extend2 = ifelse(tumor_purity > 0.7, "high", "low"))
tumor_purity <- tumor_purity %>%
  mutate(extend = case_when(
    tumor_purity > 0.8 ~ "high",
    tumor_purity >= 0.6 & tumor_purity <= 0.8 ~ "medium",
    tumor_purity < 0.6 ~ "low"
  ))

library(factoextra)
dat <- t(mVals_85)
#这个函数只能从行里面读取分组信息，所以一般情况下需要对bVals矩阵进行行列对调，然后经PCA函数处理后才可以用于这里函数中.
dat.pca <- PCA(dat, graph = F) 
###加载loading
test <- dat.pca$var$coord
list_of_dfs <- lapply(1:ncol(test), function(i) {
  data.frame(Variable = rownames(test), Value = test[, i], stringsAsFactors = FALSE)
})
# 给列表中的每个数据框命名为原矩阵中的列名
names(list_of_dfs) <- colnames(test)
dim1 <- list_of_dfs[[1]]
dim1 <- dim1[order(-abs(dim1$Value)), ]

indices <- match(test, rownames(dim1))
indices<- data.frame(Index = indices)
indices <- table(indices)
indices <- as.data.frame(indices)
ggplot(indices, aes(x = Index, y = Freq)) +
  geom_bar(stat = "identity") +
  labs(x = "Index Position", y = "Index Value", title = "Bar Plot of Indices") +
  theme_minimal()



fviz_pca_ind(dat.pca,
             geom.ind = "point", 
             col.ind = MGMT_85samples$Group, 
             addEllipses = TRUE, 
             legend.title = "Groups")

###########
tumor_purity$sample <-rownames(tumor_purity) 
tumor_purity$Group <- MGMT_85samples$Group[tumor_purity$sample %in% MGMT_85samples$sample]
#####box plot of tumor purity between primary and recurrent samples
ggplot(tumor_purity, aes(x = factor(Group), y = tumor_purity)) + 
  geom_boxplot(fill = "lightblue", color = "black") +
  geom_jitter(width = 0.2, size = 1.5, color = "red") + 
  labs(title = "Tumor purity for different groups", 
       x = "Group", 
       y = "tumor_purity") +
  theme_minimal()

wilcox_test_result <- wilcox.test(tumor_purity$tumor_purity[tumor_purity$Group %in% "primary"], tumor_purity$tumor_purity[tumor_purity$Group %in% "recurrence"])
########直接使用limma包做差异分析
designMatrix <- model.matrix(~Group+sample_group,data=MGMT_85samples)
fit <- lmFit(test, designMatrix)
fit2 <- eBayes(fit)
resultDataFrame <- topTable(fit2, coef = "Grouprecurrence", number=Inf)




########################RNA-seq data analysis
BiocManager::install("DESeq2")
library(DESeq2)
unzip("G-SAM.rnaseq.expression-vst.322.zip", list = TRUE)
current_dir <- getwd()
rna_directory <- "~/mnt/neuro-genomic-1-ro/gsam/RNA"
unzip("~/mnt/neuro-genomic-1-ro/gsam/RNA/G-SAM.rnaseq.expression-vst.322.zip", exdir = current_dir)
data <- read.csv("~/DNA Methylation - EPIC arrays/G-SAM.expression.322.vst.csv", header = TRUE, row.names = 1)

#####FALSE29  TRUE56,85个样本中的29个并没有在这322个RNA-seq测序样本中，大概是因为质量不高被删除了 
selected_columns <- MGMT_85samples$Sample_name
selected_columns %in% colnames(data)
#构建分组矩阵
primary <- colnames(data)[as.integer(substr(colnames(data),4,4)) == 1] 
primary
recurrence <- colnames(data)[as.integer(substr(colnames(data),4,4)) == 2] 
recurrence

primary_sample <- data[,primary] #原发样本的表达矩阵
recurrence_sample <- data[,recurrence]   #复发样本的表达矩阵

test <- data.frame(colnames(primary_sample),"primary")
test1 <- data.frame(colnames(recurrence_sample),"recurrence")
colnames(test1) <-c("samplenames","group") 

coldata <- rbind(test,test1)


##### test ----
designMatrix <- model.matrix(~group, data = coldata)
fit <- lmFit(data, designMatrix)
fit <- eBayes(fit)
resultDataFrame <- topTable(fit, coef = "grouprecurrence", number=Inf)
#####临界值设为<-1,>1的时候一个显著基因也没有办法找到
data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 1 ,'Up',ifelse(logFC < -1,'Down','No change')),'No change'))) %>% 
  rownames_to_column('gene')

table(data$change)
########
data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')

table(data$change)
#########volcanic plot for RNA expression visualization 
library(ggrepel)
data <- data %>%
  mutate(highlight = ifelse(gene %in% discrepancy_gene_85, "highlight", "normal"))

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
  # 添加标签
  geom_text_repel(data = filter(data, (abs(logFC) > 0 & adj.P.Val < 0.05) | gene %in% discrepancy_gene_85),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = highlight),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)") #xlab和ylab函数分别用来设置图形的X轴和Y轴标签 +
  scale_color_manual(values = c("highlight" = "red", "normal" = "#999999")) # 设置颜色
#########找到这110个下调的基因在哪里，有没有和DNA层面上甲基化升高的基因和区域相互重合
down_expression_gene <- data %>% filter(change=="Down")
up_expression_gene<- data %>% filter(change=="Up")
discrepancy_gene_85 <- c("EIF3A", "PRR7-AS1", "RSPH6A", "FBXW8", "IKZF4", "LOC102723439", 
                         "SGK494", "CHUK", "UBN1", "FRY-AS1", "RAB26", "DNMBP", "C7orf26", 
                         "BCL2L15", "SNORA19", "ZMYND17(MSS51)", "SLF1")
discrepancy_gene_85 %in% rownames(primary_sample)
discrepancy_gene_85 %in% up_expression_gene$gene
#####那17个甲基化上调的基因都不在RNA上调或者下降的列表中，现在看RNA下或者上调的基因在不在那些DMR区域内
library(biomaRt)
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
down_expression_gene <- down_expression_gene$gene
gene_info <- getBM(attributes = c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position', 'gene_biotype'),
                   filters = 'hgnc_symbol',
                   values = down_expression_gene,
                   mart = ensembl)
colnames(gene_info) <- c("gene","chr","start","end")
gene_info <- gene_info[,-5]
gene_info <- gene_info[, c("chr", "start", "end", "gene")]
gene_info$dummy <- 1
gene_info <- gene_info[-c(45:50,68),]
rownames(gene_info) <- NULL
###把手动查询的和用biomart包查询的基因位置数据框结合起来
gene_data1<- data.frame(
  chr = c("chr1","chr5","chr2","chr5","chr17","chr11","chr11","chr5","chr4","chr21","chr10","chr2","0","chr20","chr19","chr10","chr8","chr12","chr19","chr8","chr10","chr15","chr10","chr16","chr8","X","chr11","chr2","0"),
  start = c(33607472,141201641,86042253,141098278,19015313,68840383,117167721,141205273,336610,45553492,72194585,187867947,0,21479899,23945746,38383296,146080260,123636871,51288985,22536526,18820778,55463576,120451873,53344709,67331824,0,111639392,178129087,0),
  end = c(33608831,141205053,86052514,141154754,19015949,68929908,117168044,141205890,336751,45622947,72196312,188392007,0,21483572,24014998,38649026,146094874,123637924,51289417,22541522,18822265,55468937,120453047,53345901,67341212,0,111649030,178130243,0),
  gene = c("AL662907.1", "AC005753.3", "AC105053.1", "AC005753.2", "AC007952.4", "AP003071.5", "AP000892.3", "AC005753.1",
           "AC079140.2", "GATD3A", "AC022532.1", "AC007319.1", "AC005674.2", "AL133325.3", "AC139769.1", "AL117339.4",
           "AF235103.3", "AC073857.1", "AC010325.1", "AC105046.1", "AL450384.1", "AC011912.1", "AL139407.1", "AC079416.3",
           "RRS1-AS1", "AL356235.1", "AP001781.1", "AC079305.1", "AC087241.3"),
  dummy = 1  # 添加虚拟数值列
)

gene_info$chr <- paste("chr", gene_info$chr, sep = "")


gene_data <- rbind(gene_info,gene_data1)
####看一下chr10的情况
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr10"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

#####
dmr_data <- chr10_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#筛选出chr10上的转录差异表达基因
gene_chr10 <- gene_data %>% filter(chr=="chr10")
###甲基化差异表达的基因
gene_data<- data.frame(
  chr = c("chr10",  "chr10", "chr10",  "chr10",  "chr10", "chr10"),
  start = c(119033672, 74506500, 100188298,  99875571, 119060011, 73423579),
  end = c(119080817, 74529324, 100229610,  99914092, 119060138, 73433560),
  gene = c("EIF3A", "LOC102723439", "CHUK", "DNMBP", "SNORA19", "ZMYND17(MSS51)"),
  dummy = 1  # 添加虚拟数值列
)
#
circos.genomicTrackPlotRegion(gene_chr10, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
#
circos.genomicTrack(dmr_data,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, col = ifelse(value$status == "hypomethylated", "blue", "red"), pch = 16, cex = 1, ...)
                    })
#
circos.genomicLabels(gene_chr10, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)
circos.genomicLabels(gene_data, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)
#######选出tumor purity显著差异的28个样本，看看这些样本选出来的差异cpg点和之前85个有没有重合，有较大重合那就说明tumor purity影响了结果
rownames(high_purity) <- NULL
test <- high_purity$sample
test1 <- MGMT_85samples %>% filter(sample %in% test)
test1 <- test1 %>% filter(Group == "recurrence")
test2 <- test1$sample_group
test2 <- test2[test2 != "remaining"]
test3 <- MGMT_85samples %>% filter(sample_group %in% test2)
test1 <- test1[c(1:6),]
test3 <- cbind(test1,test3)


mSetSqFltNoob_28 <- mSetSqFltNoob_85[,c(test3$sample)]
mVals_28 <- getM(mSetSqFltNoob_28)
designMatrix <- model.matrix(~Group+sample_group,data=test3)
fit <- lmFit(mVals_28, designMatrix)
fit <- eBayes(fit)

resultDataFrame <- topTable(fit, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')


table(data$change)
##########使用相关性分析画每一个CPG位点的x,yplot，再用相关性分析看UP或者DOWN的哪些cpg位点是不是和tumor purity有关系
discrepent_cpg <- data
test <- discrepent_cpg %>% filter(change %in% c("Up", "Down"))
test1 <- mVals_85["cg12679899", ]

test2 <- tumor_purity[order(tumor_purity$tumor_purity), ]
test1 <- data.frame(t(test1))
colnames(test1) <- sub("^X", "", colnames(test1))

#将列名作为变量，列的值作为对应的值
test3 <- melt(test1, id.vars = NULL, measure.vars = names(test1))

test3 <- test3[test2$sample, , drop = FALSE]
rownames(test3) <- test3$variable

ggplot(test4, aes(x = variable, y = value)) +
  geom_point() +
  geom_line(group = 1) +
  labs(x = "Sample", y = "M_Value", title = "X-Y Plot of cg12679899") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

test4 <- merge(test3,test2,by.x="variable",by.y="sample")
test4 <- test4[order(test4$tumor_purity),]

cor.test(test4$tumor_purity,test4$value,method = "spearman")


#X,Yplot默认X轴按字母顺序排，所以要转换为因子变量并按照我们的要求排列test4$variable <- factor(test4$variable, levels = test4$variable)
######现在要计算1619个CPG的相关性分析，并将结果汇总

test1 <- mVals_85[test$gene, ]
list <- list()
list <- lapply(1:nrow(test1), function(i) {
  # 提取每一行并转换为数据框，保留行名作为列名
  ##这里提取的因为是数值类型的数据，转换为数据框的时候自动转换为了长格式，如果保持原样，用t函数转换
  single_row_df <- as.data.frame(test1[i, ])
  
  # 可选：保留行名
  colnames(single_row_df) <- rownames(test1)[i]
  
  return(single_row_df)
})

list <- lapply(list, function(df) {
  # 使用 merge 函数，按行名和 test2$sample 合并
  merged_df <- merge(df, test2, by.x = "row.names", by.y = "sample")
  
  # 删除自动生成的行名列（可选）
  #rownames(merged_df) <- merged_df$Row.names
  #merged_df$Row.names <- NULL
  
  return(merged_df)
})

########计算up and down一共1709个cpg的相关性分析------
list1 <- list()
list1 <- pblapply(list, function(df) {
  # 获取第二列的值
  first_column_values <- df[[2]]
  
  # 计算相关性（使用 Spearman 相关性作为示例）
  correlation_result <- cor.test(first_column_values, df$tumor_purity, method = "spearman")
  
  # 返回相关性结果
  list(
    cg_name = names(df)[2],  # 使用第二列的列名作为 cg_name
    correlation = correlation_result$estimate,
    p.value = correlation_result$p.value
  )
})




designMatrix <- model.matrix(~Group+sample_group,data=MGMT_85samples)
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

data <- data %>% filter(change %in% c("Up", "Down"))

for (i in seq_along(list1)) {
  # 获取当前 list1 元素的 cg_name
  cg_name <- list1[[i]]$cg_name
  
  # 查找 data 数据框中 gene 列中与 cg_name 匹配的行
  match_index <- which(data$gene == cg_name)
  
  if (length(match_index) > 0) {
    # 如果找到匹配的 gene，更新 rho 和 p_value 列
    data$rho[match_index] <- list1[[i]]$correlation
    data$p_value[match_index] <- list1[[i]]$p.value
  }
}





data$p <- ifelse(data$p_value < 0.05, "significant", "no")
data$trend <- ifelse(data$rho < 0, "negative", "positive")
table(data$p)


ggplot(data, aes(x=t, y=rho)) +
  geom_point(pch=16,cex=0.1)


test <- data%>%filter(p=="significant")
test <- test$gene
########把1709个按照rho值为连续颜色标签可视化在volcanic plot中
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
  # 添加标签
  geom_text_repel(data = filter(data, (abs(logFC) > 0 & adj.P.Val < 0.05)),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = rho),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)") #xlab和ylab函数分别用来设置图形的X轴和Y轴标签 +
scale_color_manual(values = c("highlight" = "red", "normal" = "#999999")) # 设置颜色

######
design <- model.matrix(~ Group + tumor_purity+sample_group, data = test)
fit <- lmFit(mVals_85, design)
fit <- eBayes(fit)
resultDataFrame <- topTable(fit, coef = "Grouprecurrence", number=Inf)
data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')


table(data$change)




test <- merge(MGMT_85samples,tumor_purity,by.x="sample",by.y="sample")





ggplot(data, aes(logFC, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(color = rho),  # 使用 rho 值作为颜色
             size = 2, 
             alpha = 0.5,
             shape = 19) +
  theme_bw(base_size = 24) + 
  scale_color_gradientn(colors = c("blue", "white", "red")) +  # 强烈对比色：蓝色到红色
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  # 添加标签
  geom_text_repel(data = filter(data, (abs(logFC) > 0 & adj.P.Val < 0.05)),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = rho),  # 保持标签颜色与点颜色一致
                  size = 2) +
  xlab("Log2FC") +
  ylab("-Log10(adj.P.Val)")

ggplot(data, aes(logFC, -log10(adj.P.Val))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "#999999") +
  geom_point(aes(color = rho),  # 使用 rho 值作为颜色
             size = 2, 
             alpha = 0.5,
             shape = 19) +
  theme_bw(base_size = 24) + 
  scale_color_gradientn(colors = c("blue", "yellow", "red")) +  # 蓝色到红色，中间使用黄色作为过渡
  theme(panel.grid = element_blank(),
        legend.position = 'right',
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold")) +
  # 添加标签
  geom_text_repel(data = filter(data, (abs(logFC) > 0 & adj.P.Val < 0.05)),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = rho),  # 保持标签颜色与点颜色一致
                  size = 2) +
  xlab("Log2FC") +
  ylab("-Log10(adj.P.Val)")





111


