###########从注释包中找出所有的基因名字
EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
gene_names <- EPICanno$UCSC_RefGene_Name
unique_gene_names <- unique(gene_names)
unique_gene_names <- unique_gene_names[unique_gene_names != ""]
split_gene_names <- unlist(strsplit(unique_gene_names, ";"))
unique_split_gene_names <- unique(split_gene_names)
unique_split_gene_names <- unique_split_gene_names[unique_split_gene_names != ""]
length(unique_split_gene_names)

#####定义一个function提取每个注释包中基因所在区域的cpg位点
cpg_of_sample <- function(parameter) 
{ mvalue_singlegene<- mVals_64[rownames(mVals) %in% rownames(EPICanno[grepl(parameter, EPICanno$UCSC_RefGene_Name),]), ]
  meanmvalue_singlegene<- colMeans(mvalue_singlegene)
}
###这是两种不同的自定义函数，下面这个可读性好很多
extract_gene_methylation <- function(gene_name) {
  # 提取基因对应的注释数据
  gene_cpgs <- EPICanno[grepl(gene_name, EPICanno$UCSC_RefGene_Name), ]
  # 提取基因对应的CpG位点名称
  gene_cpgs_name <- gene_cpgs[,"Name"]
  # 提取基因对应的甲基化值
  gene_mVals <- mVals_64[rownames(mVals_64) %in% gene_cpgs$Name, ]
  # 返回结果
  list(cpgs_name = gene_cpgs_name, mVals_singlegene = gene_mVals)
}
########制作用于差异分析的表格
discrepancy_gene_datasheet <- data.frame(matrix(nrow = length(unique_split_gene_names), ncol = ncol(mVals_64)))
rownames(discrepancy_gene_datasheet) <- unique_split_gene_names
colnames(discrepancy_gene_datasheet) <- colnames(mVals_64)
discrepancy_gene_datasheet <- tibble::tibble(discrepancy_gene_datasheet)
#######把产生的每个基因在64个样本中的平均M值填入创建的数据框中
test <- cpg_of_sample("YTHDF1")
discrepancy_gene_datasheet["EFGR", ]
#######创建一个循环结构把每个基因所在cpg位点的M平均值算出来并填入准备好的用于差异分析的表格中
library(pbapply)
install.packages("pbapply")

cpg_of_sample <- function(parameter) 
{ mvalue_singlegene<- mVals_64[rownames(mVals_64) %in% rownames(EPICanno[grepl(parameter, EPICanno$UCSC_RefGene_Name),]), ]
  column_medians <- apply(mvalue_singlegene, 2, median)
}
######把这个函数拆开一步一步来做
matching_rows <- EPICanno[grepl("YTHDF1", EPICanno$UCSC_RefGene_Name), ]
matching_cpgs <- rownames(matching_rows)
mvalue_singlegene <- mVals_64[rownames(mVals_64) %in% matching_cpgs, ]
column_medians <- apply(mvalue_singlegene, 2, median)



results <- pblapply(rownames(discrepancy_gene_datasheet1), function(gene_name) {
  meanm_value_singlegene <- cpg_of_sample(gene_name)
  return(meanm_value_singlegene)
})

discrepancy_gene_datasheet["YTHDF1", ] <- column_medians

discrepancy_gene_datasheet1 <- discrepancy_gene_datasheet[c(1:500),]
class(results)


unique_split_gene_names%in%rownames(mVals_64)


discrepancy_gene_datasheet[286,]


EPICanno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
EGFR_cpgs <- EPICanno[grepl("INTS4L1", EPICanno$UCSC_RefGene_Name),]
EGFR_cpgs_name <- EGFR_cpgs[,"Name"]
length(EGFR_cpgs_name)

EGFR_mVals <- mVals_64[rownames(mVals) %in% rownames(EGFR_cpgs), ]
nrow(EGFR_mVals)

dim(EGFR_mVals)
class(EGFR_mVals)
length(EGFR_mVals)

rownames(EGFR_cpgs) %in% rownames(mVals) 



cpg_of_sample <- function(parameter) {
  tryCatch({
    # 提取基因对应的注释数据
    gene_rows <- EPICanno[grepl(parameter, EPICanno$UCSC_RefGene_Name), ]
    # 提取对应的行名
    gene_cpgs_names <- rownames(gene_rows)
    # 提取基因对应的甲基化值
    mvalue_singlegene <- mVals_64[rownames(mVals_64) %in% gene_cpgs_names, ]
    # 确保mvalue_singlegene是一个矩阵或数据框
    if (is.null(dim(mvalue_singlegene)) || length(mvalue_singlegene) == 0) {
      warning("mvalue_singlegene is empty or not a matrix/data frame")
      return(rep(NA, ncol(mVals_64)))
    }
    # 计算列均值
    column_medians <- apply(mvalue_singlegene, 2, median)
    }, error = function(e) {
    # 在出错时返回NA值
    warning(paste("Error processing gene:", parameter, " - ", e$message))
    return(rep(NA, ncol(mVals_64)))
  })
}

# 使用pblapply进行并行处理并显示进度条
library(pbapply)
results <- pblapply(rownames(discrepancy_gene_datasheet), function(gene_name) {
  medial_value_singlegene <- cpg_of_sample(gene_name)
  return(medial_value_singlegene)
})

# 将结果填入discrepancy_gene_datasheet1
for (i in 1:length(results)) {
  if (!is.null(results[[i]])) {
    discrepancy_gene_datasheet[rownames(discrepancy_gene_datasheet)[i], ] <- results[[i]]
  }
}
library(limma)

designMatrix <- model.matrix(~Group+sample_group,data=mgmt_exist_filtered_64)
discrepancy_gene_datasheet_nonNA <- na.omit(discrepancy_gene_datasheet)
fit <- lmFit(discrepancy_gene_datasheet_nonNA, designMatrix)
fit2 <- eBayes(fit)
resultDataFrame <- topTable(fit2, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')

table(data$change)
###############################使用85个样本的数据重新计算差异化的基因与DMRs
unique_split_gene_names##从IlluminaHumanMethylationEPICanno.ilm10b4.hg19注释包中找出的所有基因
########制作用于差异分析的表格
discrepancy_gene_datasheet_85 <- data.frame(matrix(nrow = length(unique_split_gene_names), ncol = ncol(mVals_85)))
rownames(discrepancy_gene_datasheet_85) <- unique_split_gene_names
colnames(discrepancy_gene_datasheet_85) <- colnames(mVals_85)
#####计算medians
cpg_of_sample <- function(parameter) {
  tryCatch({
    # 提取基因对应的注释数据
    gene_rows <- EPICanno[grepl(parameter, EPICanno$UCSC_RefGene_Name), ]
    # 提取对应的行名
    gene_cpgs_names <- rownames(gene_rows)
    # 提取基因对应的甲基化值
    mvalue_singlegene <- mVals_85[rownames(mVals_85) %in% gene_cpgs_names, ]
    # 确保mvalue_singlegene是一个矩阵或数据框
    if (is.null(dim(mvalue_singlegene)) || length(mvalue_singlegene) == 0) {
      warning("mvalue_singlegene is empty or not a matrix/data frame")
      return(rep(NA, ncol(mVals_85)))
    }
    # 计算列均值
    column_medians <- apply(mvalue_singlegene, 2, median)
  }, error = function(e) {
    # 在出错时返回NA值
    warning(paste("Error processing gene:", parameter, " - ", e$message))
    return(rep(NA, ncol(mVals_85)))
  })
}

# 使用pblapply进行并行处理并显示进度条
results <- pblapply(rownames(discrepancy_gene_datasheet_85), function(gene_name) {
  medial_value_singlegene <- cpg_of_sample(gene_name)
  return(medial_value_singlegene)
})
# 将结果填入discrepancy_gene_datasheet1
for (i in 1:length(results)) {
  if (!is.null(results[[i]])) {
    discrepancy_gene_datasheet_85[rownames(discrepancy_gene_datasheet_85)[i], ] <- results[[i]]
  }
}

designMatrix <- model.matrix(~Group+sample_group,data=MGMT_85samples)
discrepancy_gene_datasheet_nonNA <- na.omit(discrepancy_gene_datasheet_85)
fit <- lmFit(discrepancy_gene_datasheet_nonNA, designMatrix)
fit2 <- eBayes(fit)
resultDataFrame <- topTable(fit2, coef = "Grouprecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05,
                                   ifelse(logFC > 0 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')

table(data$change)
discrepancy_gene_85 <- data[c(1:17),]
###################做DMR的分析
designMatrix <- model.matrix(~Group+sample_group,data=MGMT_85samples)
annotation_mSetSqFltNoob_85 <- getAnnotation(mSetSqFltNoob_85)
chr <- as.character(annotation_mSetSqFltNoob_85$chr)
pos <- as.integer(annotation_mSetSqFltNoob_85$pos)

myDMR <- bumphunter(mVals_85,
                    design = designMatrix,
                    coef = 2, # design的第2列
                    cutoff = 0.6,
                    nullMethod = "bootstrap",
                    B=100,
                    type = "M",
                    chr = chr,
                    pos = pos
)

head(myDMR$table)
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
library(annotatr)
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

DMR_85 <- myDMR$table
DMR_85$status <- ifelse(DMR_85$value > 0, "hypermethylated", "hypomethylated")
DMR_85 <- DMR_85[order(DMR_85$chr, DMR_85$start), ]


dmr_data <- DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)

library(circlize)
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[!cytoband$df$V1 %in% c("chrX", "chrY"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

circos.genomicTrackPlotRegion(DMR_85, 
                              ylim = c(-1.5, 1.5),
                              panel.fun = function(region, value, ...) {
                                # 遍历每个区域并单独绘制
                                for (i in 1:nrow(region)) {
                                  circos.points(region[i, , drop = FALSE], 
                                                      value[i, , drop = FALSE],
                                                      col = ifelse(value$status[i] == "hypomethylated", "blue", "red"),
)
                                }
                              })

circos.genomicTrack(DMR_85, numeric.column = 4, ylim = c(-1.5, 1.5),
                    panel.fun = function(region, value, ...) {
                      col_vector <- ifelse(value$status == "hypomethylated", "blue", "red")
                      circos.genomicPoints(region, value, col = col_vector, cex = 1, pch = 16)
                    })


test <- DMR_85[c(1:20),]

circos.genomicInitialize(DMR_85)
str(test)






data = generateRandomBed(nc = 2)
circos.genomicTrack(data, numeric.column = 4, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, ...)
                      circos.genomicPoints(region, value)
                      # 1st column in `value` while 4th column in `data`
                      circos.genomicPoints(region, value, numeric.column = 1)
                    })




#######要最开始就添加空轨道，不然基因的位点画不进去
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
#########这个是可以画出来的
circos.genomicTrack(dmr_data,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, col = ifelse(value$status == "hypomethylated", "blue", "red"), pch = 16, cex = 1, ...)
                    })
##########添加第二个轨道，加入那十九个基因的位置信息。
gene_data <- data.frame(
  chr = c("chr10", "chr5", "chr19", "chr12", "chr12", "chr10", "chr17", "chr10", "chr16", "chr13", "chr16", "chr10", "chr7", "chr1", "chr10", "chr10", "chr5"),
  start = c(119033672, 177437889, 45795713, 116910956, 56020905, 74506500, 28607963, 100188298, 4848380, 32024056, 2148144, 99875571, 6590021, 113876816, 119060011, 73423579, 94618669),
  end = c(119080817, 177447944, 45815308, 117031148, 56038435, 74529324, 28614185, 100229610, 4882401, 32031639, 2154161, 99914092, 6608726, 113887581, 119060138, 73433560, 94697621),
  gene = c("EIF3A", "PRR7-AS1", "RSPH6A", "FBXW8", "IKZF4", "LOC102723439", "SGK494", "CHUK", "UBN1", "FRY-AS1", "RAB26", "DNMBP", "C7orf26", "BCL2L15", "SNORA19", "ZMYND17(MSS51)", "SLF1"),
  dummy = 1  # 添加虚拟数值列
)
#####添加19个基因的位置信息
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
########加入前面计算的差异化cpg位点的信息
######先注释logFC＞1的位点看看在哪些地方

# 完成图形绘制
circos.clear()




############单独画每个染色体的circos图
##chr10
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr10"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)
#######要最开始就添加空轨道，不然基因的位点画不进去
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim

dmr_data <- chr10_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)

chr10_DMR_85 <- DMR_85 %>%
  filter(chr == "chr10")
circos.genomicTrack(dmr_data,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, col = ifelse(value$status == "hypomethylated", "blue", "red"), pch = 16, cex = 1, ...)
                    })
########添加基因注释
gene_data<- data.frame(
  chr = c("chr10",  "chr10", "chr10",  "chr10",  "chr10", "chr10"),
  start = c(119033672, 74506500, 100188298,  99875571, 119060011, 73423579),
  end = c(119080817, 74529324, 100229610,  99914092, 119060138, 73433560),
  gene = c("EIF3A", "LOC102723439", "CHUK", "DNMBP", "SNORA19", "ZMYND17(MSS51)"),
  dummy = 1  # 添加虚拟数值列
)
####第二个轨道
circos.genomicLabels(gene_data_chr10, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr12
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr12"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)
####
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim

chr12_DMR_85 <- DMR_85 %>%
  filter(chr == "chr12")

dmr_data <- chr12_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)

circos.genomicTrack(dmr_data,
                    panel.fun = function(region, value, ...) {
                      circos.genomicPoints(region, value, col = ifelse(value$status == "hypomethylated", "blue", "red"), pch = 16, cex = 1, ...)
                    })
gene_data<- data.frame(
  chr = c("chr12",  "chr12"),
  start = c(116910956, 56020905),
  end = c(117031148, 56038435),
  gene = c("FBXW8", "IKZF4"),
  dummy = 1  # 添加虚拟数值列
)

circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
######chr13
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr13"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr13_DMR_85 <- DMR_85 %>%
  filter(chr == "chr13")
#####
dmr_data <- chr13_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr13"),
  start = c(32024056),
  end = c(32031639),
  gene = c("FRY-AS1"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr16
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr16"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr16_DMR_85 <- DMR_85 %>%
  filter(chr == "chr16")
#####
dmr_data <- chr16_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr16"),
  start = c(4848380,2148144),
  end = c(4882401,2154161),
  gene = c("UBN1","RAB26"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr17
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr17"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr17_DMR_85 <- DMR_85 %>%
  filter(chr == "chr17")
#####
dmr_data <- chr17_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr17"),
  start = c(28607963),
  end = c(28614185),
  gene = c("SGK494"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr19
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr19"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr19_DMR_85 <- DMR_85 %>%
  filter(chr == "chr19")
#####
dmr_data <- chr19_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr19"),
  start = c(45795713),
  end = c(45815308),
  gene = c("RSPH6A"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr7
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr7"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr7_DMR_85 <- DMR_85 %>%
  filter(chr == "chr7")
#####
dmr_data <- chr7_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr7"),
  start = c(6590021),
  end = c(6608726),
  gene = c("C7orf26"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr5
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr5"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr5_DMR_85 <- DMR_85 %>%
  filter(chr == "chr5")
#####
dmr_data <- chr5_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr5"),
  start = c(177437889,94618669),
  end = c(177447944,94697621),
  gene = c("PRR7-AS1","SLF1"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)
##########chr1
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr1"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr1_DMR_85 <- DMR_85 %>%
  filter(chr == "chr1")
#####
dmr_data <- chr1_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- data.frame(
  chr = c("chr1"),
  start = c(113876816),
  end = c(113887581),
  gene = c("BCL2L15"),
  dummy = 1  # 添加虚拟数值列
)
#
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
circos.genomicLabels(gene_data, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(10), 
                     padding = 0.02)

##########把分析出来的差异cpg位点也放进轨道看看
#先用logFC>1的筛出来的67个UP位点
data <- data %>% filter(change=="Up")
# 加载annotatr包
library(annotatr)

# 加载需要的注释文件
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(minfi)
data("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# 将注释数据转换为数据框
annotation_df <- as.data.frame(annotation)
# 查看注释数据的前几行
head(annotation_df)
cpg_ids <- c("cg11571741", "cg19016652", "cg14603358", "cg08339887", "cg07497120", "cg17956145", "cg01990518", 
             "cg11234670", "cg08446900", "cg27581666", "cg13570559", "cg03702871", "cg05631447", "cg18426033", 
             "cg25270886", "ch.14.1488981R", "cg22932685", "cg07723659", "cg02347255", "cg13861904", "cg26641076", 
             "cg14869215", "cg00664609", "cg26681770", "cg27313589", "cg00597107", "cg13207490", "cg06147368", 
             "cg00444915", "cg18405140", "cg04257639", "cg18246521", "cg11099722", "cg26647549", "cg08129092", 
             "cg27454102", "cg00610577", "cg19675778", "cg04628369", "cg07316873", "cg14687145", "cg12490395", 
             "cg04267691", "cg10505257", "cg26996616", "cg03938562", "cg10277044", "cg05531737", "ch.15.90570017R", 
             "cg15886450", "cg16519911", "cg21327194", "cg04379348", "cg11641102", "cg07292807", "cg11070172", 
             "cg04730961", "cg02493798", "cg02592133", "cg01377091", "cg27472156", "cg16298768", "cg09464028", 
             "cg12129908", "cg15964468", "cg20041729", "cg23736307")
cpg_annotation <- annotation_df[annotation_df$Name %in% cpg_ids, c("chr", "pos")]
# 将起始位置和结束位置设置为相同
cpg_annotation$end <- cpg_annotation$pos
names(cpg_annotation) <- c("chr", "start", "end")
cpg_annotation <-cpg_annotation %>% mutate(cpgname=rownames(cpg_annotation))
cpg_annotation <- cpg_annotation %>%
  mutate(dummy = 1)
#再用的筛出来的1615个UP位点，没有logFC>1
data <- data%>%filter(change=="Up")
cpg_ids <- data$gene
#
cpg_annotation <- annotation_df[annotation_df$Name %in% cpg_ids, c("chr", "pos")]
#
cpg_annotation$end <- cpg_annotation$pos
names(cpg_annotation) <- c("chr", "start", "end")
cpg_annotation <-cpg_annotation %>% mutate(cpgname=rownames(cpg_annotation))
cpg_annotation <- cpg_annotation %>%
  mutate(dummy = 1)
##
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[!cytoband$df$V1 %in% c("chrX", "chrY"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)


#
#circos.par("track.height" = 0.2, "cell.padding" = c(0.02, 0.02, 0.02, 0.02), "track.margin" = c(0.01, 0.01))

circos.initializeWithIdeogram(cytoband = cytoband_filtered)
circos.genomicTrackPlotRegion(cpg_annotation, 
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
circos.genomicLabels(cpg_annotation, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.01)
#######画chr1的cpg位点的circos图
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr1"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr1_DMR_85 <- DMR_85 %>%
  filter(chr == "chr1")
#####
dmr_data <- chr1_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- cpg_annotation %>% filter(chr=="chr1")

gene_data1<- data.frame(
  chr = c("chr1"),
  start = c(113876816),
  end = c(113887581),
  gene = c("BCL2L15"),
  name="g",
  dummy = 1  # 添加虚拟数值列
)
#
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
circos.genomicTrackPlotRegion(gene_data1, 
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
circos.genomicLabels(gene_data, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)

circos.genomicLabels(gene_data1, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(1), 
                     padding = 0.02)

#######画chr5的cpg位点的circos图
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr5"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr5_DMR_85 <- DMR_85 %>%
  filter(chr == "chr5")
#####
dmr_data <- chr5_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- cpg_annotation %>% filter(chr=="chr5")

gene_data1<- data.frame(
  chr = c("chr5"),
  start = c(177437889,94618669),
  end = c(177447944,94697621),
  gene = c("PRR7-AS1","SLF1"),
  dummy = 1  # 添加虚拟数值列
)
#
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
circos.genomicTrackPlotRegion(gene_data1, 
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
circos.genomicLabels(gene_data, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)

circos.genomicLabels(gene_data1, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(1), 
                     padding = 0.02)
#######画chr7的cpg位点的circos图
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr7"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr7_DMR_85 <- DMR_85 %>%
  filter(chr == "chr7")
#####
dmr_data <- chr7_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- cpg_annotation %>% filter(chr=="chr7")

gene_data1<- data.frame(
  chr = c("chr7"),
  start = c(6590021),
  end = c(6608726),
  gene = c("C7orf26"),
  dummy = 1  # 添加虚拟数值列
)
#
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
circos.genomicTrackPlotRegion(gene_data1, 
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
circos.genomicLabels(gene_data, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)

circos.genomicLabels(gene_data1, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(1), 
                     padding = 0.02)
#######画chr10的cpg位点的circos图
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr10"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr10_DMR_85 <- DMR_85 %>%
  filter(chr == "chr10")
#####
dmr_data <- chr10_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- cpg_annotation %>% filter(chr=="chr10")

gene_data1<- data.frame(
  chr = c("chr10",  "chr10", "chr10",  "chr10",  "chr10", "chr10"),
  start = c(119033672, 74506500, 100188298,  99875571, 119060011, 73423579),
  end = c(119080817, 74529324, 100229610,  99914092, 119060138, 73433560),
  gene = c("EIF3A", "LOC102723439", "CHUK", "DNMBP", "SNORA19", "ZMYND17(MSS51)"),
  dummy = 1  # 添加虚拟数值列
)
#
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
circos.genomicTrackPlotRegion(gene_data1, 
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
circos.genomicLabels(gene_data, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)

circos.genomicLabels(gene_data1, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(1), 
                     padding = 0.02)
#######画chr12的cpg位点的circos图
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 %in% c("chr12"), ]
cytoband$df <- cytoband_filtered
circos.initializeWithIdeogram(cytoband = cytoband_filtered)

chr12_DMR_85 <- DMR_85 %>%
  filter(chr == "chr12")
#####
dmr_data <- chr12_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- cpg_annotation %>% filter(chr=="chr12")

gene_data1<- data.frame(
  chr = c("chr12",  "chr12"),
  start = c(116910956, 56020905),
  end = c(117031148, 56038435),
  gene = c("FBXW8", "IKZF4"),
  dummy = 1  # 添加虚拟数值列
)
#
circos.genomicTrackPlotRegion(gene_data, 
                              panel.fun = function(region, value, ...) {
                                # 这是一个空的轨道，仅用于初始化
                              },
                              track.height = 0.05, bg.border = NA,
                              ylim = c(0, 1))  # 设置ylim
circos.genomicTrackPlotRegion(gene_data1, 
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
circos.genomicLabels(gene_data, 
                     labels.column = 5, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(5), 
                     padding = 0.02)

circos.genomicLabels(gene_data1, 
                     labels.column = 4, 
                     side = "inside", 
                     cex = 0.6, 
                     col = "black", 
                     connection_height = mm_h(1), 
                     padding = 0.02)
##########再次缩放所有的图，chr
# 读取染色体带信息
cytoband <- read.cytoband(species = "hg19")
cytoband_filtered <- cytoband$df[cytoband$df$V1 == "chr10", ]
# 筛选出感兴趣的染色体区间（例如 30Mb 到 50Mb）
cytoband_region <- cytoband_filtered[cytoband_filtered$V2 >= 72000000 & cytoband_filtered$V3 <= 80000000, ]
# 初始化 Circos 图并绘制
circos.initializeWithIdeogram(cytoband = cytoband_region)
#

chr16_DMR_85 <- DMR_85 %>%
  filter(chr == "chr16")
chr16_DMR_85 <- chr16_DMR_85%>%filter(start>2148144,end<2154161)

#####
dmr_data <- chr10_DMR_85[, c("chr", "start", "end", "value","status")]
dmr_data <- dmr_data[order(dmr_data$chr, dmr_data$start), ]
dmr_data$chr <- as.character(dmr_data$chr)
#
gene_data<- cpg_annotation %>% filter(chr=="chr16")
gene_data<-gene_data %>% filter(start>2148144,end<2154161)

#
gene_data1<- data.frame(
  chr = c("chr19"),
  start = c(45795713),
  end = c(45815308),
  gene = c("RSPH6A"),
  dummy = 1  # 添加虚拟数值列
)
#















library(stringr)
dmr_data <- dmr_data[str_order(dmr_data$chr, numeric = TRUE), ]
table(dmr_data$chr)



discrepancy_gene_85
library(GenomicFeatures)
library(rtracklayer)
library(biomaRt)
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)
refGene <- import.ucsc("hg19", "refGene")

genes_85 <- c("EIF3A", "PRR7-AS1", "RSPH6A", "FBXW8", "IKZF4", "LOC102723439", 
           "SGK494", "CHUK", "UBN1", "FRY-AS1", "RAB26", "DNMBP", 
           "C7orf26", "BCL2L15", "SNORA19", "ZMYND17", "SLF1")
split_genes_85 <- unlist(strsplit(EPICanno$UCSC_RefGene_Name, ";"))

gene_annotations <- getBM(
  attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand"),
  filters = "hgnc_symbol",
  values = genes_85,
  mart = ensembl
)
##########先画一个DMR的热图出来
# 添加一列显著性指标，例如使用 area 或其他你认为合适的指标
library(dplyr)
library(tidyr)
dmr_data$regions <- paste(dmr_data$chr, dmr_data$start, dmr_data$end, sep = "_")
DMR_85$significance <- abs(DMR_85$area)
#######筛选出area的值
# 选择 area 值的95百分位数作为阈值
area_threshold <- quantile(DMR_85$area, 0.95)
# 筛选显著性高的 DMR
significant_DMR <- DMR_85[DMR_85$area >= area_threshold, ]
# 查看筛选结果
print(significant_DMR)
nrow(filter(DMR_85,value > 0))
nrow(filter(DMR_85,value < 0))
# 绘制火山图
test <- myDMR$table
library(ggplot2)
ggplot(test, aes(x = value, y = p.value)) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.6) +
  scale_color_discrete(name = "Chromosome") +
  theme_minimal(base_size = 15) +
  labs(
    title = "Volcano Plot of Differential Methylation Regions (DMRs)",
    x = "Effect Size (value)",
    y = "p.value"
  ) +
  theme(legend.position = "right")


DMR_85
test <- DMR_85[,c(1,4)]
teat <- t(test)
dmr_data
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

matching_rows <- grepl("SGK494", EPICanno$UCSC_RefGene_Name)

# 从EPICanno中筛选出这些行
EPICanno_SGK494 <- EPICanno[matching_rows, ]

# 查看结果
print(EPICanno_SGK494)








