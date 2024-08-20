##修剪mgmt_exist，删除没有配对的原发复发样本信息的样本

type_counts <- mgmt_exist %>%
  group_by(sample_group) %>%
  summarise(count = n(),.groups = 'drop')

mgmt_exist_filtered <- mgmt_exist %>%
  filter(sample_group %in% type_counts$sample_group[type_counts$count > 1])
###画出修剪后66个样本的MGMT预测值原发和复发前后的趋势图
ggplot(mgmt_exist_filtered, aes(x = Group, y = Estimated, group = sample_group, label=Sample_name)) +
  geom_point(position = position_jitter(width = 0.02), aes(color = sample_group)) +
  geom_line(aes(color = sample_group), position = position_dodge(width = 0.02))+
  labs(title = "Estimated Values by Sample and Recurrence Status (mgmt_exist_filtered)",
       x = "Recurrence Status",
       y = "Estimated Value") +
  geom_hline(yintercept=0.3582, linetype="dashed", 
             color = "red", size=0.25) +
  #geom_text(size=2.4) +
  geom_text(aes(label = Sample_name), vjust = -0.5, hjust = 0.5, size = 2.4, position = position_jitter(width = 0.02)) +
  theme_minimal()

# 按 sample_group 和 Group 进行分组，并计算 primary 和 recurrence 的 Estimated 值
result <- mgmt_exist_filtered %>%
  group_by(sample_group) %>%
  summarize(
    primary_estimated = Estimated[Group == "primary"],
    recurrence_estimated = Estimated[Group == "recurrence"]
  ) %>%
  filter(recurrence_estimated < primary_estimated) %>%
  ungroup()

# 找出复发状态下 Estimated 值降低的所有样本行
lower_estimated_samples <- mgmt_exist_filtered %>%
  filter(sample_group %in% result$sample_group)
###计算他们之间的差值
##创建一个新的数据框，计算每组的Estimated值的差值
result <- mgmt_exist_filtered %>%
  group_by(sample_group) %>%
  summarize(
    primary_estimated = Estimated[Group == "primary"],
    recurrence_estimated = Estimated[Group == "recurrence"],
    difference = primary_estimated - recurrence_estimated,
    .groups = 'drop'
  )

result <- data.frame(result)
result_ordered <- result[order(-result$difference), ]
result_ordered_mix <- subset(result_ordered, difference > 0.17)
result_alldecreasing <- result_ordered%>%filter(difference>0)


survive_data <- read.csv(file="~/mnt/neuro-genomic-1-ro/Clinical databases/G-SAM/dbGSAM_PUBLIC_VERSION.csv")

MGMT_progonosis <- merge(result_ordered_mix,survive_data[,
                        c("studyID","survivalDays","progressionFreeDays","survivalFromSecondSurgeryDays","initialMGMT","status")],
                        by.x = "sample_group", by.y = "studyID", all.x = TRUE)

MGMT_progonosis <- merge(result_ordered,survive_data[,
                                                     c("studyID","survivalDays","progressionFreeDays","survivalFromSecondSurgeryDays","initialMGMT","status")],
                         by.x = "sample_group", by.y = "studyID", all.x = TRUE)

MGMT_progonosis <- MGMT_progonosis[order(-MGMT_progonosis$difference), ]
MGMT_progonosis_final <- MGMT_progonosis %>% mutate(type=if_else(MGMT_progonosis$difference>0.17,"1","2"))
MGMT_progonosis_final <- MGMT_progonosis_final %>% mutate(status1=if_else(grepl("Deceased",MGMT_progonosis_final$status),1,0))


fit <- survfit(Surv(survivalDays, status1) ~ type, data = MGMT_progonosis_final)
install.packages("survival")
library(survival)
library(dplyr)
####生存曲线图
install.packages("survminer")
library(survminer)
ggsurvplot(fit,
           data = MGMT_progonosis_alldecreasing,  # 指定变量数据来源
           #fun = "cumhaz",# "event"绘制累积事件(f(y)=1-y)，# "cumhaz"绘制累积危害函数(f(y)=-log(y));# "pct"绘制生存概率(百分比)。
           linetype = 1, # 根据分层更改线型c(0,1) or  c("solid", "dashed") or "strata"
           surv.median.line = "hv", # 同时显示垂直和水平参考线 即增加中位生存时间   可选 "none"、"hv"、"h"、"v"
           palette = "aaas" ,#定义颜色 可选调色板有 "hue" "grey","npg","aaas","lancet","jco", "ucscgb","uchicago","simpsons"和"rickandmorty".
           #add.all = TRUE # 添加总患者生存曲线
           #图标题和坐标轴标签 图例标题和位置
           xlab = "Time", # 指定x轴标签
           ylab = "Survival probability", # 指定y轴标签
           title = "Survival curve of 22 decreasing MGMT sample vs 11",# 指定title标签
           legend = c(0.8,0.75), # 指定图例位置 "top"(默认),"bottom","left","right","none"
           legend.title = "Gender", # 设置图例标题
           legend.labs = c("Male", "Female") ,# 指定图例分组标签
           #坐标轴范围、刻度间距
           break.x.by = 250,# 设置x轴刻度间距
           xlim = c(0, 2500),#设置x轴范围
           #ylim = c(0, 600),#设置x轴范围
           break.y.by = .25,# 设置y轴刻度间距
           axes.offset = F, # 逻辑词，默认为TRUE。为FALSE，则生存曲线图的坐标轴从原点开始。
           #置信区间
           conf.int = TRUE,#增加置信区间，但是这东西并没啥实质作用
           ##conf.int.fill  # 设置置信区间填充的颜色
           #conf.int.style # 设置置信区间的类型，有"ribbon"(默认),"step"两种。
           conf.int.alpha = .3,# 数值，指定置信区间填充颜色的透明度；# 数值在0-1之间，0为完全透明，1为不透明。
           #P值文本大小和位置
           pval = TRUE, #log 秩检验  看两个曲线之间有无显著区别
           pval.size = 3,# 指定p值文本大小的数字，默认为 5。
           pval.coord = c(100,.3),# 长度为2的数字向量，指定p值位置x、y，如pval.coord=c(x,y)。
           #pval.method = T,#展示p统计的检验方法
           #pval.method.size = 3,# 指定检验方法 log.rank 文本的大小
           #pval.method.coord =  c(10,.3),# 指定检验方法 log.rank 文本的坐标
           #删失点
           censor = T, # 逻辑词，默认为TRUE，在图上绘制删失点。
           censor.shape = 3, # 数值或字符，用于指定删失点的形状；默认为"+"(3), 可选"|"(124)。
           censor.size = 4.5,# 指定删失点形状的大小，默认为4.5。
           #生存表
           #risk.table = TRUE, # 添加风险表""absolute" 显示处于风险中的绝对数量；# "percentage" 显示处于风险中的百分比数量# "abs_pct" 显示处于风险中的绝对数量和百分比
           risk.table.col = c(1),#"strata", # 根据分层更改风险表颜色 c(1, 2)  or c("solid", "dashed").
           fontsize = 4,# 指定风险表和累积事件表的字体大小。
           tables.y.text = T,# 逻辑词，默认显示生存表的y轴刻度标签；为FALSE则刻度标签被隐藏
           tables.y.text.col = F, # 逻辑词，默认FALSE；为TRUE，则表的y刻度标签将按strata着色。
           tables.height = .3,# 指定所有生存表的高度，数值在0-1之间，默认为0.25.
           #累积事件表
           #cumevents = T,#累计死亡人数
           #cumevents.height = .3,
           #累积删失表
           #cumcensor =T,#累计删失人数
           #cumcensor.height = .3,
           ggtheme = theme_survminer(), #图的主题
           tables.theme = theme_bw()#下面图的主题
)
####修剪样本表，我想看所有MGMT预测值下降的样本的预后与其他预后有没有关系（上面是四个下降最明显的与预后没有关系）
MGMT_progonosis_alldecreasing <- MGMT_progonosis %>% mutate(type=if_else(MGMT_progonosis$difference>0,"1","2"))
MGMT_progonosis_alldecreasing <- MGMT_progonosis_alldecreasing %>% mutate(status1=if_else(grepl("Deceased",MGMT_progonosis_alldecreasing$status),1,0))
fit <- survfit(Surv(survivalDays, status1) ~ type, data = MGMT_progonosis_alldecreasing)

###########################################volcanic map for genome
library(limma)
designMatrix <- model.matrix(~ 0 + factor(targets_GSAM$Group))
colnames(designMatrix) <- levels(factor(c("primary", "recurrence")))
contrastMatrix <- makeContrasts(primaryVsrecurrence = primary - recurrence, levels = designMatrix)
fit <- lmFit(mVals, designMatrix)
#在进行线性模型拟合时，lmFit函数将根据designMatrix中的信息，
#将mVals中的每个基因表达量数据与其对应的样本条件关联起来。
#例如，对于mVals中的第一个样本（normal.GSM4495491），
#designMatrix中的第一行（cancer为0，normal为1）告诉lmFit这是一个正常样本。
fit2 <- contrasts.fit(fit, contrastMatrix)
#这个对比矩阵设定了需要进行比较的具体参数（在这个例子中是癌症样本与正常样本的表达量差异）。
#contrasts.fit 的输出，即 fit2，包含了调整后的模型，这个模型专注于对比矩阵中定义的比较（cancer - normal），
#这允许研究者直接检验两种条件之间的差异。
fit2 <- eBayes(fit2)
#eBayes 函数通常会输出的一些关键统计结果：系数估计（Coefficients）:这些是对模型中每个变量（如不同治疗组）的效应大小的估计。
#标准误（Standard Errors）:对于每个系数估计的标准误差，表示估计值的不确定性。
#t-统计量（t-statistics）:计算得到的每个系数的 t-统计量，用于测试每个系数是否显著不同于零。
#p-值（P-values）:对数-似然（Log-odds）:对数似然比率，通常用于贝叶斯统计分析中，反映了模型拟合数据的好坏。
#F-统计量（F-statistics）:如果使用了多参数对比，eBayes 可能会计算 F-统计量，用于检验模型中多个参数同时为零的假设。
#调整的 p-值（Adjusted P-values）:当涉及多重比较时，eBayes 可能会包含调整后的 p-值，如使用 Benjamini-Hochberg 方法进行校正，
#以控制假发现率（FDR）等等

resultDataFrame <- topTable(fit2, coef="primaryVsrecurrence", number=Inf)

data <- 
  resultDataFrame %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')
table(data$change)

volcanicplot_for_ <- ggplot(data,aes(logFC, -log10(P.Value)))+
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
  geom_text_repel(data = filter(data, abs(logFC) > 1 & -log10(P.Value) > 5),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(adj.P.Val)") #xlab和ylab函数分别用来设置图形的X轴和Y轴标签

##########用预后再来分一个组呢
















 

