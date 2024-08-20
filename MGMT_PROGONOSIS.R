Global_Mvalue_allsample#mean global methylation from 87 samples
mgmt_exist#MGMT promoter methylation predict value of 87samples extracted from folders



targets_primarysample <- targets_GSAM %>%
  filter(Group == "primary")
mSetSqFltNoob_primary <- mSetSqFltNoob[,targets_primarysample$Array_position]
mVals_MGMT_progonosis <- getM(mSetSqFltNoob_primary)
####修剪出samplesheet
mgmt_progonosis_samplesheet <- mgmt_exist%>%filter(sample %in% c(colnames(mSetSqFltNoob_primary)))
mgmt_progonosis_samplesheet <- merge(mgmt_progonosis_samplesheet,MGMT_progonosis_final[,c("sample_group","survivalDays",
                              "progressionFreeDays","survivalFromSecondSurgeryDays","initialMGMT","status","type","status1")],
                              by.x="sample_group",by.y="sample_group")
mgmt_progonosis_samplesheet <- mgmt_progonosis_samplesheet %>% mutate(status1=if_else(grepl("Deceased",mgmt_progonosis_samplesheet$status),1,0))
############给甲基化与未甲基化的列添加逻辑判断列0，1
mgmt_progonosis_samplesheet <- mgmt_progonosis_samplesheet%>%
                   mutate(status1=if_else(grepl("methylated",mgmt_progonosis_samplesheet$Status),1,0))
#############画MGMT甲基化组与未甲基化组的生存曲线
fit <- survfit(Surv(survivalDays, status1) ~ Status, data = mgmt_progonosis_samplesheet)

mgmt_progonosis_survivalplot <- ggsurvplot(fit,
           data = mgmt_progonosis_samplesheet,  
           #fun = "cumhaz",# "event"绘制累积事件(f(y)=1-y)，# "cumhaz"绘制累积危害函数(f(y)=-log(y));# "pct"绘制生存概率(百分比)。
           linetype = 1, 
           surv.median.line = "hv", # 同时显示垂直和水平参考线 即增加中位生存时间   可选 "none"、"hv"、"h"、"v"
           palette = "aaas" ,
           #add.all = TRUE # 添加总患者生存曲线
           xlab = "Time", 
           ylab = "Survival probability", 
           title = "Survival curve between MGMT methylated and unmethylated ",
           legend = c(0.8,0.75), 
           legend.title = "Methyation status", 
           legend.labs = c("Methylated", "Unmethylated") ,
           break.x.by = 250,
           xlim = c(0, 2500),
           #ylim = c(0, 600),
           break.y.by = .25,
           axes.offset = F, 
           conf.int = TRUE,
           ##conf.int.fill  # 设置置信区间填充的颜色
           #conf.int.style # 设置置信区间的类型，有"ribbon"(默认),"step"两种。
           conf.int.alpha = .3,# 数值，指定置信区间填充颜色的透明度；# 数值在0-1之间，0为完全透明，1为不透明。
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

##############患者的预后与mean global 甲基化的关系
meanmethylation_progonosis <- data.frame()
meanmethylation_progonosis <- merge(mgmt_exist,Global_Mvalue_allsample[,c("mean_methy","sample")],by.x="sample",by.y="sample")
meanmethylation_progonosis <- merge(meanmethylation_progonosis,Global_Mvalue_allsample[,c("mean_methy","sample")],by.x="sample",by.y="sample")
####只保留原发患者的样本，因为不想把mean甲基化在原发和复发的可能改变加进来，虽然我之前已经证明mean methylation在原复发间没差异。
meanmethylation_progonosis <- meanmethylation_progonosis%>%filter(Group=="primary")
meanmethylation_progonosis <- merge(meanmethylation_progonosis,survive_data[,c("studyID","survivalDays","status")],by.x="sample_group",by.y="studyID")
meanmethylation_progonosis <- meanmethylation_progonosis %>% mutate(status1=
                            if_else(grepl("Censored",meanmethylation_progonosis$status),0,1))


meanmethylation_progonosis$mean_methy_group <- ifelse(meanmethylation_progonosis$mean_methy > median(meanmethylation_progonosis$mean_methy), "High", "Low")

fit <- survfit(Surv(survivalDays, status1) ~ mean_methy_group, data = meanmethylation_progonosis)
##我按大于平均值与否分为了平均甲基化值高与低组，画了生存曲线，结果没有任何统计学差异
##图片已经画了，存在这个变量中的，我只是删了代码看起来简洁
meanmethylation_progonosis_plot 
#Cox比例风险回归模型（Cox Proportional Hazards Model）
#mean_methy -0.1860：这是平均甲基化程度的回归系数（β值）。负值表示随着甲基化程度增加，风险（如死亡风险）降低。然而，该系数需要进一步解释其显著性。
#exp(coef) = 0.8303：这是风险比，即每增加一个单位的平均甲基化程度，风险减少的倍数。这里，0.8303 表示平均甲基化程度每增加一个单位，风险（如死亡风险）减少约17%（1 - 0.8303 = 0.1697）。
fit_cox <- coxph(Surv(survivalDays, status1) ~ mean_methy, data = meanmethylation_progonosis)



#############################correlation between mean methylation and level of MGMT promoter methylation
Global_Mvalue_allsample
mgmt_exist
MeanMethy_MGMT_sheet <- merge(Global_Mvalue_allsample,mgmt_exist[,c("sample","Estimated","Status")],by.x="sample",by.y="sample")
qqnorm(MeanMethy_MGMT_sheet$mean_methy, col = "red")
hist(MeanMethy_MGMT_sheet$mean_methy, breaks = 10, probability = TRUE, main = "Histogram with Normal Curve")
lines(density(MeanMethy_MGMT_sheet$mean_methy), col = "blue", lwd = 2)
#Shapiro-Wilk检验是用于检验数据正态性的常用统计检验。p值大于0.05，表示数据符合正态分布。
shapiro_test <- shapiro.test(MeanMethy_MGMT_sheet$mean_methy)
shapiro_test <- shapiro.test(MeanMethy_MGMT_sheet$Estimated)
##########相关性检验，因为不是正态分布，没用皮尔逊相关系数
spearman_result <- cor.test(MeanMethy_MGMT_sheet$mean_methy, MeanMethy_MGMT_sheet$Estimated, method = "spearman")
kendall_result <- cor.test(MeanMethy_MGMT_sheet$mean_methy, MeanMethy_MGMT_sheet$Estimated, method = "kendall")

###找出MGMT预测值下降的那些样本
result <- mgmt_exist_filtered %>%
  group_by(sample_group) %>%
  summarize(
    primary_estimated = Estimated[Group == "primary"],
    recurrence_estimated = Estimated[Group == "recurrence"],
    difference = primary_estimated - recurrence_estimated,
    .groups = 'drop'
  )
MGMT_decreasing_sheet <- result %>% filter(difference>0)
MGMT_decreasing_sheet <-data.frame(MGMT_decreasing_sheet)
####找出mean global methylation在复发前后下降的样本是哪些
MeanMethy_decreasing <- merge(mgmt_exist_filtered,Global_Mvalue_allsample[,c("sample","mean_methy")],by.x="sample",by.y="sample")
MeanMethy_decreasing <- MeanMethy_decreasing[, !(names(MeanMethy_decreasing) %in% c("Estimated", "CI_Lower", "CI_Upper", "Cutoff"))]

result_meanmethy <- MeanMethy_decreasing %>%
  group_by(sample_group) %>%
  summarize(
    primary_MeanMethy = mean_methy[Group == "primary"],
    recurrence_MeanMethy = mean_methy[Group == "recurrence"],
    difference = primary_MeanMethy - recurrence_MeanMethy,
    .groups = 'drop'
  )
MeanMethy_decreasing_sheet <- result %>% filter(difference>0)
MeanMethy_decreasing_sheet <- data.frame(MeanMethy_decreasing_sheet)
####平均甲基化下降的样本有12个，MGMT下降的有22个平均甲基化下降的12个只有7个是平均甲基化下降的
##所以平均甲基化下降的位点与MGMT预测值下降的点位并不是一致的呢
#FALSE  TRUE 5     7 
MeanMethy_decreasing_sheet$sample_group %in% MGMT_decreasing_sheet$sample_group
####使用卡方检验（Chi-squared Test）检验平均甲基化下降是否与MGMT预测值下降有关联
###补齐平均甲基化上升以及MGMT预测值上升的两个数据框
MGMT_increasing_sheet <- result %>% filter(difference<0)
MeanMethy_increasing_sheet <- result %>% filter(difference>0)

table(MeanMethy_decreasing_sheet$sample_group %in% MGMT_increasing_sheet$sample_group)

contingency_table <- matrix(c(7, 7, 5, 5), nrow = 2, 
                            dimnames = list(
                              "Decreased_Avg_Methylation" = c("Yes", "No"),
                              "Decreased_MGMT_Prediction" = c("Yes", "No")))
chi_square_test <- chisq.test(contingency_table)
#####使用相关性检验再次验证平均甲基化值与MGMT预测值有没有关联
#先看符合正态分布不
shapiro_test <- shapiro.test(result$difference)#符合
shapiro_test <- shapiro.test(result_meanmethy$difference)#不符合
###有一组不符合正态，不能用Pearson correlation coefficient，使用Spearman's rank correlation coefficient
correlation_spearman <- cor.test(result$difference, result_meanmethy$difference, method = "spearman")
























