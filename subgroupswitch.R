subgroupswitch_function <- function(fn) {
  df <- read.csv(fn) %>% 
    dplyr::rename(subtype = 1,
                  score=2) %>% 
    dplyr::arrange(-score) %>% 
    dplyr::slice_head(n=1) %>% 
    #dplyr::filter(subtype == "GBM_RTK1") %>% 
    dplyr::mutate(sample = gsub("^.+/([0-9]+_R[0-9]{2}C[0-9]{2}).+$","\\1",fn) )
  
  return(df)
}

subgroupsheet <- list.files(pattern = "cal\\.csv$", full.names = TRUE,recursive=T)
subgroupswitch_sheet = pbapply::pblapply(subgroupsheet, subgroupswitch_function)
subgroupswitch_sheet = dplyr::bind_rows(subgroupswitch_sheet)

subgroupswitch_sheet <- merge(subgroupswitch_sheet,mgmt_exist_filtered[,
                       c("sample","Group","Sample_name","sample_group","Status")],by.x="sample",by.y="sample")

###sankey plot for subgroup switch
subgroupswitch_primary_sheet <- subgroupswitch_sheet%>%filter(Group=="primary")
subgroupswitch_recurrence_sheet <- subgroupswitch_sheet%>%filter(Group=="recurrence")
table(subgroupswitch_primary_sheet$subtype)
table(subgroupswitch_recurrence_sheet$subtype)
#######制作适合sankey plot的数据表
result <- subgroupswitch_sheet %>%
  group_by(sample_group) %>%
  arrange(sample_group, Group) %>%  # 按 sample_group、sample 和 Group 排序
  mutate(previous_Group = lag(Group), previous_subtype = lag(subtype)) %>%  # 获取前一个 Group 和前一个 subtype 的值
  select(sample_group, sample, previous_subtype, subtype)  # 选择需要的列

result <- result %>% filter(!is.na(previous_subtype))
result <- data.frame(result)
#count(category1, category2) 根据 category1 和 category2 列对数据进行分组，并计算每个组合的出现次数。
conversion_counts <- result %>%
  count(previous_subtype, subtype) 

subgroupswitch<- data.frame(
  primary_subgroup = c(conversion_counts$previous_subtype),
  recurrence_subgroup = c(conversion_counts$subtype),
  freq=(c(conversion_counts$n))
)


labels <- subgroupswitch %>%
  group_by(primary_subgroup, recurrence_subgroup) %>%
  summarise(freq = sum(freq)) %>%
  ungroup()

########sankey plot
subgroupswitch_plot <- ggplot(data = subgroupswitch,
       aes(axis1 = primary_subgroup, axis2 = recurrence_subgroup, y = freq)) +
  geom_alluvium(aes(fill = primary_subgroup),width = 1/32) + 
  geom_stratum(width = 1/32) +  # 为每个状态绘制条带
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("primary_subgroup", "recurrence_subgroup"),
                   expand = c(0.15, 0.05)) +
  #scale_fill_manual(values = c("non-RTK I" = "lightgreen", "RTK I" = "lightblue")) +
  theme_void() +  # 使用简洁主题
  theme(axis.text.x = element_blank(),  # 移除x轴的文字
        axis.ticks.x = element_blank())  # 移除x轴的刻度










