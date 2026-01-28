csvdf <- transform(df, tax = factor(tax, rev(levels(tax))))
ggplot(df,
       aes(x = level, stratum = tax, alluvium = sub,
           y = freq,
           fill = tax, label = tax)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  ggtitle("vaccination survey responses at three points in time")

library(reshape2)
data=read.csv("kegg.csv", sep = ',')
data2 <- melt(data)
df=read.csv("ags_city.csv", sep = ',')
merge <- merge(data2,df,by="sample_name", all = F)



df <- read.csv("ko.txt", sep = "\t")
merge <- merge(data,df,by="ko", all = F)
write.csv(merge, "kegg_rpkm_patheay_l2.csv")
data <- as.matrix(data)

# 提取矩阵对角线上的所有值
diagonal_values <- diag(data)

# 提取对应的行名
row_names <- rownames(data)

# 创建包含对角线值和行名的数据框
diagonal_data <- data.frame(Value = diagonal_values, RowName = row_names)
write.csv(merge, "virus-hots-GTDB.csv")

# 打印结果
print(diagonal_data)
brewer.pal(8, 'Dark2')
[1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
[8] "#666666"

library(ggsci)
library("scales")
pal= pal_npg("nrc")(12)
show_col(pal)


df1=read.csv("city_family_lifestyle.csv")
df2=read.csv("kegg_rpkm_patheay_l2.csv")
df2=aggregate(.~l2,df2,sum)

merge <- merge(df1,df2,by="mag")
write.csv(df2, "kegg_rpkm_pathway_l2_sum.csv")

library(reshape2)     #  首先加载一下reshape2包
df1=read.csv("kegg.csv")
df2=read.csv("ags.csv")

df = CAT_KEGG
DF2 <- melt(df1, na.rm = FALSE)
df2=read.csv("groups.txt", sep = "\t")
merge <- merge(DF2,df2,by="sample_name")
ggplot(merge,aes(CITY,phylum))+
  geom_col(aes(fill=variable), position = "fill")+
  #scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(6))+(values = colorRampPalette(brewer.pal(8, "Dark2"))(11))+
  scale_fill_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(11))+
  theme_bw()
write.csv(merge, "kegg_ags.csv")

df1=read_excel("virus_phylum_family_abundance.xlsx", sheet = "phylum")
df2=aggregate(.~phylum,df1,sum)
write.csv(df2, "virus_phylum_abundance_sum_barplot.csv")

#relative abundance
otu <- read.csv("london.csv", header = TRUE, row.names = 1, sep = ',')
otu[otu > 0] <- 1
write.csv(otu, "london2.csv")


total_rpkm <- colSums(otu)
rel_abundance <- data.frame(t(apply(otu, 1, function(x) x/total_rpkm*100)))
write.csv(rel_abundance, "host_family_rel_abundance_sum_barplot.csv")
