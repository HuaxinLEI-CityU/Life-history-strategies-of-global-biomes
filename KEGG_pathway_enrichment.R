library("clusterProfiler")
# import KO_ID
KO_list1=read.csv("ko_enriched_in_marine.csv")[[1]]
KO_list2=read.csv("ko_enriched_in_surface.csv")[[1]]
KO_list3=read.csv("ko_enriched_in_soil.csv")[[1]]
KO_list4=read.csv("ko_enriched_in_wwtp.csv")[[1]]

#enrichment
result1=enrichKEGG(gene =KO_list1,
                   organism = "ko",
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.01
)
result2=enrichKEGG(KO_list2,
                   organism = "ko",
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.01
)
result3=enrichKEGG(KO_list3,
                   organism = "ko",
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.01
)
result4=enrichKEGG(KO_list4,
                   organism = "ko",
                   pvalueCutoff = 0.01,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.01
)
#save results
dotplot(result1, showCategory = 20)
dotplot(result2, showCategory = 15)
dotplot(result3, showCategory = 15)
dotplot(result4, showCategory = 15)

write.csv(file="skml_KO_enrichment_marine.csv",data.frame(result1),row.names=F)
write.csv(file="KO_enrichment_surface.csv",data.frame(result2),row.names=F)
write.csv(file="KO_enrichment_soul.csv",data.frame(result3),row.names=F)
write.csv(file="KO_enrichment_wwtp.csv",data.frame(result4),row.names=F)
