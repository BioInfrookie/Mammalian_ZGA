library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggvenn)

rm(list = ls())
options(stringsAsFactors = F)

############ YY 
PP_table_yy <- read.table("/home/qzh/Analysis/qd/PP/ICSI_v2/FPKM_dir/YY_PP_FPKM.txt",sep = "\t",header = T)
PP_table_yy <- PP_table_yy[,-2]

PP_table_yy <- PP_table_yy[which((PP_table_yy$PP_MII+PP_table_yy$PP_GV) < 4 ),]

PP_table_yy <- PP_table_yy[!duplicated(PP_table_yy$Gene_id), ]

rownames(PP_table_yy) <- PP_table_yy[,1]
PP_table_yy <- PP_table_yy[,-1]
PP_table_yy <- PP_table_yy[which(rowSums(PP_table_yy) > 0),]
gene_id <- rownames(PP_table_yy)

######### PN5
PP_table_yy_PN5 <- PP_table_yy[PP_table_yy$PP_PN5 > 5 & PP_table_yy$PP_MII < 3,]
PP_PN5_ZGA <- rownames(PP_table_yy_PN5)
PP_table_yy_PN5 <- data.frame(t(log2(t(PP_table_yy_PN5))))
PP_table_yy_PN5 <- PP_table_yy_PN5[order(PP_table_yy_PN5$PP_PN5,decreasing = T),]
PP_table_yy_PN5$gene <- rownames(PP_table_yy_PN5)
PP_table_yy_PN5$Stage <- "PN5"

######### M2C
PP_table_yy_2C <- PP_table_yy[PP_table_yy$PP_2C > 5 & PP_table_yy$PP_PN5 < 3,]
PP_2C_ZGA <- rownames(PP_table_yy_2C)
PP_table_yy_2C <- data.frame(t(log2(t(PP_table_yy_2C))))
PP_table_yy_2C <- PP_table_yy_2C[order(PP_table_yy_2C$PP_2C,decreasing = T),]
PP_table_yy_2C$gene <- rownames(PP_table_yy_2C)
PP_table_yy_2C$Stage <- "2C"

######### 4C
PP_table_yy_4C <- PP_table_yy[PP_table_yy$PP_4C > 5 & PP_table_yy$PP_2C < 3,]
PP_4C_ZGA <- rownames(PP_table_yy_4C)
PP_table_yy_4C <- data.frame(t(log2(t(PP_table_yy_4C))))
PP_table_yy_4C <- PP_table_yy_4C[order(PP_table_yy_4C$PP_4C,decreasing = T),]
PP_table_yy_4C$gene <- rownames(PP_table_yy_4C)
PP_table_yy_4C$Stage <- "4C"

######### 8C
PP_table_yy_8C <- PP_table_yy[PP_table_yy$PP_8C > 5 & PP_table_yy$PP_4C < 3,]
PP_8C_ZGA <- rownames(PP_table_yy_8C)
PP_table_yy_8C <- data.frame(t(log2(t(PP_table_yy_8C))))
PP_table_yy_8C <- PP_table_yy_8C[order(PP_table_yy_8C$PP_8C,decreasing = T),]
PP_table_yy_8C$gene <- rownames(PP_table_yy_8C)
PP_table_yy_8C$Stage <- "8C"

Merge_PP_log_table <- data.frame(rbind(
  PP_table_yy_PN5,
  PP_table_yy_2C,
  PP_table_yy_4C,
  PP_table_yy_8C))


Merge_PP_log_table <- Merge_PP_log_table[!duplicated(Merge_PP_log_table$gene), ]
# PP_ZGA_gene <- Merge_PP_log_table[,c(1,2,9)]
# write.table(PP_ZGA_gene,file = "/home/qzh/Analysis/qd/PP/ICSI_v2/PP_ZGA_gene.list",
#             sep = "\t",col.names = T,row.names = F,quote = F)

PP_ZGA_gene <- Merge_PP_log_table[,c(7,8)]

Merge_PP_log_table <- read.table("/home/qzh/Analysis/qd/PP/ICSI_v2/FPKM_dir/YY_PP_FPKM.txt",sep = "\t",header = T)
Merge_PP_log_table <- Merge_PP_log_table[,-2]
Merge_PP_log_table <- Merge_PP_log_table[!duplicated(Merge_PP_log_table), ]

rownames(Merge_PP_log_table) <- Merge_PP_log_table[,1]
Merge_PP_log_table <- Merge_PP_log_table[,-1]
Merge_PP_log_table <- data.frame(t(log2(t(Merge_PP_log_table))))
Merge_PP_log_table <- Merge_PP_log_table[which(rownames(Merge_PP_log_table) %in% PP_ZGA_gene$gene),]

order_map <- match(PP_ZGA_gene$gene, 
                   rownames(Merge_PP_log_table))
Merge_PP_log_table <- Merge_PP_log_table[order_map,]

Merge_PP_log_table <- Merge_PP_log_table %>%
  mutate_all(~ replace(., is.infinite(.), -6))
Merge_PP_log_table[Merge_PP_log_table>=6] <- 6
Merge_PP_log_table[Merge_PP_log_table<=-6] <- -6
PP_ZGA_log_table <- Merge_PP_log_table

heat_plot <- pheatmap::pheatmap(PP_ZGA_log_table,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                width = 2.1,
                                height = 2.5,
                                show_rownames = F,
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")
heat_plot