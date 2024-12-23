library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggvenn)

rm(list = ls())
options(stringsAsFactors = F)

# 
MM_table_yy <- read.table("/home/qzh/Analysis/qd/MM_B6/FPKM_dir/YY_MM_B6_FPKM.txt",sep = "\t",header = T)
MM_table_yy_exp <- MM_table_yy[MM_table_yy$MM_GV > 5 |
                                 MM_table_yy$MM_MII > 5 |
                                 MM_table_yy$MM_PN5 > 5 |
                                 MM_table_yy$MM_E2C > 5 |
                                 MM_table_yy$MM_L2C > 5 |
                                 MM_table_yy$MM_4C > 5 |
                                 MM_table_yy$MM_8C > 5 ,]

MM_table_yy_exp$Stage <- "Expression"
# write.table(MM_table_yy_exp[,c(1,2,10)],file = "/home/qzh/Analysis/qd/MM_B6/MM_Expression_gene.list",
#             sep = "\t",col.names = T,row.names = F,quote = F)
# ZGA gene
MM_table_yy <- read.table("/home/qzh/Analysis/qd/MM_B6/FPKM_dir/YY_MM_B6_FPKM.txt",sep = "\t",header = T)
MM_table_yy <- MM_table_yy[which(MM_table_yy$MM_GV < 1 ),]
MM_table_yy <- MM_table_yy[which(MM_table_yy$MM_MII < 1 ),]

MM_table_yy <- MM_table_yy[!duplicated(MM_table_yy$Gene_name), ]

rownames(MM_table_yy) <- MM_table_yy[,1]
MM_table_yy <- MM_table_yy[,-1]
MM_table_yy <- MM_table_yy[which(rowSums(MM_table_yy) > 0),]
gene_id <- rownames(MM_table_yy)

######### PN5
MM_table_yy_PN5 <- MM_table_yy[MM_table_yy$MM_PN5 > 5 & MM_table_yy$MM_MII < 1
                               & MM_table_yy$MM_GV < 1,]
MM_PN5_ZGA <- rownames(MM_table_yy_PN5)
MM_table_yy_PN5 <- data.frame(t(log2(t(MM_table_yy_PN5))))
MM_table_yy_PN5 <- MM_table_yy_PN5[order(MM_table_yy_PN5$MM_PN5,decreasing = T),]
MM_table_yy_PN5$gene <- rownames(MM_table_yy_PN5)
MM_table_yy_PN5$Stage <- "PN5"

######### E2C
MM_table_yy_E2C <- MM_table_yy[MM_table_yy$MM_E2C > 5 & MM_table_yy$MM_PN5 < 1,]
MM_E2C_ZGA <- rownames(MM_table_yy_E2C)
MM_table_yy_E2C <- data.frame(t(log2(t(MM_table_yy_E2C))))
MM_table_yy_E2C <- MM_table_yy_E2C[order(MM_table_yy_E2C$MM_E2C,decreasing = T),]
MM_table_yy_E2C$gene <- rownames(MM_table_yy_E2C)
MM_table_yy_E2C$Stage <- "E2C"

######### L2C
MM_table_yy_L2C <- MM_table_yy[MM_table_yy$MM_L2C > 5 & MM_table_yy$MM_E2C < 1,]
MM_L2C_ZGA <- rownames(MM_table_yy_L2C)
MM_table_yy_L2C <- data.frame(t(log2(t(MM_table_yy_L2C))))
MM_table_yy_L2C <- MM_table_yy_L2C[order(MM_table_yy_L2C$MM_L2C,decreasing = T),]
MM_table_yy_L2C$gene <- rownames(MM_table_yy_L2C)
MM_table_yy_L2C$Stage <- "L2C"

######### 4C
MM_table_yy_4C <- MM_table_yy[MM_table_yy$MM_4C > 5 & MM_table_yy$MM_L2C < 1,]
MM_4C_ZGA <- rownames(MM_table_yy_4C)
MM_table_yy_4C <- data.frame(t(log2(t(MM_table_yy_4C))))
MM_table_yy_4C <- MM_table_yy_4C[order(MM_table_yy_4C$MM_4C,decreasing = T),]
MM_table_yy_4C$gene <- rownames(MM_table_yy_4C)
MM_table_yy_4C$Stage <- "4C"

######### 8C
MM_table_yy_8C <- MM_table_yy[MM_table_yy$MM_8C > 5 & MM_table_yy$MM_4C < 1,]
MM_8C_ZGA <- rownames(MM_table_yy_8C)
MM_table_yy_8C <- data.frame(t(log2(t(MM_table_yy_8C))))
MM_table_yy_8C <- MM_table_yy_8C[order(MM_table_yy_8C$MM_8C,decreasing = T),]
MM_table_yy_8C$gene <- rownames(MM_table_yy_8C)
MM_table_yy_8C$Stage <- "8C"

Merge_MM_log_table <- data.frame(rbind(
  MM_table_yy_PN5,
  MM_table_yy_E2C,
  MM_table_yy_L2C,
  MM_table_yy_4C,
  MM_table_yy_8C))

# Merge_MM_gene <- Merge_MM_gene[!duplicated(Merge_MM_gene$Gene_name), ]
# write.table(Merge_MM_gene,"/home/qzh/Analysis/qd/MM_B6/ZGA_gene/ZGA_gene_id_list.txt",
#             sep = "\t",quote = F,col.names = T,row.names = F)

Merge_MM_log_table <- Merge_MM_log_table[!duplicated(Merge_MM_log_table$gene), ]
MM_ZGA_gene <- Merge_MM_log_table[,c(8,9)]

Merge_MM_log_table <- read.table("/home/qzh/Analysis/qd/MM_B6/FPKM_dir/YY_MM_B6_FPKM.txt",sep = "\t",header = T)
Merge_MM_log_table <- Merge_MM_log_table[,-c(1,6,10)]
Merge_MM_log_table <- Merge_MM_log_table[!duplicated(Merge_MM_log_table$Gene_name), ]

rownames(Merge_MM_log_table) <- Merge_MM_log_table[,1]
Merge_MM_log_table <- Merge_MM_log_table[,-1]
Merge_MM_log_table <- data.frame(t(log2(t(Merge_MM_log_table))))
Merge_MM_log_table <- Merge_MM_log_table[which(rownames(Merge_MM_log_table) %in% MM_ZGA_gene$gene),]

order_map <- match(MM_ZGA_gene$gene, 
                   rownames(Merge_MM_log_table))
Merge_MM_log_table <- Merge_MM_log_table[order_map,]

Merge_MM_log_table <- Merge_MM_log_table %>%
  mutate_all(~ replace(., is.infinite(.), -6))
Merge_MM_log_table[Merge_MM_log_table>=6] <- 6
Merge_MM_log_table[Merge_MM_log_table<=-6] <- -6
MM_ZGA_log_table_3 <- Merge_MM_log_table

heat_plot <- pheatmap::pheatmap(MM_ZGA_log_table_3,
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