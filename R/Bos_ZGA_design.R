library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggvenn)

rm(list = ls())
options(stringsAsFactors = F)

# Expression
BB_table_yy <- read.table("/home/qzh/Analysis/qd/BB_v1/FPKM_dir/YY_BB_FPKM.txt",
                          sep = "\t",header = T)
BB_table_yy_exp <- BB_table_yy[BB_table_yy$BB_GV > 5 |
                                 BB_table_yy$BB_MII > 5 |
                                 BB_table_yy$BB_4C > 5 |
                                 BB_table_yy$BB_8C > 5 |
                                 BB_table_yy$BB_16C > 5 |
                                 BB_table_yy$BB_Blast > 5 ,]

BB_table_yy_exp$Stage <- "Expression"
# write.table(BB_table_yy_exp[,c(1,2,9)],file = "/home/qzh/Analysis/qd/BB_v1/BB_Expression_gene.list",
#             sep = "\t",col.names = T,row.names = F,quote = F)

############ YY 
BB_table_yy <- read.table("/home/qzh/Analysis/qd/BB_v1/FPKM_dir/YY_BB_FPKM.txt",
                          sep = "\t",header = T)
BB_table_yy <- BB_table_yy[,-2]
BB_table_yy <- BB_table_yy[which(BB_table_yy$BB_GV + BB_table_yy$BB_MII < 3 ),]
BB_table_yy <- BB_table_yy[!duplicated(BB_table_yy$Gene_id), ]

rownames(BB_table_yy) <- BB_table_yy[,1]
BB_table_yy <- BB_table_yy[,-1]
BB_table_yy <- BB_table_yy[which(rowSums(BB_table_yy) > 0),]
gene_id <- rownames(BB_table_yy)

######### 4C
BB_table_yy_4C <- BB_table_yy[BB_table_yy$BB_4C > 5 ,]
BB_4C_ZGA <- rownames(BB_table_yy_4C)
BB_table_yy_4C <- data.frame(t(log2(t(BB_table_yy_4C))))
BB_table_yy_4C <- BB_table_yy_4C[order(BB_table_yy_4C$BB_4C,decreasing = T),]
BB_table_yy_4C$gene <- rownames(BB_table_yy_4C)
BB_table_yy_4C$Stage <- "4C"

######### 8C
BB_table_yy_8C <- BB_table_yy[BB_table_yy$BB_8C > 5 & BB_table_yy$BB_4C < 3,]
BB_8C_ZGA <- rownames(BB_table_yy_8C)
BB_table_yy_8C <- data.frame(t(log2(t(BB_table_yy_8C))))
BB_table_yy_8C <- BB_table_yy_8C[order(BB_table_yy_8C$BB_8C,decreasing = T),]
BB_table_yy_8C$gene <- rownames(BB_table_yy_8C)
BB_table_yy_8C$Stage <- "8C"

######### 16C
BB_table_yy_16C <- BB_table_yy[BB_table_yy$BB_16C > 5 & BB_table_yy$BB_8C < 3,]
BB_16C_ZGA <- rownames(BB_table_yy_16C)
BB_table_yy_16C <- data.frame(t(log2(t(BB_table_yy_16C))))
BB_table_yy_16C <- BB_table_yy_16C[order(BB_table_yy_16C$BB_16C,decreasing = T),]
BB_table_yy_16C$gene <- rownames(BB_table_yy_16C)
BB_table_yy_16C$Stage <- "16C"
# 
######### bla
BB_table_yy_bla <- BB_table_yy[BB_table_yy$BB_Blast >= 5 & BB_table_yy$BB_16C <= 3,]
BB_bla_ZGA <- rownames(BB_table_yy_bla)
BB_table_yy_bla <- data.frame(t(log2(t(BB_table_yy_bla))))
BB_table_yy_bla <- BB_table_yy_bla[order(BB_table_yy_bla$BB_Blast,decreasing = T),]
BB_table_yy_bla$gene <- rownames(BB_table_yy_bla)
BB_table_yy_bla$Stage <- "Blastocyst"

Merge_BB_log_table <- data.frame(rbind(
  BB_table_yy_4C,
  BB_table_yy_8C,
  BB_table_yy_16C,
  BB_table_yy_bla))

# Merge_BB_gene <- Merge_BB_log_table[,c(1,2,10)]
# Merge_BB_gene <- Merge_BB_gene[!duplicated(Merge_BB_gene$Gene_id), ]
# write.table(Merge_BB_gene,"/home/qzh/Analysis/qd/BB_v1/BB_ZGA_gene_id_list.txt",
#             sep = "\t",quote = F,col.names = T,row.names = F)

Merge_BB_log_table <- Merge_BB_log_table[!duplicated(Merge_BB_log_table$gene), ]
BB_ZGA_gene <- Merge_BB_log_table[,c(7,8)]

Merge_BB_log_table <- read.table("/home/qzh/Analysis/qd/BB/FPKM_dir/YY_BB_FPKM.txt",
                                 sep = "\t",header = T)
Merge_BB_log_table <- Merge_BB_log_table[,-2]
Merge_BB_log_table <- Merge_BB_log_table[!duplicated(Merge_BB_log_table$Gene_id), ]

rownames(Merge_BB_log_table) <- Merge_BB_log_table[,1]
Merge_BB_log_table <- Merge_BB_log_table[,-1]
Merge_BB_log_table <- data.frame(t(log2(t(Merge_BB_log_table))))
Merge_BB_log_table <- Merge_BB_log_table[which(rownames(Merge_BB_log_table) %in% BB_ZGA_gene$gene),]

order_map <- match(BB_ZGA_gene$gene, 
                   rownames(Merge_BB_log_table))
Merge_BB_log_table <- Merge_BB_log_table[order_map,]

Merge_BB_log_table <- Merge_BB_log_table %>%
  mutate_all(~ replace(., is.infinite(.), -6))
Merge_BB_log_table[Merge_BB_log_table>=6] <- 6
Merge_BB_log_table[Merge_BB_log_table<=-6] <- -6

heat_plot <- pheatmap::pheatmap(Merge_BB_log_table,
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