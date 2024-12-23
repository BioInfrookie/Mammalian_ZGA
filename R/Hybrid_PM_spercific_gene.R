library(dplyr)
# Order by PM_M
rm(list = ls())
options(stringsAsFactors = F)

# Mouse_ZGA
MM_ZGA_gene_table <- read.table("/home/qzh/Analysis/qd/MM_B6/MM_ZGA_gene.list",
                                sep = "\t",header = T)
# PP_ZGA
PP_ZGA_gene_table <- read.table("/home/qzh/Analysis/qd/PP/ICSI_v2/PP_ZGA_gene.list",
                                sep = "\t",header = T)

# Single sample
##  PM_P
PM_P_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_P_FPKM.txt",
                         sep = "\t",header = T)
PM_P_table <- PM_P_table[,-2]
PM_P_table <- PM_P_table[!duplicated(PM_P_table$Gene_id), ]

rownames(PM_P_table) <- PM_P_table[,1]
PM_P_table <- PM_P_table[,-1]

PM_P_Sample_name <- lapply(names(PM_P_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_P_table) <- PM_P_Sample_name

PM_P_table <- PM_P_table[which(rownames(PM_P_table) %in% PP_ZGA_gene_table$Gene_id),]

# 1C
PM_P_1C_table <- PM_P_table[,c(1:15)]
count_greater_than_zero <- apply(PM_P_1C_table > 1, 2, sum)
PM_P_1C_table <- PM_P_1C_table[, count_greater_than_zero > 0, drop = FALSE]

# 2C
PM_P_2C_table <- PM_P_table[,c(16:27)]
count_greater_than_zero <- apply(PM_P_2C_table > 1, 2, sum)
counts <- apply(PM_P_2C_table, 2, function(x) sum(x > 1))
sorted_cols <- names(PM_P_2C_table)[order(-counts)]
PM_P_2C_table <- PM_P_2C_table[, sorted_cols]

# 4C
PM_P_4C_table <- PM_P_table[,c(28:36)]
count_greater_than_zero <- apply(PM_P_4C_table > 1, 2, sum)
counts <- apply(PM_P_4C_table, 2, function(x) sum(x > 1))
sorted_cols <- names(PM_P_4C_table)[order(-counts)]
PM_P_4C_table <- PM_P_4C_table[, sorted_cols]

# 8C
PM_P_8C_table <- PM_P_table[,c(37:41)]
count_greater_than_zero <- apply(PM_P_8C_table > 1, 2, sum)
counts <- apply(PM_P_8C_table, 2, function(x) sum(x > 1))
sorted_cols <- names(PM_P_8C_table)[order(-counts)]
PM_P_8C_table <- PM_P_8C_table[, sorted_cols]

PM_P_sort_table <- data.frame(cbind(
  PM_P_1C_table,
  PM_P_2C_table,
  PM_P_4C_table,
  PM_P_8C_table))

order_map <- match(PP_ZGA_gene_table$Gene_id, 
                   rownames(PM_P_sort_table))
PM_P_sort_table <- PM_P_sort_table[order_map,]
PM_P_sort_table <- PM_P_sort_table[which(rowSums(PM_P_sort_table)>0),]

PM_P_sort_table <- data.frame(t(log2(t(PM_P_sort_table))))

# PM_P_sort_table <- PM_P_sort_table %>%
#   mutate_all(~ replace(., is.infinite(.), -6))
PM_P_sort_table[PM_P_sort_table>=6] <- 6
PM_P_sort_table[PM_P_sort_table<=-6] <- -6

heat_plot <- pheatmap::pheatmap(PM_P_sort_table,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                show_colnames = F,
                                width = 2.1,
                                height = 2.5,
                                filename="/home/qzh/Analysis/qd/PM/Fig2A.PM_P_activate_single_embryo_heatmap.pdf",
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                gaps_col =c(15,26,34),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")

# PM_M
PM_M_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_M_FPKM.txt",
                         sep = "\t",header = T)
PM_M_table <- PM_M_table[,-1]
PM_M_table <- PM_M_table[!duplicated(PM_M_table$Gene_name), ]

rownames(PM_M_table) <- PM_M_table[,1]
PM_M_table <- PM_M_table[,-1]

PM_M_Sample_name <- lapply(names(PM_M_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_M_table) <- PM_M_Sample_name

PM_M_table <- PM_M_table[which(rownames(PM_M_table) %in% MM_ZGA_gene_table$gene),]
PM_M_table <- PM_M_table[,which(names(PM_M_table) %in% names(PM_P_sort_table))]

order_map <- match(MM_ZGA_gene_table$gene, 
                   rownames(PM_M_table))
PM_M_table <- PM_M_table[order_map,]

order_map <- match(names(PM_P_sort_table), 
                   names(PM_M_table))
PM_M_table <- PM_M_table[,order_map]
PM_M_table <- PM_M_table[which(rowSums(PM_M_table)>0),]
PM_M_table <- data.frame(t(log2(t(PM_M_table))))

PM_M_table <- PM_M_table %>%
  mutate_all(~ replace(., is.infinite(.), -6))
PM_M_table[PM_M_table>=6] <- 6
PM_M_table[PM_M_table<=-6] <- -6

heat_plot <- pheatmap::pheatmap(PM_M_table,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                show_colnames = F,
                                width = 2.1,
                                height = 2.5,
                                filename = "/home/qzh/Analysis/qd/PM/Fig2A.PM_M_activate_single_embryo_heatmap.pdf",
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                gaps_col =c(15,26,34),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")

# Merge sample
##  PM_M
PM_M_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_M_FPKM.txt",
                         sep = "\t",header = T)
PM_M_table <- PM_M_table[,-1]
PM_M_table <- PM_M_table[!duplicated(PM_M_table$Gene_name), ]

rownames(PM_M_table) <- PM_M_table[,1]
PM_M_table <- PM_M_table[,-1]

PM_M_Sample_name <- lapply(names(PM_M_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_M_table) <- PM_M_Sample_name

PM_M_table <- PM_M_table[which(rownames(PM_M_table) %in% MM_ZGA_gene_table$gene),]
PM_M_table <- PM_M_table[,which(names(PM_M_table) %in% names(PM_P_sort_table))]

# 1C
PM_M_average <- data.frame(row.names(PM_M_table))
PM_M_1C_table <- PM_M_table[,c(1:15)]
PM_M_average$PM_M_PN5 <- apply(PM_M_1C_table,1,mean)

# 2C
PM_M_2C_table <- PM_M_table[,c(16:26)]
PM_M_average$PM_M_2C <- apply(PM_M_2C_table,1,mean)

# 4C
PM_M_4C_table <- PM_M_table[,c(27:34)]
PM_M_average$PM_M_4C <- apply(PM_M_4C_table,1,mean)

# 8C
PM_M_8C_table <- PM_M_table[,c(35:39)]
PM_M_average$PM_M_8C <- apply(PM_M_8C_table,1,mean)

names(PM_M_average) <- c("Gene_name","PM_M_PN5","PM_M_2C","PM_M_4C","PM_M_8C")
# write.table(PM_M_average,file = "/home/qzh/Analysis/qd/PM/PM_M_activate_average_FPKM.txt",
#             sep = "\t",col.names = T,row.names = F,quote = F)

order_map <- match(MM_ZGA_gene_table$gene, 
                   PM_M_average$Gene_name)
PM_M_average <- PM_M_average[order_map,]

PM_M_average <- PM_M_average[,-1]
PM_M_average <- PM_M_average[which(rowSums(PM_M_average)>0),]
PM_M_average <- data.frame(t(log2(t(PM_M_average))))

# PM_M_average <- PM_M_average %>%
#   mutate_all(~ replace(., is.infinite(.), -6))
PM_M_average[PM_M_average>=6] <- 6
PM_M_average[PM_M_average<=-6] <- -6

heat_plot <- pheatmap::pheatmap(PM_M_average,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                show_colnames = T,
                                width = 2.1,
                                height = 2.5,
                                filename = "/home/qzh/Analysis/qd/PM/Fig2A.PM_M_activate_average_heatmap.pdf",
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")

# PM_P
PM_P_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_P_FPKM.txt",
                         sep = "\t",header = T)
PM_P_table <- PM_P_table[,-2]
PM_P_table <- PM_P_table[!duplicated(PM_P_table$Gene_id), ]

rownames(PM_P_table) <- PM_P_table[,1]
PM_P_table <- PM_P_table[,-1]

PM_P_Sample_name <- lapply(names(PM_P_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_P_table) <- PM_P_Sample_name

PM_P_table <- PM_P_table[which(rownames(PM_P_table) %in% PP_ZGA_gene_table$Gene_id),]
PM_P_table <- PM_P_table[,which(names(PM_P_table) %in% names(PM_P_sort_table))]

# 1C
PM_P_average <- data.frame(row.names(PM_P_table))
PM_P_1C_table <- PM_P_table[,c(1:15)]
PM_P_average$PM_P_PN5 <- apply(PM_P_1C_table,1,mean)

# 2C
PM_P_2C_table <- PM_P_table[,c(16:26)]
PM_P_average$PM_M_2C <- apply(PM_P_2C_table,1,mean)

# 4C
PM_P_4C_table <- PM_P_table[,c(27:34)]
PM_P_average$PM_M_4C <- apply(PM_P_4C_table,1,mean)

# 8C
PM_P_8C_table <- PM_P_table[,c(35:39)]
PM_P_average$PM_M_8C <- apply(PM_P_8C_table,1,mean)

names(PM_P_average) <- c("Gene_name","PM_P_PN5","PM_P_2C","PM_P_4C","PM_P_8C")

order_map <- match(PP_ZGA_gene_table$Gene_id, 
                   PM_P_average$Gene_name)
PM_P_average <- PM_P_average[order_map,]

PM_P_average <- PM_P_average[,-1]
PM_P_average <- PM_P_average[which(rowSums(PM_P_average)>0),]
PM_P_average <- data.frame(t(log2(t(PM_P_average))))

PM_P_average <- PM_P_average %>%
  mutate_all(~ replace(., is.infinite(.), -6))
PM_P_average[PM_P_average>=6] <- 6
PM_P_average[PM_P_average<=-6] <- -6

heat_plot <- pheatmap::pheatmap(PM_P_average,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                show_colnames = T,
                                width = 2.1,
                                height = 2.5,
                                filename = "/home/qzh/Analysis/qd/PM/Fig2A.PM_P_activate_average_heatmap.pdf",
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")


# Design ZGA gene 
##  PM_M
PM_M_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_M_FPKM.txt",
                         sep = "\t",header = T)
MM_table_yy <- read.table("/home/qzh/Analysis/qd/MM_B6/FPKM_dir/YY_MM_B6_FPKM.txt",
                          sep = "\t",header = T)

MM_table_yy <- MM_table_yy[,-2]
MM_table_yy <- MM_table_yy[!duplicated(MM_table_yy$Gene_id), ]
PM_M_table <- PM_M_table[,-2]
PM_M_table <- PM_M_table[!duplicated(PM_M_table$Gene_id), ]

rownames(PM_M_table) <- PM_M_table[,1]
PM_M_table <- PM_M_table[,-1]

PM_M_Sample_name <- lapply(names(PM_M_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_M_table) <- PM_M_Sample_name

PM_M_table <- PM_M_table[,which(names(PM_M_table) %in% names(PM_P_sort_table))]

# MII
MM_table_MII <- MM_table_yy[,c(3,4)]

# 1C
PM_M_1C_table <- PM_M_table[,c(1:15)]
PM_M_average <- data.frame(apply(PM_M_1C_table,1,mean))

# 2C
PM_M_2C_table <- PM_M_table[,c(16:26)]
PM_M_average$PM_M_2C <- data.frame(apply(PM_M_2C_table,1,mean))

# 4C
PM_M_4C_table <- PM_M_table[,c(27:34)]
PM_M_average$PM_M_4C <- apply(PM_M_4C_table,1,mean)

# 8C
PM_M_8C_table <- PM_M_table[,c(35:39)]
PM_M_average$PM_M_8C <- apply(PM_M_8C_table,1,mean)

PM_M_average <- data.frame(cbind(
  MM_table_MII,
  PM_M_average))

names(PM_M_average) <- c("MM_GV","MM_MII","PM_M_PN5","PM_M_2C","PM_M_4C","PM_M_8C")
rownames(PM_M_average) <- rownames(PM_M_table)
PM_M_table <- PM_M_average

PM_M_table <- PM_M_table[which(rowSums(PM_M_table) > 0),]
gene_id <- rownames(PM_M_table)

######### PN5
PM_M_table_PN5 <- PM_M_table[PM_M_table$PM_M_PN5 > 5 ,]
MM_PN5_ZGA <- rownames(PM_M_table_PN5)
PM_M_table_PN5 <- data.frame(t(log2(t(PM_M_table_PN5))))
PM_M_table_PN5 <- PM_M_table_PN5[order(PM_M_table_PN5$PM_M_PN5,decreasing = T),]
PM_M_table_PN5$gene <- rownames(PM_M_table_PN5)
PM_M_table_PN5$Stage <- "PN5"

######### 2C
PM_M_table_2C <- PM_M_table[PM_M_table$PM_M_2C > 5 & PM_M_table$PM_M_PN5 < 1,]
MM_2C_ZGA <- rownames(PM_M_table_2C)
PM_M_table_2C <- data.frame(t(log2(t(PM_M_table_2C))))
PM_M_table_2C <- PM_M_table_2C[order(PM_M_table_2C$PM_M_2C,decreasing = T),]
PM_M_table_2C$gene <- rownames(PM_M_table_2C)
PM_M_table_2C$Stage <- "2C"

######### 4C
PM_M_table_4C <- PM_M_table[PM_M_table$PM_M_4C > 5 & PM_M_table$PM_M_2C < 1,]
MM_4C_ZGA <- rownames(PM_M_table_4C)
PM_M_table_4C <- data.frame(t(log2(t(PM_M_table_4C))))
PM_M_table_4C <- PM_M_table_4C[order(PM_M_table_4C$PM_M_4C,decreasing = T),]
PM_M_table_4C$gene <- rownames(PM_M_table_4C)
PM_M_table_4C$Stage <- "4C"

######### 8C
PM_M_table_8C <- PM_M_table[PM_M_table$PM_M_8C > 5 & PM_M_table$PM_M_4C < 1,]
MM_8C_ZGA <- rownames(PM_M_table_8C)
PM_M_table_8C <- data.frame(t(log2(t(PM_M_table_8C))))
PM_M_table_8C <- PM_M_table_8C[order(PM_M_table_8C$PM_M_8C,decreasing = T),]
PM_M_table_8C$gene <- rownames(PM_M_table_8C)
PM_M_table_8C$Stage <- "8C"


Merge_PM_M_log_table <- data.frame(rbind(
  PM_M_table_PN5,
  PM_M_table_2C,
  PM_M_table_4C,
  PM_M_table_8C))

Merge_PM_M_log_table <- Merge_PM_M_log_table[!duplicated(Merge_PM_M_log_table$gene), ]
MM_ZGA_gene <- Merge_PM_M_log_table[,c(7,8)]

# write.table(MM_ZGA_gene,file = "/home/qzh/Analysis/qd/PM/YY_PM_M_ZGA.list",
#             sep = "\t",col.names = T,row.names = F,quote = F)
Merge_PM_M_log_table <- Merge_PM_M_log_table[,-c(7,8)]

Merge_PM_M_log_table <- Merge_PM_M_log_table %>%
  mutate_all(~ replace(., is.infinite(.), -6))
Merge_PM_M_log_table[Merge_PM_M_log_table>=6] <- 6
Merge_PM_M_log_table[Merge_PM_M_log_table<=-6] <- -6
MM_ZGA_log_table_3 <- Merge_PM_M_log_table

heat_plot <- pheatmap::pheatmap(MM_ZGA_log_table_3,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                width = 2.1,
                                height = 2.5,
                                show_rownames = F,
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")

# Abnormal
## PM_M
MM_expression_gene_table <- read.table("/home/qzh/Analysis/qd/MM_B6/MM_Expression_gene.list",
                                       sep = "\t",header = T)

Merge_PM_M_log_table <- data.frame(rbind(
  PM_M_table_PN5,
  PM_M_table_2C,
  PM_M_table_4C,
  PM_M_table_8C))

abnormal_PM_M_ZGA_gene <- Merge_PM_M_log_table[which(!(Merge_PM_M_log_table$gene %in% 
                                                         MM_expression_gene_table$Gene_id)),]

abnormal_PM_M_ZGA_gene <- abnormal_PM_M_ZGA_gene[!duplicated(abnormal_PM_M_ZGA_gene$gene), ]
abnormal_PM_M_gene <- abnormal_PM_M_ZGA_gene[,c(7,8)]
abnormal_PM_M_ZGA_gene <- abnormal_PM_M_ZGA_gene[,-c(7,8)]

abnormal_PM_M_ZGA_gene <- abnormal_PM_M_ZGA_gene %>%
  mutate_all(~ replace(., is.infinite(.), -6))
abnormal_PM_M_ZGA_gene[abnormal_PM_M_ZGA_gene>=6] <- 6
abnormal_PM_M_ZGA_gene[abnormal_PM_M_ZGA_gene<=-6] <- -6

heat_plot <- pheatmap::pheatmap(abnormal_PM_M_ZGA_gene,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                width = 2.1,
                                height = 2.5,
                                border_color="white",
                                filename = "/home/qzh/Analysis/qd/PM/PM_M_abnormal_ZGA_heatmap.pdf",
                                legend_breaks=seq(-6,6,3),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")

# PM_M
PM_M_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_M_FPKM.txt",
                         sep = "\t",header = T)
MM_table_yy <- read.table("/home/qzh/Analysis/qd/MM_B6/FPKM_dir/YY_MM_B6_FPKM.txt",
                          sep = "\t",header = T)

MM_table_yy <- MM_table_yy[,-2]
MM_table_yy <- MM_table_yy[!duplicated(MM_table_yy$Gene_id), ]
PM_M_table <- PM_M_table[,-2]
PM_M_table <- PM_M_table[!duplicated(PM_M_table$Gene_id), ]

rownames(PM_M_table) <- PM_M_table[,1]
PM_M_table <- PM_M_table[,-1]

PM_M_Sample_name <- lapply(names(PM_M_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_M_table) <- PM_M_Sample_name

PM_M_table <- PM_M_table[,which(names(PM_M_table) %in% names(PM_P_sort_table))]
PM_M_table <- PM_M_table[which(rownames(PM_M_table) %in% rownames(abnormal_PM_M_ZGA_gene)),]

order_map <- match(rownames(abnormal_PM_M_ZGA_gene), 
                   rownames(PM_M_table))
PM_M_table <- PM_M_table[order_map,]

order_map <- match(names(PM_P_sort_table), 
                   names(PM_M_table))
PM_M_table <- PM_M_table[,order_map]

PM_M_table <- data.frame(t(log2(t(PM_M_table))))

# MB_B_average <- MB_B_average %>%
#   mutate_all(~ replace(., is.infinite(.), -6))
PM_M_table[PM_M_table>=6] <- 6
PM_M_table[PM_M_table<=-6] <- -6

heat_plot <- pheatmap::pheatmap(PM_M_table,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                show_colnames = F,
                                width = 2.1,
                                height = 2.5,
                                filename = "/home/qzh/Analysis/qd/PM/PM_M_abnormal_single_ZGA_heatmap.pdf",
                                border_color=F,
                                legend_breaks=seq(-6,6,3),
                                gaps_col =c(15,26,34),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")
# ggsave(heat_plot,file="/home/qzh/Analysis/qd/PM/PM_M_abnormal_single_ZGA_heatmap.pdf",
#        width = 7,height = 5)

# PM_P
PM_P_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_P_FPKM.txt",
                         sep = "\t",header = T)
PP_table_yy <- read.table("/home/qzh/Analysis/qd/PP/ICSI_v2/FPKM_dir/YY_PP_FPKM.txt",
                          sep = "\t",header = T)

PP_table_yy <- PP_table_yy[,-2]
PP_table_yy <- PP_table_yy[!duplicated(PP_table_yy$Gene_id), ]
PM_P_table <- PM_P_table[,-2]
PM_P_table <- PM_P_table[!duplicated(PM_P_table$Gene_id), ]

rownames(PM_P_table) <- PM_P_table[,1]
PM_P_table <- PM_P_table[,-1]

PM_P_Sample_name <- lapply(names(PM_P_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_P_table) <- PM_P_Sample_name

PM_P_table <- PM_P_table[,which(names(PM_P_table) %in% names(PM_P_sort_table))]

# MII
PP_table_MII <- PP_table_yy[,c(2,3)]

# 1C
PM_P_1C_table <- PM_P_table[,c(1:15)]
PM_P_average <- data.frame(apply(PM_P_1C_table,1,mean))

# 2C
PM_P_2C_table <- PM_P_table[,c(16:26)]
PM_P_average$PM_P_2C <- data.frame(apply(PM_P_2C_table,1,mean))

# 4C
PM_P_4C_table <- PM_P_table[,c(27:34)]
PM_P_average$PM_P_4C <- apply(PM_P_4C_table,1,mean)

# 8C
PM_P_8C_table <- PM_P_table[,c(35:39)]
PM_P_average$PM_P_8C <- apply(PM_P_8C_table,1,mean)

PM_P_average <- data.frame(cbind(
  PP_table_MII,
  PM_P_average))

names(PM_P_average) <- c("PP_GV","PP_MII","PM_P_PN5","PM_P_2C","PM_P_4C","PM_P_8C")
rownames(PM_P_average) <- rownames(PM_P_table)
PM_P_table <- PM_P_average

PM_P_table <- PM_P_table[which((PM_P_table$PP_GV+PM_P_table$PP_MII) < 4 ),]

PM_P_table <- PM_P_table[which(rowSums(PM_P_table) > 0),]
gene_id <- rownames(PM_P_table)

######### PN5
PM_P_table_PN5 <- PM_P_table[PM_P_table$PM_P_PN5 > 5 & PM_P_table$PP_MII < 3,]
MM_PN5_ZGA <- rownames(PM_P_table_PN5)
PM_P_table_PN5 <- data.frame(t(log2(t(PM_P_table_PN5))))
PM_P_table_PN5 <- PM_P_table_PN5[order(PM_P_table_PN5$PM_P_PN5,decreasing = T),]
PM_P_table_PN5$gene <- rownames(PM_P_table_PN5)
PM_P_table_PN5$Stage <- "PN5"

######### 2C
PM_P_table_2C <- PM_P_table[PM_P_table$PM_P_2C > 5 & PM_P_table$PM_P_PN5 < 3,]
MM_2C_ZGA <- rownames(PM_P_table_2C)
PM_P_table_2C <- data.frame(t(log2(t(PM_P_table_2C))))
PM_P_table_2C <- PM_P_table_2C[order(PM_P_table_2C$PM_P_2C,decreasing = T),]
PM_P_table_2C$gene <- rownames(PM_P_table_2C)
PM_P_table_2C$Stage <- "2C"

######### 4C
PM_P_table_4C <- PM_P_table[PM_P_table$PM_P_4C > 5 & PM_P_table$PM_P_2C < 3,]
MM_4C_ZGA <- rownames(PM_P_table_4C)
PM_P_table_4C <- data.frame(t(log2(t(PM_P_table_4C))))
PM_P_table_4C <- PM_P_table_4C[order(PM_P_table_4C$PM_P_4C,decreasing = T),]
PM_P_table_4C$gene <- rownames(PM_P_table_4C)
PM_P_table_4C$Stage <- "4C"

######### 8C
PM_P_table_8C <- PM_P_table[PM_P_table$PM_P_8C > 5 & PM_P_table$PM_P_4C < 3,]
MM_8C_ZGA <- rownames(PM_P_table_8C)
PM_P_table_8C <- data.frame(t(log2(t(PM_P_table_8C))))
PM_P_table_8C <- PM_P_table_8C[order(PM_P_table_8C$PM_P_8C,decreasing = T),]
PM_P_table_8C$gene <- rownames(PM_P_table_8C)
PM_P_table_8C$Stage <- "8C"


Merge_PM_P_log_table <- data.frame(rbind(
  PM_P_table_PN5,
  PM_P_table_2C,
  PM_P_table_4C,
  PM_P_table_8C))

Merge_PM_P_log_table <- Merge_PM_P_log_table[!duplicated(Merge_PM_P_log_table$gene), ]
PM_ZGA_gene <- Merge_PM_P_log_table[,c(7,8)]

# write.table(PM_ZGA_gene,file = "/home/qzh/Analysis/qd/PM/YY_PM_P_ZGA.list",
#             sep = "\t",col.names = T,row.names = F,quote = F)
Merge_PM_P_log_table <- Merge_PM_P_log_table[,-c(7,8)]

Merge_PM_P_log_table <- Merge_PM_P_log_table %>%
  mutate_all(~ replace(., is.infinite(.), -6))
Merge_PM_P_log_table[Merge_PM_P_log_table>=6] <- 6
Merge_PM_P_log_table[Merge_PM_P_log_table<=-6] <- -6

heat_plot <- pheatmap::pheatmap(Merge_PM_P_log_table,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize_col=13,
                                fontsize_row = 14,
                                show_rownames = F,
                                border_color="white",
                                legend_breaks=seq(-6,6,3),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")
heat_plot

PP_ZGA_gene_table <- read.table("/home/qzh/Analysis/qd/PP/ICSI_v2/PP_ZGA_gene.list",
                                sep = "\t",header = T)
PP_expression_gene_table <- read.table("/home/qzh/Analysis/qd/PP/PP_Expression_gene.list",
                                       sep = "\t",header = T)

Merge_PM_P_log_table <- data.frame(rbind(
  PM_P_table_PN5,
  PM_P_table_2C,
  PM_P_table_4C,
  PM_P_table_8C))

abnormal_PM_P_ZGA_gene <- Merge_PM_P_log_table[which(!(Merge_PM_P_log_table$gene %in% 
                                                         PP_expression_gene_table$Gene_id)),]

abnormal_PM_P_ZGA_gene <- abnormal_PM_P_ZGA_gene[!duplicated(abnormal_PM_P_ZGA_gene$gene), ]
abnormal_PM_P_gene <- abnormal_PM_P_ZGA_gene[,c(7,8)]
# write.table(abnormal_PM_P_gene,file = "/home/qzh/Analysis/qd/PM/YY_PM_P_abnormal_ZGA.list",
#             sep = "\t",col.names = T,row.names = F,quote = F)
abnormal_PM_P_ZGA_gene <- abnormal_PM_P_ZGA_gene[,-c(7,8)]

abnormal_PM_P_ZGA_gene <- abnormal_PM_P_ZGA_gene %>%
  mutate_all(~ replace(., is.infinite(.), -6))
abnormal_PM_P_ZGA_gene[abnormal_PM_P_ZGA_gene>=6] <- 6
abnormal_PM_P_ZGA_gene[abnormal_PM_P_ZGA_gene<=-6] <- -6

heat_plot <- pheatmap::pheatmap(abnormal_PM_P_ZGA_gene,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                border_color=F,
                                width = 2.1,
                                height = 2.5,
                                filename = "/home/qzh/Analysis/qd/PM/PM_P_abnormal_ZGA_heatmap.pdf",
                                legend_breaks=seq(-6,6,3),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")
# ggsave(heat_plot,file="/home/qzh/Analysis/qd/PM/PM_P_abnormal_ZGA_heatmap.pdf",
#        width = 2,height = 5)

# PM_P
PM_P_table <- read.table("/home/qzh/Analysis/qd/PM/FPKM_dir/YY_PM_P_FPKM.txt",
                         sep = "\t",header = T)
PP_table_yy <- read.table("/home/qzh/Analysis/qd/PP/ICSI_v2/FPKM_dir/YY_PP_FPKM.txt",
                          sep = "\t",header = T)

PP_table_yy <- PP_table_yy[,-2]
PP_table_yy <- PP_table_yy[!duplicated(PP_table_yy$Gene_id), ]
PM_P_table <- PM_P_table[,-2]
PM_P_table <- PM_P_table[!duplicated(PM_P_table$Gene_id), ]

rownames(PM_P_table) <- PM_P_table[,1]
PM_P_table <- PM_P_table[,-1]

PM_P_Sample_name <- lapply(names(PM_P_table), 
                           function(col_name) unlist(strsplit(col_name, "_only"))[1])
names(PM_P_table) <- PM_P_Sample_name

PM_P_table <- PM_P_table[,which(names(PM_P_table) %in% names(PM_P_sort_table))]

PM_P_table <- PM_P_table[which(rownames(PM_P_table) %in% rownames(abnormal_PM_P_ZGA_gene)),]

order_map <- match(rownames(abnormal_PM_P_ZGA_gene), 
                   rownames(PM_P_table))
PM_P_table <- PM_P_table[order_map,]

order_map <- match(names(PM_P_sort_table), 
                   names(PM_P_table))
PM_P_table <- PM_P_table[,order_map]

PM_P_table <- data.frame(t(log2(t(PM_P_table))))

# MB_B_average <- MB_B_average %>%
#   mutate_all(~ replace(., is.infinite(.), -6))
PM_P_table[PM_P_table>=6] <- 6
PM_P_table[PM_P_table<=-6] <- -6

heat_plot <- pheatmap::pheatmap(PM_P_table,
                                cluster_rows = F,
                                cluster_cols = F,
                                fontsize = 8,
                                show_rownames = F,
                                show_colnames = F,
                                width = 2.1,
                                height = 2.5,
                                border_color=F,
                                filename = "/home/qzh/Analysis/qd/PM/PM_P_abnormal_single_ZGA_heatmap.pdf",
                                legend_breaks=seq(-6,6,3),
                                gaps_col =c(15,26,34),
                                # breaks=bk,
                                color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                angle_col = "90")
# ggsave(heat_plot,file="/home/qzh/Analysis/qd/PM/PM_P_abnormal_single_ZGA_heatmap.pdf",
#        width = 7,height = 5)
