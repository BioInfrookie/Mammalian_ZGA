## PCA
library(corrplot)
library(gmodels)
library(ggpubr)
library(ggrepel)
library(ggdendro)
library(fdrtool)

rm(list = ls())
options(stringsAsFactors = F)

MM_ZGA_gene_table <- read.table("/home/qzh/Analysis/qd/MM_B6/ZGA_gene/ZGA_gene_id_list_bla.txt",
                                sep = "\t",header = T)

# DEseq
library(DESeq2)
Eif1ad_ASO_table  <- read.table("/home/qzh/Analysis/qd/MM_Eif1ad_KO/02-result/ASO3_Eif1ad_family_count.txt",
                                sep = "\t",header = T)
Xiewei_ZGA_table <- read.table("/home/qzh/Analysis/qd/MM_Eif1ad_KO/02-result/zga_xie_obox_2023.txt",
                               sep = " ",header = T)
ID_trans <- read.table("/home/qzh/ref/mm10/Ensemble_id_trans.txt",header = T,sep = "\\",quote = "")

exprSet <- Eif1ad_ASO_table
rownames(exprSet) <- exprSet[,1]
exprSet <- exprSet[,-1]
exprSet <- exprSet[which(rowSums(exprSet) > 0),]
exprSet <- round(exprSet,digits = 0)
condition <- factor(c(rep("WT",3),rep("Mut",3)),levels = c("WT","Mut"))
colData <- data.frame(sampleName = colnames(exprSet),
                      fileName = colnames(exprSet),
                      condition = condition)

dds <- DESeqDataSetFromMatrix(exprSet,colData,design= ~condition)
dds <- DESeq(dds)
resLFC <- lfcShrink(dds, coef="condition_Mut_vs_WT", type="apeglm")
DEG <- as.data.frame(resLFC)
DEG$gene <- row.names(DEG)
DEG <- na.omit(DEG)
DEG$change = as.factor(ifelse(abs(DEG$log2FoldChange) >=1 & DEG$pvalue <= 0.05,
                              ifelse(DEG$log2FoldChange >= 1,'UP','DOWN'),'NOT'))
this_tile <- paste0('Down gene numbers: ',nrow(DEG[DEG$change == 'DOWN',]),
                    '\nUp gene numbers: ',nrow(DEG[DEG$change == 'UP',]))

MA <- ggplot(DEG, aes(x = log2FoldChange,y = -log10(pvalue), color = change)) +
  geom_point(data=DEG[which(DEG$change=="NOT"),], 
             aes(x = log2FoldChange,y = -log10(pvalue)),size=0.5,fill="grey") +
  geom_point(data=DEG[which(DEG$change=="UP"),], 
             aes(x = log2FoldChange,y = -log10(pvalue)),size=1,fill="red") +
  geom_point(data=DEG[which(DEG$change=="DOWN"),], 
             aes(x = log2FoldChange,y = -log10(pvalue)),size=1,fill="darkgreen") +
  theme_classic() +
  ggtitle(this_tile) + 
  theme(legend.position="none") +
  geom_vline(aes(xintercept=1),linetype="dashed") +
  geom_vline(aes(xintercept=-1),linetype="dashed") +
  xlab("log2 Foldchange") + ylab("-log10(P-value)") +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_color_manual(values=c(NOT = "grey", UP = "red",
                              DOWN="darkgreen")) + 
  theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    axis.line = element_line(linewidth =0),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA))
MA
# ggsave(MA,filename="/home/qzh/Analysis/qd/MM_Eif1ad_KO/Mouse_Eif1ad_KD_heatmap_Volcano.pdf",
#        width = 3.2,height = 3.5)
DEG_table <- DEG
DEG_table <- DEG_table[order(DEG_table$log2FoldChange,decreasing = F),]

# write.table(DEG_table,file = "/home/qzh/Analysis/qd/MM_Eif1ad_KO/02-result/count/ASO2_vs_SCR_L2C_diff_Single_embryo_v3.txt",
#             sep = "\t",col.names = T,row.names = F,quote = F)

DEG_table$ZGA_diff <- as.factor(ifelse(DEG_table$gene %in% Xiewei_ZGA_table$gene_id & DEG_table$change != "NOT",
                                       ifelse(DEG_table$change == "UP",'UP','DOWN'),'NOT'))

this_tile <- paste0('Down ZGA gene numbers: ',nrow(DEG_table[DEG_table$ZGA_diff == 'DOWN',]),
                    '\nUp ZGA gene numbers: ',nrow(DEG_table[DEG_table$ZGA_diff == 'UP',]))

MA <- ggplot(DEG_table, aes(x = log2FoldChange,y = -log10(pvalue), color = ZGA_diff)) +
  geom_point(data=DEG_table[which(DEG_table$ZGA_diff=="NOT"),], 
             aes(x = log2FoldChange,y = -log10(pvalue)),size=0.5,fill="grey") +
  geom_point(data=DEG_table[which(DEG_table$ZGA_diff=="UP"),], 
             aes(x = log2FoldChange,y = -log10(pvalue)),size=1,fill="red") +
  geom_point(data=DEG_table[which(DEG_table$ZGA_diff=="DOWN"),], 
             aes(x = log2FoldChange,y = -log10(pvalue)),size=1,fill="darkgreen") +
  theme_classic() +
  ggtitle(this_tile) + 
  theme(legend.position="none") +
  geom_vline(aes(xintercept=1),linetype="dashed") +
  geom_vline(aes(xintercept=-1),linetype="dashed") +
  geom_hline(aes(yintercept=-log10(0.05)),linetype="dashed") +
  xlab("log2 Foldchange") + ylab("-log10(P-value)") +
  theme_classic() +
  theme(legend.position = "none") + 
  scale_color_manual(values=c(NOT = "grey", UP = "red",
                              DOWN="darkgreen")) + 
  theme(
    axis.title.x = element_text(size=9),
    axis.title.y = element_text(size=9),
    axis.line = element_line(linewidth =0),
    panel.border = element_rect(color = "black", size = 0.8, fill = NA))
MA
ggsave(MA,filename="/home/qzh/Analysis/qd/MM_Eif1ad_KO/Fig4H.Mouse_Eif1ad_KD_Xiewei_list_Volcano.pdf",
       width = 3,height = 3.3)

# Heatmap
Xiewei_Diff_table <- read.table("/home/qzh/Analysis/qd/MM_Eif1ad_KO/02-result/Xiewei_OboxKO_DEseq_diff.txt",
                                sep = "\t",header = T)
Xiewei_ZGA_table <- read.table("/home/qzh/Analysis/qd/MM_Eif1ad_KO/02-result/zga_xie_obox_2023.txt",
                               sep = " ",header = T)

Xiewei_Down_table <- Xiewei_Diff_table[which(Xiewei_Diff_table$Group %in% "down-regulated"),]
Xiewei_Down_minor_table <- Xiewei_Down_table[which(Xiewei_Down_table$class %in% "minor"),]
Xiewei_Down_minor_table <- merge(Xiewei_Down_minor_table,Xiewei_ZGA_table,by.x="gene.name",by.y="gene_name")
Xiewei_Down_major_table <- Xiewei_Down_table[which(Xiewei_Down_table$class %in% "major"),]
Xiewei_Down_major_table <- merge(Xiewei_Down_major_table,Xiewei_ZGA_table,by.x="gene.name",by.y="gene_name")

Down_ZGA_gene <- DEG_table[which(DEG_table$ZGA_diff %in% "DOWN"),]
Minor_ZGA_gene <- Xiewei_ZGA_table[which(Xiewei_ZGA_table$class %in% c("minor")),]
Major_ZGA_gene <- Xiewei_ZGA_table[which(Xiewei_ZGA_table$class %in% c("major")),]

Minor_ZGA_Down <- Minor_ZGA_gene[which(Minor_ZGA_gene$gene_id %in% Down_ZGA_gene$gene),]
Major_ZGA_Down <- Major_ZGA_gene[which(Major_ZGA_gene$gene_id %in% Down_ZGA_gene$gene),]

# 
MM_2C_table <- read.table("/home/qzh/Analysis/qd/YY590_RNAseq_Eif1ad_gene_family/02-result/FPKM_dir/MM_L2C_FPKM.txt",
                          sep = "\t",header = T)
MM_2C_table <- MM_2C_table[,-2]
MM_2C_table <- MM_2C_table[!duplicated(MM_2C_table$Gene_id), ]

rownames(MM_2C_table) <- MM_2C_table[,1]
MM_2C_table <- MM_2C_table[,-1]

MM_2C_Sample_name <- lapply(names(MM_2C_table), 
                            function(col_name) unlist(strsplit(col_name, "_only"))[1])
MM_2C_Sample_name <- lapply(MM_2C_Sample_name, 
                            function(col_name) unlist(strsplit(col_name, "Smart3_"))[2])
names(MM_2C_table) <- MM_2C_Sample_name

MM_2C_table <- MM_2C_table[,which(names(MM_2C_table) %in% names(Eif1ad_table_yy))]
order_map <- match(names(Eif1ad_table_yy), 
                   names(MM_2C_table))
MM_2C_table <- MM_2C_table[,order_map]

MM_L2C_average <- MM_2C_table+1

Minor_MM_L2C_average <- MM_L2C_average[which(rownames(MM_L2C_average) %in% Minor_ZGA_gene$gene_id),]
Minor_MM_L2C_down <- Minor_MM_L2C_average[which(rownames(Minor_MM_L2C_average) %in% Minor_ZGA_Down$gene_id),]
Minor_MM_L2C_Nodown <- Minor_MM_L2C_average[which(!(rownames(Minor_MM_L2C_average) %in% Minor_ZGA_Down$gene_id)),]

Major_MM_L2C_average <- MM_L2C_average[which(rownames(MM_L2C_average) %in% Major_ZGA_gene$gene_id),]
Major_MM_L2C_down <- Major_MM_L2C_average[which(rownames(Major_MM_L2C_average) %in% Major_ZGA_Down$gene_id),]
Major_MM_L2C_Nodown <- Major_MM_L2C_average[which(!(rownames(Major_MM_L2C_average) %in% Major_ZGA_Down$gene_id)),]

MM_L2C_average <- rbind(Minor_MM_L2C_down,Minor_MM_L2C_Nodown,
                        Major_MM_L2C_down,Major_MM_L2C_Nodown)

MM_L2C_average_Zscore <- data.frame(t(scale(t(MM_L2C_average))))
MM_L2C_average_Zscore <- MM_L2C_average_Zscore %>%
  mutate_all(~ replace(., is.na(.), 0))
names(MM_L2C_average_Zscore) <- c("WT_L2C_rep1","WT_L2C_rep2","WT_L2C_rep3",
                                  "ASO_L2C_rep1","ASO_L2C_rep2","ASO_L2C_rep3")

heat_plot_r1 <- pheatmap::pheatmap(MM_L2C_average_Zscore,
                                   cluster_rows = F,
                                   cluster_cols = F,
                                   fontsize_col=13,
                                   fontsize_row = 14,
                                   show_rownames = F,
                                   show_colnames = T,
                                   border_color="white",
                                   width = 2.1,
                                   height = 5,
                                   gaps_row =c(64),
                                   gaps_col =c(3),
                                   color = colorRampPalette(c("#2b7cc0","#090504","#e6e64b"))(100),
                                   angle_col = "90")
ggsave(heat_plot,file="/home/qzh/Analysis/qd/MM_Eif1ad_KO/Fig4F.Mouse_Eif1ad_KD_Xiewei_list_heatmap_Zscore.pdf",
       width = 2.1,height = 5)

# Fold Change
ASO3_DEG_table <- read.table("/home/qzh/Analysis/qd/YY646_ASO3_Eif1ad/02-result/count/ASO3_vs_SCR_L2C_diff_Single_embryo.txt",
                             sep = "\t",header = T)
ASO_DEG_table <- ASO_DEG_table[,c(6,2)]

order_map <- match(rownames(MM_L2C_average_Zscore),
                   ASO_DEG_table$gene)
ASO_DEG_table <- ASO_DEG_table[order_map,]

rownames(ASO_DEG_table) <- ASO_DEG_table[,1]
ASO_DEG_table <- ASO_DEG_table[,-1]
ASO_DEG_table <- ASO_DEG_table %>%
  mutate_all(~ replace(., is.na(.), 0))
names(ASO_DEG_table) <- c("ASO3_L2C_FC")

heat_plot_r2 <- pheatmap::pheatmap(ASO_DEG_table,
                                   cluster_rows = F,
                                   cluster_cols = F,
                                   fontsize_col= 8,
                                   fontsize_row = 8,
                                   show_rownames = F,
                                   show_colnames = T,
                                   border_color="white",
                                   width = 1.1,
                                   height = 5,
                                   gaps_row =c(57),
                                   gaps_col =c(1),
                                   color = colorRampPalette(c("#5E4FA2","#3288BD","#3288BD","#66C2A5","#ABDDA4","#FEE08B","#9E0142"))(200),
                                   angle_col = "90")

heat_plot <- cowplot::plot_grid(heat_plot_r1$gtable, heat_plot_r2$gtable)
ggsave(heat_plot,file="/home/qzh/Analysis/qd/MM_Eif1ad_KO/Fig4.Mouse_Eif1ad_KD_Xiewei_list_heatmap_Zscore_and_FoldChange.pdf",
       width = 4,height = 5)

Xiewei_Diff_table <- read.table("/home/qzh/Analysis/qd/MM_Eif1ad_KO/02-result/zga_xie_obox_2023.txt",
                                sep = " ",header = T)
Minor_ZGA_gene <- Xiewei_Diff_table[which(Xiewei_Diff_table$class %in% "minor"),]
Major_ZGA_gene <- Xiewei_Diff_table[which(Xiewei_Diff_table$class %in% "major"),]

IP<-euler(list(Minor_ZGA=toupper(Minor_ZGA_gene$gene_id),
               Obox_KO_Minor=toupper(Xiewei_Down_minor_table$gene_id),
               Eif1ad_Aso_Minor=toupper(Minor_ZGA_Down$gene_id)))

Venn_plot_Minor <- plot(IP,
                        fills = list(fill=c("#bfbebe","#f5be38","#276fad"),
                                     alpha=1),
                        quantities = list(c("Obox_KO_Minor", 
                                            "Eif1ad_Aso_Minor",
                                            "Minor_ZGA"),
                                          col="black",
                                          type = "counts",
                                          cex=1.3),
                        # main = list(label="PN5",cex=1),
                        labels = list(col="black",font=1,cex=1.2,lineheight=1),
                        edges = list(col="white",alpha=1))

IP<-euler(list(Major_ZGA=unique(Major_ZGA_gene$gene_id),
               Obox_KO_Major=unique(Xiewei_Down_major_table$gene_id),
               Eif1ad_Aso_Major=unique(Major_ZGA_Down$gene_id)))
Venn_plot_Major <- plot(IP,
                        fills = list(fill=c("#bfbebe","#f5be38","#276fad"),
                                     alpha=1),
                        quantities = list(c("Major_ZGA",
                                            "Obox_KO_Major",
                                            "Eif1ad_Aso_Major"),
                                          col="black",
                                          type = "counts",
                                          cex=1.3),
                        # main = list(label="PN5",cex=1),
                        labels = list(col="black",font=1,cex=1.2,lineheight=1),
                        edges = list(col="white",alpha=1))

Venn_plot <- cowplot::plot_grid(
  Venn_plot_Minor,
  Venn_plot_Major,nrow = 2
)

ggsave(Venn_plot,filename="/home/qzh/Analysis/qd/MM_Eif1ad_KO/Fig4F2.Mouse_Eif1ad_KD_Veen.pdf",
       width = 4,height = 5)


MM_2C_table <- read.table("/home/qzh/Analysis/qd/YY590_RNAseq_Eif1ad_gene_family/02-result/FPKM_dir/MM_L2C_FPKM.txt",
                          sep = "\t",header = T)
MM_2C_table <- MM_2C_table[,-1]
MM_2C_table <- MM_2C_table[!duplicated(MM_2C_table$Gene_name), ]

rownames(MM_2C_table) <- MM_2C_table[,1]
MM_2C_table <- MM_2C_table[,-1]

MM_2C_Sample_name <- lapply(names(MM_2C_table), 
                            function(col_name) unlist(strsplit(col_name, "_only"))[1])
MM_2C_Sample_name <- lapply(MM_2C_Sample_name, 
                            function(col_name) unlist(strsplit(col_name, "Smart3_"))[2])
names(MM_2C_table) <- MM_2C_Sample_name

MM_2C_table <- MM_2C_table[,which(names(MM_2C_table) %in% names(Eif1ad_table_yy))]
order_map <- match(names(Eif1ad_table_yy), 
                   names(MM_2C_table))
MM_2C_table <- MM_2C_table[,order_map]
names(MM_2C_table) <- c("WT_L2C_rep1","WT_L2C_rep2","WT_L2C_rep3",
                        "ASO_L2C_rep1","ASO_L2C_rep2","ASO_L2C_rep3")
# MM_2C_table$Gene <- rownames(MM_2C_table)

Obox_gene_table1 <- c("Obox3","Obox4-ps35")
Obox_gene_table2 <- c("Obox1","Obox2","Obox5","Obox7","Duxf3")
Eif1ad_gene_table <- c("Eif1a","Eif1ax","Eif1ad","Eif1ad2","Eif1ad3","Eif1ad4","Eif1ad5",
                       "Eif1ad6","Eif1ad7","Eif1ad8","Eif1ad9",
                       "Eif1ad10","Eif1ad11","Eif1ad12","Eif1ad13","Eif1ad14",
                       "Eif1ad15","Eif1ad16","Eif1ad17","Eif1ad18","Eif1ad19")

Merge_display_gene <- c(rev(Eif1ad_gene_table),
                        rev(Obox_gene_table2),
                        rev(Obox_gene_table1))
# Filter data
MM_filter_table <- MM_2C_table[which(rownames(MM_2C_table) %in% Merge_display_gene),]
order_map <- match(Merge_display_gene, 
                   rownames(MM_filter_table))
MM_filter_table <- MM_filter_table[order_map,]
MM_filter_table$Gene_name <- rownames(MM_filter_table)
MM_filter_table[which(MM_filter_table$Gene_name %in% "Obox4-ps35"),"Gene_name"] <- "Obox4"
MM_filter_table[which(MM_filter_table$Gene_name %in% "Duxf3"),"Gene_name"] <- "Dux"
MM_filter_table$Gene_name <- factor(MM_filter_table$Gene_name, levels=MM_filter_table$Gene_name)

MM_filter_table_melt<-reshape2::melt(
  MM_filter_table,
  id.vars=c("Gene_name"),#要保留的主字段
  variable.name = "Cell_type",#转换后的分类字段名称（维度）
  value.name = "FPKM" #转换后的度量值名称
)
level_name <- unique(MM_filter_table_melt$Cell_type)
MM_filter_table_melt$Cell_type <- factor(MM_filter_table_melt$Cell_type, levels=level_name)
MM_filter_table_melt$log2_FPKM <- log2(MM_filter_table_melt$FPKM+1)

p1 <- ggplot(MM_filter_table_melt,aes(x=Cell_type,y=Gene_name,colour=log2_FPKM)) +
  scale_color_gradientn(values = seq(0,1,0.2),colours = c("#ffffff","#3c8abd","#083061")) +
  theme_bw()+
  geom_point(size = MM_filter_table_melt$log2_FPKM/4, stroke = 5) +
  theme(panel.grid = element_blank(),
    axis.text.x = element_text(angle =90,hjust =1,vjust = 0.5,
                               size=14, color = "black"),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size=13, color = "black")) +
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid")) 
p1
ggsave(p1,filename="/home/qzh/Analysis/qd/MM_Eif1ad_KO/Mouse_Eif1ad_KD_Obox_Eif1ad_heatmap.pdf",
       width = 4,height = 9)
