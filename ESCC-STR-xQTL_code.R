#Figure 1
## B
#write.csv(QTL_Num,file = "G:/STR图对应的数据和代码整理/fig1-b_QTL_Num画图数据.csv",row.names = F)
QTL_Num <- read.csv(file ="./fig1-b_QTL_Num画图数据.csv")
QTL_Num$TR_xQTL <- factor(QTL_Num$TR_xQTL,levels = c("aQTL","eQTL","sQTL","mQTL"))
library(ggplot2)
library(ggbreak)

pdf("E:/swy/EnsembleTR/QTLs_counts_State_new.pdf",width = 6,height = 8)
ggplot(QTL_Num, aes(x = TR_xQTL, y = Num, fill = TR_xQTL)) + 
  geom_bar(stat = "identity") + 
  geom_text(aes(label = Num), vjust = -0.2, size = 5)+# 使用原始值
  scale_fill_manual(values = c("#f1c0be", "#acbedc","#DEAAD4","#F1C89A")) +  # 手动设置颜色
  theme_minimal() + 
  scale_y_break(c(6500,50000), scales = "fixed", expand = expansion(add = c(0, 500)))+
  labs(title = "STR-xQTL Counts", x = "QTL Type", y = "xQTLs") +
  theme_minimal() +  # 使用无边框的主题
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()


## C
#write.csv(qq_sampled,file = "/share/pub/sunwy/sunwy/EnsembleTR_result/STR_SNP/LD_result/fig1-c_QQ_plot.csv",row.names = F)
library(ggplot2)
qq_sampled  <- read.csv("./fig1-c_QQ_plot.csv")
#Q-Q plot
p <- ggplot(qq_sampled, aes(x = expected, y = observed, color = group, shape = group)) +
  geom_point(size = 1, alpha = 1) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 0.8) +
  labs(x = "Expected -log10(p)", 
       y = "Observed -log10(p)",
       title = "Q-Q Plot of STR-xQTLs") +
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "top"
  ) +
  scale_color_manual(values = c("aQTL" = "#f1c0be", 
                                "eQTL" = "#acbedc",
                                "sQTL" = "#DEAAD4",
                                "mQTL" = "#F1C89A")) +
  scale_shape_manual(values = c("sQTL" = 16, "aQTL" = 17, "eQTL" = 15, "mQTL" = 18))
# 显示或保存图像
print(p)
ggsave("QQ-plot_new.pdf", plot = p, width = 3.5, height = 4.2)


## D
library(dplyr)
library(ggplot2)
eqtl_HORMAD1 <- read.csv("./fig1_d_eQTL数据.csv")
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4")
sort(unique(eqtl_HORMAD1$geno))
eqtl_HORMAD1$geno <- factor(as.character(eqtl_HORMAD1$geno), levels = sort(as.numeric(unique(eqtl_HORMAD1$geno))))

pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_HORMAD1_susieR.pdf",width =4,height = 5)
ggplot(data = eqtl_HORMAD1, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of HORMAD1 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()
####sQTL####
dat <- read.csv("./fig1_d_sQTL数据.csv")
sort(unique(dat$geno))
cell_type_cols <- c("#ED7E8B","#FFD19C","#B797C6","#699ECA", "#52854C","#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4")
dat$geno <- factor(as.numeric(dat$geno), levels = sort(as.numeric(unique(dat$geno))))

pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/sqtl_SEMA4G.pdf",width =3,height = 3)
ggplot(data = dat, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "Normalized \n splicing Ratio")+  # 修改坐标轴名称
  ggtitle("TR-sqtl of SEMA4G")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()
####mQTL####
mqtl_TAF1B <- read.csv("./fig1_d_mQTL数据.csv")
sort(unique(mqtl_TAF1B$geno))
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4")
mqtl_TAF1B$geno <- factor(as.character(mqtl_TAF1B$geno), levels = sort(as.numeric(unique(mqtl_TAF1B$geno))))

pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/mqtl_TAF1B.pdf",width =4,height = 5)
ggplot(data = mqtl_TAF1B, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "STR Length", y = "DNA methylation")+  # 修改坐标轴名称
  ggtitle("STR-mQTL of TAF1B in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## E
library(VennDiagram)
library(dplyr)
load("./fig1_e_韦恩图数据.RData")
venn.plot <- venn.diagram(
  x = STR_genes_list,
  category.names = c("aSTR", "mSTR", "eSTR", "sSTR"),  # ✅ 强制顺序
  force.unique = TRUE,
  filename = NULL,
  fill = c("#f1c0be","#F1C89A", "#acbedc", "#DEAAD4" ),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Four QTL STR-Gene Overlap"
)
pdf("STR_xQTL_Veen_new.pdf",width = 5,height = 5)
grid::grid.draw(venn.plot)
dev.off()



## F
library(ggplot2)
summary_200KB <- read.csv(file = "./fig1_f_LD画图数据.csv")
pdf(file = "/share/pub/sunwy/sunwy/EnsembleTR_result/STR_SNP/LD_result/STR-SNP_LD_R2_compare_0_200kb.pdf",width = 6,height = 4)
ggplot(summary_200KB, aes(x = distance, y = mean, color = type, shape = type)) +
  geom_point(size = 2) +
  # 如果想要平滑趋势线，可以加：
  geom_smooth(se = FALSE, method = "loess", span = 0.6) +
  theme_bw() +
  labs(
    x = "Distance (kb)",
    y = expression(LD~(r^2)),
    title = "LD decay: SNP-SNP vs STR-SNP"
  ) + 
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1),  # 添加黑色边框
        plot.title = element_text(hjust = 0.5))  # 居中标题

dev.off()


## G
library(ggridges)
library(ggplot2)
four_phenotype_h2  <- read.csv(file = "./fig1_g_STR遗传力数据.csv")
# 计算每组median
median_h2_phenotype <- four_phenotype_h2 %>%
  group_by(phenotype) %>%
  summarise(median_h2 = mean(h2))

four_phenotype_h2$phenotype <- factor(
  four_phenotype_h2$phenotype,
  levels = c("APA", "Gene", "Splicing", "CpG")
)

median_h2_phenotype$phenotype <- factor(
  median_h2_phenotype$phenotype,
  levels = c("APA", "Gene", "Splicing", "CpG")
)
pdf(file = "/share/pub/changxj/sunwy/STR_Heritability/STR对四种分子表型的遗传力.pdf",width = 5,height = 3)
ggplot(four_phenotype_h2, aes(x = h2, fill = phenotype, color = phenotype)) +
  geom_density(alpha = 0.6, adjust = 1.5, bw = 0.1, size = 0.8) +
  geom_vline(
    data = median_h2_phenotype,
    aes(xintercept = median_h2, color = phenotype),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  )+
  
  # ✅ 加均值文字
  geom_text(
    data = median_h2_phenotype,
    aes(x = 0.3, y = 3, label = round(median_h2, 3)),
    color = "black",
    vjust = 2,
    size = 5
  ) +
  facet_grid(phenotype ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(
    "APA" = "#f1c0be",
    "Gene" = "#acbedc",
    "Splicing" = "#DEAAD4",
    "CpG" = "#F1C89A"
  )) +
  scale_color_manual(values = c(
    "APA" = "#f1c0be",
    "Gene" = "#acbedc",
    "Splicing" = "#DEAAD4",
    "CpG" = "#F1C89A"
  ))+
  labs(
    x = "STR-associated heritability of four molecular phenotypes",
    y = NULL
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    text = element_text(color = "black"),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.2, "lines")
  )
dev.off()


#Figure 2
## A
all_genome_or <- read.csv("./fig2_a_data.csv")
all_genome_or$phenotype <- factor(all_genome_or$phenotype,levels = rev(c("aSTR","eSTR","sSTR","mSTR")))
all_genome_or$row_id <- as.numeric(all_genome_or$Element)
library(ggplot2)
bg_data <- data.frame(
  Element = levels(all_genome_or$Element)
)
bg_data$row_id <- seq_len(nrow(bg_data))
pdf(file = "E:/swy/EnsembleTR/合并四种分子表型QTL显著的STR对应基因组元件结果图_CI_log.pdf",width = 6,height = 6)
pd <- position_dodge(width = 0.6)
ggplot(all_genome_or, aes(x = logOR, y = Element, color = phenotype)) +
  
  # 👉 斑马线背景
  geom_rect(
    data = bg_data,
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = row_id - 0.5,
      ymax = row_id + 0.5
    ),
    inherit.aes = FALSE,
    fill = rep(c("grey95", "white"), length.out = length(unique(all_genome_or$Element)))
  ) +
  
  geom_point(position = pd, size = 3) +
  
  geom_errorbarh(
    aes(xmin = log_CI_low, xmax = log_CI_high),
    position = pd,
    height = 0.2
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  scale_color_manual(values = c(
    aSTR = "#E58B87",
    eSTR = "#89A3CD",
    sSTR = "#D28AC4",
    mSTR = "#E79F57"
  )) +
  
  labs(
    x = "Odds Ratio",
    y = "",
    color = "Phenotype"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 13, colour = "black"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )
dev.off()

## B
all_motif_or <- read.csv(file = "./fig2_b_data1_all_motif_or.csv")
lm_results <- read.csv(file = "./fig2_b_data2_lm_results.csv")
library(ggplot2)
all_motif_or$phenotype <- factor(
  all_motif_or$phenotype,
  levels = c("aSTR","eSTR","sSTR","mSTR")
)
ggplot(all_motif_or, aes(x = motif_length, y = OR, color = phenotype)) +
  geom_point(size = 3) +
  annotate(geom = 'rect',xmin = -Inf, xmax = Inf,ymin = 1, ymax = 2,fill = "#F1C89A",alpha = 0.15 ) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "#F1C89A",
    size = 0.8
  ) +
  
  geom_smooth(
    method = lm,
    se = FALSE
  ) +
  
  scale_color_manual(values = c(
    "#f1c0be",  # aSTR
    "#acbedc",  # eSTR
    "#DEAAD4",  # mSTR
    "#E79F57"   # sSTR
  )) +
  
  scale_y_continuous(
    breaks = seq(0,2,0.5),
    expand = c(0,0.01),
    limits = c(0,2)
  ) +
  
  labs(
    x = "Motif length",
    y = "Odds Ratio",
    title = "Motif size features of STRs"
  ) +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 15, colour = "black"),
    axis.line = element_line(size = 0.5),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

## C
library(ggplot2)
eSTR_motif_or_results_top10 <- read.csv(file = "./fig2_c_data2_eSTR_motif_or_results_top10.csv")
eSTR_motif_or_results_top10$Motif <- factor(eSTR_motif_or_results_top10$Motif, 
                                            levels = eSTR_motif_or_results_top10$Motif[order(eSTR_motif_or_results_top10$OR, decreasing = TRUE)])

pdf("E:/swy/EnsembleTR/eQTL/q005/motif_unit_OR.pdf",width = 6,height = 4)
ggplot(eSTR_motif_or_results_top10, aes(x = Motif, y = OR, fill = sig)) + 
  geom_bar(stat = "identity") + 
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1,
           fill = "#e66c2a", alpha = 0.2, color = '#ffd2ba') +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "#e66c2a", size = 0.8) +
  #scale_y_continuous(breaks = seq(0, 3, 1), expand = c(0, 0.01), limits = c(0, 3)) +
  
  # 自定义显著性颜色
  scale_fill_manual(values = c("Significant" = "#acbedc", "Not significant" = "grey80")) +
  
  labs(title = "STR-eQTL Motif unit", x = "Motif unit", y = "Odds Ratio") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.8, size = 14, angle = 30),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),
    # 加上刻度线：
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()
aSTR_motif_or_results_top10 <- read.csv(file = "./fig2_c_data1_aSTR_motif_or_results_top10.csv")
# 修改Motif的因子顺序，使得柱形图按OR值大小排序
aSTR_motif_or_results_top10$Motif <- factor(aSTR_motif_or_results_top10$Motif, 
                                            levels = aSTR_motif_or_results_top10$Motif[order(aSTR_motif_or_results_top10$OR, decreasing = TRUE)])
#####画图用的互补碱基的motif#####
pdf("E:/swy/EnsembleTR/aQTL/q005/motif_unit_OR_new03.02.pdf",width = 6,height = 4)
ggplot(aSTR_motif_or_results_top10, aes(x = Motif, y = OR, fill = sig)) + 
  geom_bar(stat = "identity") + 
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1,
           fill = "#e66c2a", alpha = 0.2, color = '#ffd2ba') +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "#e66c2a", size = 0.8) +
  #scale_y_continuous(breaks = seq(0, 3, 1), expand = c(0, 0.01), limits = c(0, 3)) +
  
  # 自定义显著性颜色
  scale_fill_manual(values = c("Significant" = "#f1c0be", "Not significant" = "grey80")) +
  
  labs(title = "STR-aQTL Motif unit", x = "Motif unit", y = "Odds Ratio") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0, size = 14, angle = 90),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),
    # 加上刻度线：
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()
library(ggplot2)
sSTR_motif_OR_results_top10 <- read.csv(file = "./fig2_c_data3_sSTR_motif_or_results_top10.csv")
# 修改Motif的因子顺序，使得柱形图按OR值大小排序
sSTR_motif_OR_results_top10$Motif <- factor(sSTR_motif_OR_results_top10$Motif, 
                                            levels = sSTR_motif_OR_results_top10$Motif[order(sSTR_motif_OR_results_top10$OR, decreasing = TRUE)])
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/sQTL_result/sSTR_anno/motif_unit_OR.pdf",width = 7,height = 5)
ggplot(sSTR_motif_OR_results_top10, aes(x = Motif, y = OR, fill = Motif)) + 
  geom_bar(stat = "identity",fill = "#DEAAD4") + 
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf,    #设置y>1背景图层的颜色、范围、透明度       
           ymin =1, ymax=1,           
           fill = "#e66c2a", alpha = 0.2, color = '#ffd2ba')+  
  geom_hline(yintercept = 1,linetype = "dashed",     #设置y=1的虚线        
             color = "#e66c2a",size=0.8)+
  scale_y_continuous(breaks = seq(0,6,1),                     
                     expand = c(0,0.01), limits = c(0,6))+
  
  theme_minimal() + 
  labs(title = "STR-sQTL Motif unit", x = "Motif unit", y = "Odds Ratio") +
  theme_minimal() +  # 使用无边框的主题
  theme(
    axis.text.x = element_text(hjust = 0.8, size = 14,angle = 30),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),  # 仅保留坐标轴线
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()
mSTR_motif_or_results_top10 <- read.csv(file = "./fig2_c_data4_mSTR_motif_or_results_top10.csv")
library(ggplot2)
# 修改Motif的因子顺序，使得柱形图按OR值大小排序
mSTR_motif_or_results_top10$Motif <- factor(mSTR_motif_or_results_top10$Motif, 
                                            levels = mSTR_motif_or_results_top10$Motif[order(mSTR_motif_or_results_top10$OR, decreasing = TRUE)])
pdf(file = "E:/swy/EnsembleTR/mQTL/CpG_STR_motif_unit.pdf",width = 6,height = 4)
ggplot(mSTR_motif_or_results_top10, aes(x = Motif, y = OR, fill = sig)) + 
  geom_bar(stat = "identity") + 
  annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1,
           fill = "#e66c2a", alpha = 0.2, color = '#ffd2ba') +  
  geom_hline(yintercept = 1, linetype = "dashed", color = "#e66c2a", size = 0.8) +
  #scale_y_continuous(breaks = seq(0, 3, 1), expand = c(0, 0.01), limits = c(0, 3)) +
  
  # 自定义显著性颜色
  scale_fill_manual(values = c("Significant" = "#F1C89A", "Not significant" = "grey80")) +
  
  labs(title = "STR-mQTL Motif unit", x = "Motif unit", y = "Odds Ratio") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.8, size = 14, angle = 30),
    axis.text.y = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),
    # 加上刻度线：
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(3, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.title = element_blank(),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()


## D
library(ggplot2)
library(dplyr)
library(forcats)
OR_long_plot <- read.csv(file = "./fig2_d_TF_OR_STR画图数据.csv")
TF_order <- OR_long_plot %>%
  group_by(TF_name) %>%
  slice_max(logOR, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(phenotype, desc(logOR))

OR_long_plot$TF_name <- factor(OR_long_plot$TF_name,
                               levels = TF_order$TF_name)

pdf("/share/pub/changxj/sunwy/TF_bed/top10_TF_group.pdf")
ggplot(OR_long_plot, aes(x = phenotype, y = TF_name, fill = logOR)) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = "#379fb4", mid = "white", high = "#b379b4", midpoint = 1
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(hjust = 0.5, colour = "black",size = 12),
    axis.text.y = element_text(colour = "black",size = 12),
    panel.border = element_blank()
  )
dev.off()

## F
library(ggplot2)
plot_df <- read.csv(file = "./fig2_f_bmediatR_画图数据.csv")
ggplot(
  plot_df,
  aes(
    x = Type,
    y = Num,
    fill = Model
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8,
    position = position_dodge(width = 0.9)
  ) +
  geom_text(
    aes(label = Label),
    position = position_dodge(width = 0.9),
    hjust = -0.1,
    size = 4
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "co-local" = "#A29EC1",
      "partial med" = "#609ACA",
      "complete med" = "#EC6A90"
    )
  ) +
  theme_classic() +
  labs(
    x = "",
    y = "Number of triplets"
  )

## G
library(dplyr)
library(ggplot2)
eqtl_TREX1 <- read.csv(file = "./fig2_g_eqtl_TREX1画图数据.csv")
eqtl_TREX1$geno_new <- factor(as.numeric(eqtl_TREX1$geno_new), levels = sort(as.numeric(unique(eqtl_TREX1$geno_new))))
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4","#009E73","#FFDB6D")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_TREX1_bmeditaR.pdf",width =4,height = 5)
ggplot(data = eqtl_TREX1, aes(x = geno_new, y = exp, fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "Sum number of STR units", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of TREX1 in ESCC")+  theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14)
  )
dev.off()

## H
library(ggplot2)
mqtl_TREX1 <- read.csv(file = "./fig2_h_mqtl_TREX1画图数据.csv")
mqtl_TREX1$geno_new <- factor(as.character(mqtl_TREX1$geno_new), levels = sort(as.character(unique(mqtl_TREX1$geno_new))))
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4","#009E73","#FFDB6D")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/mqtl_TREX1_bmeditaR.pdf",width = 4,height = 5)
ggplot(data = mqtl_TREX1, aes(x = geno_new, y = exp, fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "CpG")+  # 修改坐标轴名称
  ggtitle("TR-mQTL of TREX1 in ESCC")+theme_bw() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14)
  )
dev.off()

## I
library(ggplot2)
TREX1_tr_CpG <- read.csv(file = "./fig2_i_TREX1_tr_CpG.csv")
# 3. 残差 vs TR长度的相关分析
cor.test(TREX1_tr_CpG$resid_TREX1, as.numeric(TREX1_tr_CpG$TR_TREX1))
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA")
TREX1_tr_CpG$geno_new <- factor(as.character(TREX1_tr_CpG$geno_new), levels = sort(as.character(unique(TREX1_tr_CpG$geno_new))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/resid_mqtl_TREX1_bmeditaR.pdf",width = 4,height = 5)
ggplot(TREX1_tr_CpG, aes(x = geno_new, y = resid_TREX1,fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "Residual Exp")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of HLA-DRB5 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))

dev.off()


#Figure 3
## A
library(VennDiagram)
library(dplyr)
load("./fig3_a_韦恩图画图数据.RData")
venn.plot <- venn.diagram(
  x = STR_genes_list,
  category.names = c("aGene", "mGene", "eGene", "sGene"),
  force.unique = TRUE,
  filename = NULL,
  fill = c("#f1c0be", "#F1C89A", "#acbedc", "#DEAAD4"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.2,
  main = "Four QTL xGene Overlap"
)

pdf("E:/swy/EnsembleTR/四种分子表型QTL对应基因的Veen.pdf",width = 5,height = 5)
grid::grid.draw(venn.plot)
dev.off()

## B
library(ggplot2)
GO_bp_long <- read.csv(file = "./fig3_b_GO_bp画图数据.csv")
pdf("E:/swy/EnsembleTR/new_GO_bp_dotplot_plot2.pdf",  width = 11, height = 7)
ggplot(GO_bp_long, aes(x = QTLtype, y = Term, color = log_p, size = log_combined)) +
  geom_point() +
  scale_color_gradient(low = "grey", high = "#E54924") +  # 颜色梯度
  scale_size_continuous(range = c(2, 8)) +  # 点的大小
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5,size = 14),  # 旋转 X 轴文本
    axis.title.x = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 添加黑色框线
  ) +
  labs(x = "QTL Type", y = "GO Biological Process", 
       color = "-log10(p)", size = "Combined.Score", 
       title = "GO Biological Processes across STR-xQTL")

dev.off()


## C

library(ggplot2)
library(dplyr)
custom_colors <- c("#f1c0be", "#acbedc","#DEAAD4","#F1C89A")
CD4_df <- read.csv(file = "./fig3_c_CD4_画图数据.csv")
pdf("E:/swy/EnsembleTR/Immune_cell_marker/CD4cell_P.pdf",width = 5,height = 6)
ggplot(CD4_df, aes(x = QTL_type, y = log10P, fill = QTL_type)) + 
  geom_bar(stat = "identity") + 
  #geom_text(aes(label = log10P), vjust = -0.2, size = 5)+# 使用原始值
  scale_fill_manual(values = custom_colors) +  # 手动设置颜色
  theme_bw() + 
  labs(title = "CD4+ T cell",x = "QTL Type", y = "-log10(P)") + # 使用无边框的主题
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()
CD8_df <- read.csv(file = "./fig3_c_CD8_画图数据.csv")
pdf("E:/swy/EnsembleTR/Immune_cell_marker/CD8cell_P.pdf",width = 5,height = 6)
ggplot(CD8_df, aes(x = QTL_type, y = log10P, fill = QTL_type)) + 
  geom_bar(stat = "identity") + 
  #geom_text(aes(label = log10P), vjust = -0.2, size = 5)+# 使用原始值
  scale_fill_manual(values = custom_colors) +  # 手动设置颜色
  theme_bw() + 
  labs(title = "CD8+ T cell",x = "QTL Type", y = "-log10(P)") + # 使用无边框的主题
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()
NKT_df <- read.csv(file = "./fig3_c_NKT_画图数据.csv")
pdf("E:/swy/EnsembleTR/Immune_cell_marker/NKT_cell_P.pdf",width = 5,height = 6)
ggplot(NKT_df, aes(x = QTL_type, y = log10P, fill = QTL_type)) + 
  geom_bar(stat = "identity") + 
  #geom_text(aes(label = log10P), vjust = -0.2, size = 5)+# 使用原始值
  scale_fill_manual(values = custom_colors) +  # 手动设置颜色
  theme_bw() + 
  labs(title = "NK-T cell",x = "QTL Type", y = "-log10(P)") + # 使用无边框的主题
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()


## D
library(pheatmap)
QTL_num_ordered <- read.csv(file = "./fig3_d_左图.csv")
heat_colors <- c("white", "#f1c0be","#acbedc","#DEAAD4","#F1C89A")
pdf(file = "E:swy/EnsembleTR/maf_QTL_Gene_type_heatmap.pdf",height = 4,width = 2.3)
pheatmap(QTL_num_ordered,
         color = heat_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         legend = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         border_color = "grey90",
         breaks = c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5))  # 关键：确保颜色与数值一一对应
dev.off()
library(ComplexHeatmap)
library(circlize)
library(maftools)
load("./fig3_d瀑布图.RData")
pdf(file = "E:swy/EnsembleTR/maf_QTL_Gene_type_瀑布图.pdf",height = 6,width = 8)
oncoplot(
  maf = QTL_maf_object,
  genes = QTL_Gene,  # 使用所有基因
  draw_titv = FALSE,
  sortByAnnotation = TRUE,
  removeNonMutated = TRUE
)
dev.off()

## E
library(ggplot2)
all_pathway <- read.csv(file = "./fig3_e药物富集数据.csv")
all_pathway$QTL <- factor(all_pathway$QTL,levels = c("aQTL","eQTL","sQTL","mQTL"))
pdf("E:/swy/Drug/four_QTL_Drug_cancer_pathway.pdf",width = 8,height = 5)
ggplot(all_pathway, aes(x = QTL, y = PATHWAY_NAME, fill = proportion)) +  
  geom_tile(color = "white", width = 0.95, height = 0.95) +
  scale_fill_gradientn(
    colours = c("#4575B4", "#E0F3F8", "#FFFFBF", "#FEE090", "#D73027"),
    values = scales::rescale(c(0, 0.2, 0.4, 0.6, 0.8)),
    limits = c(0, 0.8),
    name = "The proportion\nof significant\nassociated genes"
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 11, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid = element_blank()
  )
dev.off()


#Figure 4
## C
library(dplyr)
library(ggplot2)
aqtl_PRH1 <- read.csv("./fig4_c_aQTL_PRH1.csv")
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4")
aqtl_PRH1$geno <- factor(as.character(aqtl_PRH1$geno), levels = sort(as.character(unique(aqtl_PRH1$geno))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/aqtl_PRH1_susieR.pdf",width =4,height = 5)
ggplot(data = aqtl_PRH1, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "APA level")+  # 修改坐标轴名称
  ggtitle("TR-aQTL of PRH1")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()

## D
library(ggplot2)
PRH1_susieR <- read.csv(file = "./fig4_d_down_susieR_PRH1.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/aSTR_PRH1_susieR.pdf",width =6,height =4)
ggplot(PRH1_susieR, aes(x = variant_id, y = pip, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of PRH1",
       x = "Genomic Position", y = "Casual Probability (PIP)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()   # 去除次网格线
  )
dev.off()
PRH1_susieR <- read.csv(file = "./fig4_d_up_susieR_p_PRH1.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/aSTR_PRH1_susieR_log10p.pdf",width =6,height =4)
ggplot(PRH1_susieR, aes(x = variant_id, y = log10p, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  #scale_x_continuous(limits = c(x_min, x_max), labels = scales::comma) +
  labs(title = "Fine-mapping of PRH1",
       x = "Genomic Position", y = "-log10(P)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size=8, colour = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dev.off()

## F
library(ggplot2)
eqtl_ULK4 <- read.csv(file = "./fig4_f_eQTL_ULK4.csv")
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4")
eqtl_ULK4$geno <- factor(as.character(eqtl_ULK4$geno), levels = sort(as.character(unique(eqtl_ULK4$geno))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_ULK4_susieR.pdf",width =4,height = 5)
ggplot(data = eqtl_ULK4, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of ULK4 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## G
library(ggplot2)
ULK4_susieR <-read.csv(file = "./fig4_g_down_susieR_ULK4.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eSTR_ULK4_susieR.pdf",width =4,height =3)
ggplot(ULK4_susieR, aes(x = variant_id, y = pip, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of ULK4",
       x = "Genomic Position", y = "Casual Probability (PIP)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()   # 去除次网格线
  )
dev.off()
ULK4_susieR <- read.csv(file = "./fig4_g_up_susieR_p_ULK4.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eSTR_ULK4_susieR_log10p.pdf",width =4,height =3)
ggplot(ULK4_susieR, aes(x = variant_id, y = log10p, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of ULK4",
       x = "Genomic Position", y = "-log10(P)") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()   # 去除次网格线
  )
dev.off()

## H
library(ggplot2)
all_site_to_TR <- read.csv(file = "./fig4_h_distance画图数据.csv")
all_site_to_TR$type <- factor(
  all_site_to_TR$type,
  levels = c("fine_mapped_aSTR", "fine_mapped_eSTR", "fine_mapped_sSTR", "fine_mapped_mSTR")
)
type_cols <- c(
  fine_mapped_aSTR = "#E8A0A0",
  fine_mapped_eSTR = "#8FAAD3",
  fine_mapped_sSTR = "#C994C7",
  fine_mapped_mSTR = "#F2BE8D"
)
ggplot(all_site_to_TR,
       aes(x = bin_center,
           y = relative_enrichment,
           color = type)) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    se = FALSE,
    linewidth = 1.0
  ) +
  facet_grid(type ~ ., scales = "free_y") +
  scale_color_manual(values = type_cols) +
  labs(
    title = "",
    x = "Distance to Regulatory site (Kb)",
    y = "Fold Enrichment"
  ) +
  scale_x_continuous(
    breaks = seq(-400, 400, by = 200),
    limits = c(-500, 500)
  ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.y = element_text(size = 9, angle = 0, hjust = 0),
    axis.text.x = element_text(size = 9, colour = "black"),
    axis.text.y = element_text(size = 8, colour = "black"),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.line = element_line(size = 0.8),
    plot.title = element_text(hjust = 0.5, size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1, "lines")
  )

## I
library(ggplot2)
library(tidyr)
GO_bp_long <- read.csv(file = "./fig4_i_fine_mapped_gene富集数据.csv")
GO_bp_long$Term <- factor(GO_bp_long$Term, levels = rev(unique(GO_bp_long$Term)))
pdf("E:/swy/EnsembleTR/fine_mapped_STR_GO_bp_dotplot_plot.pdf",  width = 11, height = 7)
ggplot(GO_bp_long, aes(x = QTLtype, y = Term, color = log_p, size = Combined.Score)) +
  geom_point() +
  scale_color_gradient(low = "grey", high = "#E54924") +  # 颜色梯度
  scale_size_continuous(range = c(2, 8)) +  # 点的大小
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5,size = 14),  # 旋转 X 轴文本
    axis.title.x = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 添加黑色框线
  ) +
  labs(x = "QTL Type", y = "GO Biological Process", 
       color = "-log10(p)", size = "Combined.Score", 
       title = "GO Biological Processes across fine-mapped STR-xQTL")

dev.off()


#Figure 5
## C
library(dplyr)
library(ggplot2)
eqtl_ABHD11_AS1_new <- read.csv(file = "./fig5_c_eQTL_ABHD11_AS1.csv")
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4")
eqtl_ABHD11_AS1_new$geno <- factor(as.character(eqtl_ABHD11_AS1_new$geno), levels = sort(as.character(unique(eqtl_ABHD11_AS1_new$geno))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_ABHD11_AS1_susieR.pdf",width =4,height = 5)
ggplot(data = eqtl_ABHD11_AS1_new, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of ABHD11_AS1 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## E
library(ggplot2)
input <- read.csv(file = "./fig5_e_STR_SNP各自p值数据.csv")
highlight_variant <- "chr7:73740687"
# 设置颜色阶梯
pdf("E:/swy/EnsembleTR/ABHD11_AS1_coloc_eSTR_pp4.pdf",width = 8,height = 5)
ggplot(input, aes(x = log10p_GWAS, y = log10p_sQTL, color = R2)) +
  geom_point(size = 2) +
  geom_text(
    data = subset(input, SNP_ID == highlight_variant),
    aes(x = log10p_GWAS, y = log10p_sQTL, label = highlight_variant),
    vjust = -1, color = "red", size = 3
  ) +
  
  scale_color_gradientn(colors = c("blue", "green", "yellow", "red"),
                        limits = c(0, 1),
                        name = expression(r^2)) +
  labs(x = expression(-log[10](P)~"of ESCC GWAS"),
       y = expression(-log[10](P)~"of eQTL"),
       title = "PP4 = 0.87304411") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

## F
library(ggplot2)
ABHD11_AS1_susieR <- read.csv(file = "./fig5_f_down_susieR_ABHD11_AS1.csv")
ABHD11_AS1_susieR$variant <-  factor(ABHD11_AS1_susieR$variant, levels = sort(unique(ABHD11_AS1_susieR$variant)))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eSTR_ABHD11_AS1_susieR.pdf",width =4,height =3)
ggplot(ABHD11_AS1_susieR, aes(x = variant_id, y = pip, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of ABHD11_AS1",
       x = "Genomic Position", y = "Casual Probability (PIP)") +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  )
dev.off()
ABHD11_AS1_susieR <- read.csv(file = "./fig5_f_up_susieR_p_ABHD11_AS1.csv")
ABHD11_AS1_susieR$variant_id <- as.numeric(sub(".*:", "", ABHD11_AS1_susieR$variant))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eSTR_ABHD11_AS1_susieR_log10p.pdf",width =4,height =3)
ggplot(ABHD11_AS1_susieR, aes(x = variant_id, y = log10p, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of ABHD11_AS1",
       x = "Genomic Position", y = "-log10(P)") +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  )
dev.off()

## G
library(ggplot2)
library(dplyr)
plot_df <- read.csv(file = "./fig5_g_STR实验数据_新.xlsx")
plot_df$group <- factor(plot_df$group, levels = c("(AC)10", "(AC)20"))

# 计算均值和SEM
stat_df <- plot_df %>%
  group_by(group) %>%
  summarise(
    mean = mean(value),
    sem = sd(value) / sqrt(n()),
    .groups = "drop"
  )

p_value <- t.test(value ~ group, data = plot_df)$p.value
p_label <- ifelse(p_value < 0.001, "***",
                  ifelse(p_value < 0.01, "**",
                         ifelse(p_value < 0.05, "*", "ns")))
pdf("ABHD11_luciferase_barplot.pdf", width = 2.4, height = 2.4)
ggplot(plot_df, aes(x = group, y = value, fill = group)) +
  geom_col(
    data = stat_df,
    aes(x = group, y = mean, fill = group),
    width = 0.55,
    color = "black",
    linewidth = 0.5
  ) +
  geom_errorbar(
    data = stat_df,
    aes(x = group, ymin = mean - sem, ymax = mean + sem),
    width = 0.18,
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.08,
    size = 1.6,
    alpha = 0.8
  ) +
  scale_fill_manual(values = c("(AC)10" = "grey80", "(AC)20" = "#EF767A")) +
  scale_color_manual(values = c("(AC)10" = "#7ECF92", "(AC)20" = "#D95F5F")) +
  scale_y_continuous(
    limits = c(0, 1.5),
    breaks = seq(0, 1.5, 0.5),
    expand = c(0, 0)
  ) +
  # 显著性横线
  annotate("segment", x = 1, xend = 1, y = 1.25, yend = 1.38, linewidth = 0.5) +
  annotate("segment", x = 2, xend = 2, y = 1.25, yend = 1.38, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 2, y = 1.38, yend = 1.38, linewidth = 0.5) +
  annotate("text", x = 1.5, y = 1.42, label = p_label, size = 4) +
  labs(
    x = "",
    y = "Relative luminescence"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 11, color = "black"),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5)
  )

dev.off()

## H
library(survminer)
library(survival)
ESCC_2_STR_exp_clin <- read.csv(file = "./ABHD11_AS1按照STR是否扩增分高低两组的临床信息数据.csv")
ESCC_2_STR_exp_clin$new_group <- ifelse(ESCC_2_STR_exp_clin$Group < 20, "Low",
                                        ifelse(ESCC_2_STR_exp_clin$Group > 20, "High", "Equal"))
ESCC_2_STR_exp_clin <- ESCC_2_STR_exp_clin[ESCC_2_STR_exp_clin$new_group !="Low",]
surv_obj <- Surv(ESCC_2_STR_exp_clin$day, ESCC_2_STR_exp_clin$Status)
diff <- survdiff(surv_obj ~ new_group, data = ESCC_2_STR_exp_clin)
pvalue <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
print(paste0("###### 生存差异 p-value = ", round(pvalue, 4), " ######"))
fit <- survfit(surv_obj ~ new_group, data = ESCC_2_STR_exp_clin)
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/ABHD11_AS1_survival_opmit_High.pdf",width = 6,height = 6)
ggsurvplot(
  fit,
  data = ESCC_2_STR_exp_clin,
  pval = TRUE, # 是否显示P值，TRUE可以换为自己想要表达的字符串
  pval.method = TRUE, # 是否展示P值的检验方法
  legend.title='STRs amplification proportion', # 自定义图例的标题
  # legend.labs=c('Male','Female'), # 自定义分组变量的名字
  conf.int = F, # 是否展示置信区间
  risk.table = TRUE, # 是否增加 risk table
  risk.table.col = "strata", # 可以根据分组改变risk table的字体颜色
  linetype = "strata", # 可以根据分组改变分组曲线的样式
  #surv.median.line = "hv", # 指出中位生存率
  ggtheme = theme_bw(), # 设置主题
  palette = c("#ED7E8B","#FFD19C") # 设置曲线的颜色
)
dev.off()


## I
library(ggplot2)
ESCC_2_STR_exp_clin_prop <- read.csv(file = "./ABHD11_AS1_TNM分期数据.csv")
ESCC_2_STR_exp_clin_prop$new_group <- factor(ESCC_2_STR_exp_clin_prop$new_group,levels = c("Equal","High"))
####卡方检验计算p值
TNM_tab <- xtabs(n ~ new_group + TNM_stage,
                 data = ESCC_2_STR_exp_clin_prop)
TNM_p_value <- chisq.test(TNM_tab, correct = FALSE)$p.value
TNM_p_label <- paste0("P = ", signif(TNM_p_value, 3))
TNM_p_label
#save(os_baidu_1_2_cli_prop,file = "/share/pub/sunwy/sunwy/EnsembleTR_result/TNM_data.RData")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/ABHD11_AS1_TNM分期堆积直方图.pdf",width = 4,height = 5)
ggplot(ESCC_2_STR_exp_clin_prop, aes(x = new_group, y = prop, fill = TNM_stage)) +
  geom_bar(stat = "identity", position = "fill", color = "white", size = 0.5) +  # 加边框线
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("I/II" = "#B5A8CA", "III/IV" = "#826BA2")) +
  labs(x = "new group", y = "Proportion", fill = "TNM Stage",title = "TNM Stage in two Groups") +
  #annotate("text", x = 1.5, y = 1.05, label = "p-value = 0.1048", size = 5, fontface = "italic") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),  # 仅保留坐标轴线
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14)
  )
dev.off()

## J
library(ggplot2)
met_ESCC_2_STR_exp_clin  <- read.csv(file = "./ABHD11_AS1淋巴转移数据.csv")
met_ESCC_2_STR_exp_clin$Lymphatic_metastasis_new <- as.character(met_ESCC_2_STR_exp_clin$Lymphatic_metastasis_new)
met_ESCC_2_STR_exp_clin$new_group <- factor(met_ESCC_2_STR_exp_clin$new_group,levels = c("Equal","High"))
####卡方检验计算p值
met_tab <- xtabs(n ~ new_group + Lymphatic_metastasis_new,
                 data = met_ESCC_2_STR_exp_clin)
met_p_value <- chisq.test(met_tab, correct = FALSE)$p.value
met_p_label <- paste0("P = ", signif(met_p_value, 3))
met_p_label

pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/ABHD11_AS1_淋巴结转移堆积直方图.pdf",width = 4,height = 5)
ggplot(met_ESCC_2_STR_exp_clin, aes(x = new_group, y = prop, fill = Lymphatic_metastasis_new)) +
  geom_bar(stat = "identity", position = "fill", color = "white", size = 0.5) +  # 加边框线
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = c("0" = "#deebf7", "1" = "#487DB2")) +  # 自定义渐变色
  labs(x = "new group", y = "Proportion", fill = "Lymphatic_metastasis_new",title = "Lymphatic metastasis in two Groups") +
  #annotate("text", x = 1.5, y = 1.05, label = "p-value = ", size = 5, fontface = "italic") +
  theme_minimal()+  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),  # 仅保留坐标轴线
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14)
  )
dev.off()

#Figure S2
## A
library(ggplot2)
library(dplyr)
all_data <- read.csv(file = "./figS2_a_纳入STR后遗传力差值数据.csv")
all_data$phenotype <- factor(
  all_data$phenotype,
  levels = c("APA", "Expression", "Splicing", "DNA methylation")
)

mean_df <- all_data %>%
  group_by(phenotype) %>%
  summarise(
    mean_h2 = mean(h2_chazhi, na.rm = TRUE),
    .groups = "drop"
  )

pdf(file = "/share/pub/changxj/sunwy/all_Heritability/遗传力差值_密度分布.pdf",width = 5,height = 4)
ggplot(all_data, aes(x = h2_chazhi, fill = phenotype, color = phenotype)) +
  geom_density(alpha = 0.3, size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
  geom_vline(
    data = mean_df,
    aes(xintercept = mean_h2, color = phenotype),
    linetype = "solid",
    size = 0.8
  ) +
  # 标注均值
  geom_text(
    data = mean_df,
    aes(
      x = mean_h2,
      y = Inf,
      label = paste0("mean value = ", round(mean_h2, 3)),
      color = phenotype
    ),
    inherit.aes = FALSE,
    vjust = 1.5,
    hjust = -0.05,
    size = 3.5
  ) +
  
  xlim(0, 0.15) +
  
  scale_fill_manual(values = c(
    "APA" = "#E58B87",
    "Expression" = "#89A3CD",
    "Splicing" = "#D28AC4",
    "DNA methylation" = "#E79F57"
  )) +
  
  scale_color_manual(values = c(
    "APA" = "#E58B87",
    "Expression" = "#89A3CD",
    "Splicing" = "#D28AC4",
    "DNA methylation" = "#E79F57"
  )) +
  
  labs(
    x = expression(Delta*h^2~"(Total - SNP joint)"),
    y = "Density"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    strip.background = element_blank(),    
    strip.text = element_text(size = 12, colour = "black"),  
    axis.text = element_text(size = 13, colour = "black"),
    legend.position = "none"             
  ) +
  
  facet_wrap(~phenotype, ncol = 1, strip.position = "right")
dev.off()


## B
library(ggridges)
library(ggplot2)
library(dplyr)
SNP_mean_h2_phenotype <- read.csv(file = "./figS2_b_SNP遗传力数据.csv")
SNP_mean_h2_phenotype <- SNP_four_phenotype_h2 %>%
  group_by(phenotype) %>%
  summarise(mean_h2 = mean(h2))

SNP_four_phenotype_h2$phenotype <- factor(
  SNP_four_phenotype_h2$phenotype,
  levels = c("APA", "Gene", "Splicing", "CpG")
)

SNP_mean_h2_phenotype$phenotype <- factor(
  SNP_mean_h2_phenotype$phenotype,
  levels = c("APA", "Gene", "Splicing", "CpG")
)
pdf(file = "/share/pub/changxj/sunwy/STR_Heritability/SNP对四种分子表型的遗传力.pdf",width = 5,height = 3)
ggplot(SNP_four_phenotype_h2, aes(x = h2, fill = phenotype, color = phenotype)) +
  geom_density(alpha = 0.6, adjust = 1.5, bw = 0.1, size = 0.8) +
  geom_vline(
    data = SNP_mean_h2_phenotype,
    aes(xintercept = mean_h2, color = phenotype),
    linetype = "dashed",
    linewidth = 0.8,
    show.legend = FALSE
  )+
  
  # ✅ 加均值文字
  geom_text(
    data = SNP_mean_h2_phenotype,
    aes(x = 0.3, y = 3, label = round(mean_h2, 3)),
    color = "black",
    vjust = 2,
    size = 5
  ) +
  facet_grid(phenotype ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = c(
    "APA" = "#f1c0be",
    "Gene" = "#acbedc",
    "Splicing" = "#DEAAD4",
    "CpG" = "#F1C89A"
  )) +
  scale_color_manual(values = c(
    "APA" = "#f1c0be",
    "Gene" = "#acbedc",
    "Splicing" = "#DEAAD4",
    "CpG" = "#F1C89A"
  ))+
  labs(
    x = "SNP-associated heritability of four molecular phenotypes",
    y = NULL
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    text = element_text(color = "black"),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    legend.position = "none",
    panel.spacing = unit(0.2, "lines")
  )
dev.off()


#Figure S3
##A
library(dplyr)
library(ggplot2)
all_hmm_or <- read.csv(file = "./figS3_a_染色质状态富集数据.csv")
all_hmm_or$Element <- factor(all_hmm_or$Element, levels = state_order)
all_hmm_or$phenotype <- factor(all_hmm_or$phenotype,levels = rev(c("aSTR","eSTR","sSTR","mSTR")))
pd <- position_dodge(width = 0.6)
pdf(file = "E:/swy/EnsembleTR/四种分子表型QTL显著的STR对应染色质状态富集结果图_CI_log.pdf",width = 7,height = 7)
ggplot(all_hmm_or,
       aes(x = logOR,
           y = Element,
           color = phenotype)) +
  
  geom_point(
    position = pd,
    size = 3
  ) +
  
  geom_errorbarh(
    aes(xmin = log_CI_low, xmax = log_CI_high),
    position = pd,
    height = 0.2
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  
  scale_color_manual(values = c(
    aSTR = "#E58B87",
    eSTR = "#89A3CD",
    sSTR = "#D28AC4",
    mSTR = "#E79F57"
  )) +
  
  labs(
    x = "log(OR)",
    y = "",
    color = "Phenotype"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text = element_text(size = 13,colour = "black"),
    axis.title = element_text(size = 14),
    legend.text = element_text(size = 13)
  )
dev.off()

## B
phenotype_site <- read.csv(file = "./figS3_b_每个分子表型对应STR画图数据.csv")
phenotype_site$type <- factor(phenotype_site$type,levels = c("DNA Methylation","Splicing","Expression","APA"))
library(ggplot2)
pdf("E:/swy/EnsembleTR/每个位点对应平均STR数目柱状图.pdf",width = 6 ,height = 4)
ggplot(phenotype_site, aes(x = type, y = phenotyoe_mean, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +  # 横向条形图
  scale_fill_manual(values = c("APA" = "#f1c0be", "Expression" = "#acbedc", "Splicing" = "#EBCCE4","DNA Methylation" = "#F1C89A")) +
  labs(
    x = NULL,
    y = "Mean number of associated sites",
    title = NULL
  ) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()

## C
xSTR_site <- read.csv(file = "./figS3_c_xSTR对应分子表型数据.csv")
xSTR_site$type <- factor(xSTR_site$type,levels = c("mSTR","sSTR","eSTR","aSTR"))
library(ggplot2)
pdf("E:/swy/EnsembleTR/每个STR对应平均表型数目柱状图.pdf",width = 6 ,height = 4)
ggplot(xSTR_site, aes(x = type, y = phenotyoe_mean, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +  # 横向条形图
  scale_fill_manual(values = c("aSTR" = "#f1c0be", "eSTR" = "#acbedc", "sSTR" = "#EBCCE4","mSTR" = "#F1C89A")) +
  labs(
    x = NULL,
    y = "Mean number of associated sites",
    title = NULL
  ) +
  theme_bw() + 
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14,colour = "black"),
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.title.x = element_text(size = 15,colour = "black"),
    axis.title.y = element_text(size = 15,colour = "black"),
    axis.line = element_line(size = 0.5),  # 仅保留坐标轴线
    # 加上刻度线：
    #axis.ticks = element_line(color = "black", size = 0.5),
    #axis.ticks.length = unit(5, "pt"),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14,colour = "black"),
    panel.grid.major = element_blank(),  # 删除大网格线
    panel.grid.minor = element_blank(),  # 删除小网格线
    panel.background = element_blank()   # 删除背景色
  ) 
dev.off()

## D
library(ggtern)
plot_data_CpG <- read.csv(file = "./figS3_d_CpG三元关系对展示数据.csv")
genes_to_label <- c("ENSG00000213689_chr3:48385339-48385370-2-TG_3:48444085","ENSG00000250072_chr5:149074232-149074245-2-TG_5:149055268")
plot_data_CpG$Model <- factor(plot_data_CpG$Model, levels = c("complete", "co-local", "partial", "other"))
pdf(file = "E:/swy/EnsembleTR/bmeditaR/Three_CpG_figure.pdf",width = 6,height = 6)
ggtern(data = plot_data_CpG,
       aes(x = `complete med`, y = `co-local`, z = `partial med`,
           color = Model, fill = Model, size = log10_qvalue)) +
  # 加粗三角形边框
  tern_limits(T = 1, L = 1, R = 1) +
  theme_showarrows() +
  theme_showsecondary() +
  theme_rgbw() +
  theme(
    tern.axis.line.T = element_line(size = 1, colour = "black"),
    tern.axis.line.L = element_line(size = 1, colour = "black"),
    tern.axis.line.R = element_line(size = 1, colour = "black"),
    tern.axis.text.T = element_blank(),  # 去坐标刻度
    tern.axis.text.L = element_blank(),
    tern.axis.text.R = element_blank()
  ) +
  # 实心圆圈 + 边框粗
  geom_point(shape = 21, stroke = 0.3, alpha = 0.5) +
  geom_point(
    data = subset(plot_data_CpG, id %in% genes_to_label),
    aes(x = `complete med`, y = `co-local`, z = `partial med`),
    shape = 21,
    size = 3,
    fill = "black",
    color = "black",
    stroke = 1
  ) +
  geom_label(
    data = subset(plot_data_CpG, id %in% genes_to_label),
    aes(label = id),
    size = 3,
    label.size = 0.2,
    color = "black",
    fill = "white"
  )+
  # 手动颜色设置
  scale_color_manual(
    values = c(
      "complete" = "#DB498E",
      "co-local" = "#A29EC1",
      "partial" = "#609ACA",
      "other" = "gray"
    ),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "complete" = "#EC6A90",
      "co-local" = "#A29EC1",
      "partial" = "#609ACA",
      "other" = "#EBEBEB"
    ),
    drop = TRUE
  ) +
  
  # 图例美化
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12)
  ) +
  guides(
    size = guide_legend(title = "-log10(q.value)")
  )
dev.off()

## E
TREX1_tr_CpG <- read.csv(file = "./figS3_e_TREX1_tr_CpG散点图数据.csv")
TREX1_tr_CpG$color <- as.factor(TREX1_tr_CpG$geno_new)
library(ggpubr)
library(ggplot2)
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/NEW_TREX1_points.pdf",width = 6,height = 4)
ggplot(TREX1_tr_CpG, aes(x = DNA_methylation, y = Gene_expression)) +
  geom_point(aes(colour = color)) +
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+theme_classic()+
  scale_color_manual(values = c("#009E73","#FFDB6D", "#4E84C4", "#CC79A7", "#D55E00", "#D55E00", "#D55E00"))
dev.off()



#Figure S4
## A
library(ggtern)
library(dplyr)
library(ggplot2)
plot_data <- read.csv(file = "./figS4_a_Splicing三元关系对展示数据.csv")
plot_data$Model <- factor(plot_data$Model, levels = c("complete", "co-local", "partial", "other"))
pdf(file = "E:/swy/EnsembleTR/bmeditaR/Three_splicing_figure1.pdf",width = 6,height = 6)
ggtern(data = plot_data,
       aes(x = `complete med`, y = `co-local`, z = `partial med`,
           color = Model, fill = Model, size = log10_qvalue)) +
  # 加粗三角形边框
  tern_limits(T = 1, L = 1, R = 1) +
  theme_showarrows() +
  theme_showsecondary() +
  theme_rgbw() +
  theme(
    tern.axis.line.T = element_line(size = 1, colour = "black"),
    tern.axis.line.L = element_line(size = 1, colour = "black"),
    tern.axis.line.R = element_line(size = 1, colour = "black"),
    tern.axis.text.T = element_blank(),  # 去坐标刻度
    tern.axis.text.L = element_blank(),
    tern.axis.text.R = element_blank()
  ) +
  # 实心圆圈 + 边框粗
  geom_point(shape = 21, stroke = 0.3, alpha = 0.5) +
  geom_point(
    data = subset(plot_data, id %in% genes_to_label),
    aes(x = `complete med`, y = `co-local`, z = `partial med`),
    shape = 21,
    size = 3,
    fill = "black",
    color = "black",
    stroke = 1
  ) +
  geom_label(
    data = subset(plot_data, id %in% genes_to_label),
    aes(label = id),
    size = 3,
    label.size = 0.2,
    color = "black",
    fill = "white"
  )+
  # 手动颜色设置
  scale_color_manual(
    values = c(
      "complete" = "#DB498E",
      "co-local" = "#A29EC1",
      "partial" = "#609ACA",
      "other" = "gray"
    ),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "complete" = "#EC6A90",
      "co-local" = "#A29EC1",
      "partial" = "#609ACA",
      "other" = "#EBEBEB"
    ),
    drop = TRUE
  ) +
  
  # 图例美化
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12)
  ) +
  guides(
    size = guide_legend(title = "-log10(q.value)")
  )
dev.off()


## B
MUC20P1_tr_Splicing <- read.csv(file = "./figS4_bMUC20P1_points画图数据.csv")
library(RColorBrewer)
n <- length(unique(MUC20P1_tr_Splicing$color))  # 替换成你的分组列
# 自动生成足够的颜色
cell_type_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
# 建立多元线性回归模型
library(ggpubr)
library(ggplot2)
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/MUC20P1_points.pdf",width = 6,height = 4)
ggplot(MUC20P1_tr_Splicing, aes(x = Splicing_level, y = Gene_expression)) +
  geom_point(aes(colour = color)) +
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+theme_classic()+
  scale_color_manual(values = cell_type_cols)
dev.off()


## C
eqtl_MUC20P1 <- read.csv(file = "./figS4_c_eQTL_MUC20P1.csv")
library(RColorBrewer)
length(unique(eqtl_MUC20P1$geno_new))
n <- length(unique(eqtl_MUC20P1$geno_new))  # 替换成你的分组列
# 自动生成足够的颜色
cell_type_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
eqtl_MUC20P1$geno_new <- factor(as.character(eqtl_MUC20P1$geno_new), levels = sort(as.numeric(unique(eqtl_MUC20P1$geno_new))))
library(ggplot2)
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_MUC20P1_bmeditaR.pdf",width =4,height = 5)
ggplot(data = eqtl_MUC20P1, aes(x = geno_new, y = exp, fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of MUC20P1 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## D
splicing_MUC20P1 <- read.csv(file = "./figS4_d_sQTL_MUC20P1.csv")
library(ggplot2)
library(RColorBrewer)
n <- length(unique(splicing_MUC20P1$geno_new)) # 替换分组列
cell_type_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)
splicing_MUC20P1$geno_new <- factor(as.character(splicing_MUC20P1$geno_new), levels = sort(as.numeric(unique(splicing_MUC20P1$geno_new))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/splicing_MUC20P1_bmeditaR.pdf",width = 4,height = 5)
ggplot(data = splicing_MUC20P1, aes(x = geno_new, y = exp, fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "Splicing")+  # 修改坐标轴名称
  ggtitle("TR-splicing of HLA-DRB5 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## E
MUC20P1_tr_Splicing_clean <- read.csv(file = "./figS4_e_resid_eQTL_MUC20P1.csv")
cor.test(MUC20P1_tr_Splicing_clean$resid_MUC20P1, as.numeric(MUC20P1_tr_Splicing_clean$TR_MUC20P1))
library(ggplot2)
library(RColorBrewer)
n <- length(unique(MUC20P1_tr_Splicing$color))  # 替换成你的分组列
cell_type_cols <- colorRampPalette(brewer.pal(12, "Set3"))(n)

MUC20P1_tr_Splicing_clean$geno_new <- factor(as.character(MUC20P1_tr_Splicing_clean$geno_new), levels = sort(as.numeric(unique(MUC20P1_tr_Splicing_clean$geno_new))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/resid_splicing_MUC20P1_bmeditaR.pdf",width = 4,height = 5)
ggplot(MUC20P1_tr_Splicing_clean, aes(x = geno_new, y = resid_MUC20P1,fill = geno_new)) +
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "Residual Exp")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of MUC20P1 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))

dev.off()


## F
plot_data  <- read.csv(file = "./figS4_f_APA三元关系对画图数据.csv")
library(ggtern)
library(ggplot2)
genes_to_label <- c("ENSG00000198502_chr6:32557839-32557864-2-CA_XM_011514562.2|HLA-DRB5|chr6|-",
                    "ENSG00000214941_chr17:16199217-16199244-2-AC_NM_001042697.2|ZSWIM7|chr17|-")

plot_data$Model <- factor(plot_data$Model, levels = c("complete", "co-local", "partial", "other"))
pdf(file = "E:/swy/EnsembleTR/bmeditaR/Three_apa_figure.pdf",width = 6,height = 6)
ggtern(data = plot_data,
       aes(x = `complete med`, y = `co-local`, z = `partial med`,
           color = Model, fill = Model, size = log10_qvalue)) +
  # 加粗三角形边框
  tern_limits(T = 1, L = 1, R = 1) +
  theme_showarrows() +
  theme_showsecondary() +
  theme_rgbw() +
  theme(
    tern.axis.line.T = element_line(size = 1, colour = "black"),
    tern.axis.line.L = element_line(size = 1, colour = "black"),
    tern.axis.line.R = element_line(size = 1, colour = "black"),
    tern.axis.text.T = element_blank(),  # 去坐标刻度
    tern.axis.text.L = element_blank(),
    tern.axis.text.R = element_blank()
  ) +
  # 实心圆圈 + 边框粗
  geom_point(shape = 21, stroke = 0.3, alpha = 1) +
  geom_point(
    data = subset(plot_data, id %in% genes_to_label),
    aes(x = `complete med`, y = `co-local`, z = `partial med`),
    shape = 21,
    size = 3,
    fill = "black",
    color = "black",
    stroke = 1
  ) +
  geom_label(
    data = subset(plot_data_CpG, id %in% genes_to_label),
    aes(label = id),
    size = 3,
    label.size = 0.2,
    color = "black",
    fill = "white"
  )+
  # 手动颜色设置
  scale_color_manual(
    values = c(
      "complete" = "#DB498E",
      "co-local" = "#A29EC1",
      "partial" = "#609ACA",
      "other" = "gray"
    ),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "complete" = "#EC6A90",
      "co-local" = "#A29EC1",
      "partial" = "#609ACA",
      "other" = "#EBEBEB"
    ),
    drop = TRUE
  ) +
  
  # 图例美化
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12)
  ) +
  guides(
    size = guide_legend(title = "-log10(q.value)")
  )
dev.off()

## G
library(ggplot2)
HLA_DRB5_tr_apa <- read.csv(file = "./figS4_g_HLA_DRB5_points.csv")
HLA_DRB5_tr_apa$color <- as.factor(HLA_DRB5_tr_apa$geno_new)
library(ggpubr)
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/HLA_DRB5_points.pdf",width = 6,height = 4)
ggplot(HLA_DRB5_tr_apa, aes(x = apa_level, y = Gene_expression)) +
  geom_point(aes(colour = color)) +
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+theme_classic()+
  scale_color_manual(values = c("#FF7DA0","#FFA5A5", "#377EB8"))
dev.off()

## H
library(ggplot2)
eqtl_HLA_DRB5 <- read.csv(file = "./figS4_h_eqtl_HLA_DRB5.csv")
cell_type_cols <- c("#ED7E8B","#FFD19C","#B797C6")
eqtl_HLA_DRB5$geno_new <- factor(as.character(eqtl_HLA_DRB5$geno_new), levels = sort(as.character(unique(eqtl_HLA_DRB5$geno_new))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_HLA_DRB5_bmeditaR.pdf",width =4,height = 5)
ggplot(data = eqtl_HLA_DRB5, aes(x = geno_new, y = exp, fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of HLA-DRB5 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()

## I
library(ggplot2)
aqtl_HLA_DRB5 <- read.csv(file = "./figS4_i_aqtl_HLA_DRB5.csv")
length(unique(aqtl_HLA_DRB5$geno_new))
cell_type_cols <- c("#D65190","#731A73","#FF8C00","#699ECA","#FFDB6D", "#52854C", "#4E84C4","#009E73","#FFDB6D")
aqtl_HLA_DRB5$geno_new <- factor(as.character(aqtl_HLA_DRB5$geno_new), levels = sort(as.character(unique(aqtl_HLA_DRB5$geno_new))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/aqtl_HLA_DRB5_bmeditaR.pdf",width = 4,height = 5)
ggplot(data = aqtl_HLA_DRB5, aes(x = geno_new, y = exp, fill = geno_new)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "APA")+  # 修改坐标轴名称
  ggtitle("TR-aQTL of HLA-DRB5 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()

## J
HLA_DRB5_tr_apa_clean <- read.csv(file = "./figS4_j_Residual_eQTL_HLA_DRB5.csv")
cor.test(HLA_DRB5_tr_apa_clean$resid_HLA_DRB5, as.numeric(HLA_DRB5_tr_apa_clean$TR_HLA_DRB5))
library(ggplot2)
cell_type_cols <- c("#ED7E8B","#FFD19C","#B797C6")
HLA_DRB5_tr_apa_clean$geno_new <- factor(as.character(HLA_DRB5_tr_apa_clean$geno_new), levels = sort(as.character(unique(HLA_DRB5_tr_apa_clean$geno_new))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/resid_aqtl_HLA_DRB5_bmeditaR.pdf",width = 4,height = 5)
ggplot(HLA_DRB5_tr_apa_clean, aes(x = geno_new, y = resid_HLA_DRB5,fill = geno_new)) +
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "Residual Exp")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of HLA-DRB5 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))

dev.off()

#Figure S5
## A
library(ggplot2)
library(tidyr)
no_GO_bp_long <- read.csv(file = "./figS5_a_非显著QTL对应基因的富集气泡图数据.csv")
no_GO_bp_long$Term <- factor(no_GO_bp_long$Term, levels = rev(unique(no_GO_bp_long$Term)))
no_GO_bp_long$QTLtype <- factor(no_GO_bp_long$QTLtype,levels = c("aQTL","eQTL","sQTL","mQTL"))
pdf("E:/swy/EnsembleTR/NO_new_GO_bp_dotplot_plot.pdf",  width = 11, height = 7)
ggplot(no_GO_bp_long, aes(x = QTLtype, y = Term, color = log_p, size = log_combined)) +
  geom_point() +
  scale_color_gradient(low = "grey", high = "#E54924") +  # 颜色梯度
  scale_size_continuous(range = c(2, 8)) +  # 点的大小
  theme_minimal() +
  theme(
    axis.text.x = element_text(hjust = 0.5,size = 14),  # 旋转 X 轴文本
    axis.title.x = element_text(size = 15),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  # 添加黑色框线
  ) +
  labs(x = "QTL Type", y = "GO Biological Process", 
       color = "-log10(p)", size = "Combined.Score", 
       title = "GO Biological Processes across STR-nQTL")

dev.off()

## B
library(ggplot2)
library(patchwork)
results <- read.csv(file = "./figS5_b_四种分子表型对应基因突变OR值数据.csv")
p1 <-  ggplot(results, aes(x = QTL, y = OR)) +
  geom_col(fill = "#CEC1DC", alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = CI_low, ymax = CI_high),
                width = 0.2, size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "#ED7E8B") +
  geom_text(aes(label = signif, y = OR + 0.15), vjust = 0, size = 6) +
  theme_minimal(base_size = 14) +
  labs(y = "OR (Fisher test)", x = "") +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))

p1

## C
library(ggplot2)
#eGene绘制棒棒糖图
top_10_eGene_long <- read.csv(file = "./figS5_c_eGene_cnv.csv")
top_10_eGene_long$symbol <- factor(top_10_eGene_long$symbol,
                                   levels = unique(top_10_eGene_long$symbol))
p1 <- ggplot(top_10_eGene_long, aes(x = symbol, y = Frequency, color = CNV_type)) +
  geom_segment(aes(xend = symbol, y = 0, yend = Frequency), size = 0.8) + # 棒
  geom_point(size = 3) +  # 圆点
  scale_color_manual(values = c("Gain" = "#E41A1C", "Loss" = "#377EB8")) + # 红蓝
  labs(x = "eQTL top10 Gene", y = "Frequency", color = "CNV Type") +theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,size=14,angle = 60),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
p1
#aGene绘制棒棒糖图
top_10_aGene_long <- read.csv(file = "./figS5_c_aGene_cnv.csv")
top_10_aGene_long$symbol <- factor(top_10_aGene_long$symbol,
                                   levels = unique(top_10_aGene_long$symbol))
p2 <- ggplot(top_10_aGene_long, aes(x = symbol, y = Frequency, color = CNV_type)) +
  geom_segment(aes(xend = symbol, y = 0, yend = Frequency), size = 0.8) + # 棒
  geom_point(size = 3) +  # 圆点
  scale_color_manual(values = c("Gain" = "#E41A1C", "Loss" = "#377EB8")) + # 红蓝
  labs(x = "aQTL top10 Gene", y = "Frequency", color = "CNV Type") +theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,size=14,angle = 60),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
p2
#sGene绘制棒棒糖图
top_10_sGene_long <- read.csv(file = "./figS5_c_sGene_cnv.csv")
top_10_sGene_long$symbol <- factor(top_10_sGene_long$symbol,
                                   levels = unique(top_10_sGene_long$symbol))
p3 <- ggplot(top_10_sGene_long, aes(x = symbol, y = Frequency, color = CNV_type)) +
  geom_segment(aes(xend = symbol, y = 0, yend = Frequency), size = 0.8) + # 棒
  geom_point(size = 3) +  # 圆点
  scale_color_manual(values = c("Gain" = "#E41A1C", "Loss" = "#377EB8")) + # 红蓝
  labs(x = "sQTL top10 Gene", y = "Frequency", color = "CNV Type") +theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,size=14,angle = 60),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
p3
#mGene绘制棒棒糖图
top_10_mGene_long <- read.csv(file = "./figS5_c_mGene_cnv.csv")
top_10_mGene_long$symbol <- factor(top_10_mGene_long$symbol,
                                   levels = unique(top_10_mGene_long$symbol))
p4 <- ggplot(top_10_mGene_long, aes(x = symbol, y = Frequency, color = CNV_type)) +
  geom_segment(aes(xend = symbol, y = 0, yend = Frequency), size = 0.8) + # 棒
  geom_point(size = 3) +  # 圆点
  scale_color_manual(values = c("Gain" = "#E41A1C", "Loss" = "#377EB8")) + # 红蓝
  labs(x = "mQTL top10 Gene", y = "Frequency", color = "CNV Type") +theme_bw() +
  theme(axis.text.x = element_text(hjust = 1,size=14,angle = 60),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
p4



#Figure S6
## B
library(ggplot2)
dat <- read.csv(file = "figS6_b_sQTL_LPAR1.csv")
cell_type_cols <- c("#bcd1d2","#CEC1DC","#EDEAC9")
dat$geno <- factor(as.numeric(dat$geno), levels = sort(as.numeric(unique(dat$geno))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/sqtl_LPAR1_susieR.pdf",width =6,height = 5)
ggplot(data = dat, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  scale_y_log10()+
  theme_classic() +
  labs(x = "TR Length", y = "Normalized \n splicing Ratio")+  # 修改坐标轴名称
  ggtitle("TR-sqtl of LPAR1")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## C
library(ggplot2)
LPAR1_susieR <- read.csv(file = "./figS6_c_down_susieR_LPAR1.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/sSTR_LPAR1_susieR.pdf",width =6,height =4)
ggplot(LPAR1_susieR, aes(x = variant_id, y = pip, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of LPAR1",
       x = "Genomic Position", y = "Casual Probability (PIP)") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()
LPAR1_susieR <- read.csv(file = "./figS6_c_up_susieR_p_LPAR1.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/sSTR_LPAR1_susieR_log10p.pdf",width =6,height =4)
ggplot(LPAR1_susieR, aes(x = variant_id, y = log10p, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  #scale_x_continuous(limits = c(x_min, x_max), labels = scales::comma) +
  labs(title = "Fine-mapping of LPAR1",
       x = "Genomic Position", y = "-log10(P)") +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "none",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## E
library(dplyr)
eqtl_KRTAP4 <- read.csv("./figS6_e_mQTL_KRTAP4.csv") 
cell_type_cols <- c("#bcd1d2","#CEC1DC","#EDEAC9")
eqtl_KRTAP4$geno <- factor(as.character(eqtl_KRTAP4$geno), levels = sort(as.character(unique(eqtl_KRTAP4$geno))))
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eqtl_KRTAP4_susieR.pdf",width =4,height = 5)
ggplot(data = eqtl_KRTAP4, aes(x = geno, y = exp, fill = geno)) + 
  geom_boxplot() + 
  scale_fill_manual(values = cell_type_cols) + 
  theme_classic() +
  labs(x = "TR Length", y = "log2(Exp+1)")+  # 修改坐标轴名称
  ggtitle("TR-eQTL of KRTAP4 in ESCC")+theme_bw() +
  theme(axis.text.x = element_text(hjust = 0.5,size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        axis.line = element_line(size=1),
        plot.title = element_text(hjust = 0.5,size = 15),
        legend.position = "right",
        legend.title = element_blank(),
        legend.background = element_rect(fill = 'transparent'),
        legend.text = element_text(size = 14))
dev.off()


## F
KRTAP4_susieR <- read.csv(file = "./figS6_f_down_susieR_KRTAP4.csv")
library(ggplot2)
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eSTR_KRTAP4_susieR.pdf",width =4,height =3)
ggplot(KRTAP4_susieR, aes(x = variant_id, y = pip, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of KRTAP4",
       x = "Genomic Position", y = "Casual Probability (PIP)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()   # 去除次网格线
  )
dev.off()
KRTAP4_susieR <- read.csv(file = "./figS6_f_up_susieR_p_KRTAP4.csv")
pdf("/share/pub/sunwy/sunwy/EnsembleTR_result/eSTR_KRTAP4_susieR_log10p.pdf",width =4,height =3)
ggplot(KRTAP4_susieR, aes(x = variant_id, y = log10p, color = type)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("SNP" = "gray", "STR" = "darkorange")) +
  labs(title = "Fine-mapping of KRTAP4",
       x = "Genomic Position", y = "-log10(P)") +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size=8,colour = "black"),
    axis.title.x = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 9),
    legend.position = "none",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 9),
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()   # 去除次网格线
  )
dev.off()


## G
library(ggridges)
library(dplyr)
library(ggplot2)
all_distence <- read.csv(file = "./figS6_g_QTL显著分子表型距离STR画图数据.csv")
all_distence$group <- factor(all_distence$group,levels = c("Exp","APA","Splicing","Methylation"))
pdf("E:/swy/EnsembleTR/four_distence_to_STR.pdf",width = 8,height = 4)
ggplot(all_distence, aes(x = distence, y = group,
                         fill = group, color = group)) +
  geom_density_ridges(rel_min_height = 0.001, alpha = 0.5) +
  stat_density_ridges(aes(color = group),
                      quantile_lines = TRUE,
                      quantiles = 2,
                      linewidth = 1,
                      alpha = 0.5) +
  scale_fill_manual(values = c( "#933F42", "#CFC5BB", "#6386A5","#757877")) +
  scale_color_manual(values = c( "#933F42", "#CFC5BB", "#6386A5","#757877")) +
  scale_y_discrete("") + theme_classic() + 
  theme_minimal() +  # 使用无边框的主题
  theme(
    axis.text.x = element_text(hjust = 0.5, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.title.y = element_text(size = 15),
    axis.line = element_line(size = 1),  # 仅保留坐标轴线
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.background = element_rect(fill = 'transparent'),
    legend.text = element_text(size = 14)
  ) + labs(x = "Distance", y = "",title = "Distance to STRs")

dev.off()



#Figure S7
## A
library(dplyr)
library(ggplot2)
all_genome_or <- read.csv(file = "./figS7_精细定位之后的基因组元件富集数据.csv")
state_order <- rev(c(
  "3' UTR","Promoter","Exon","Downstream",
  "5' UTR","Intron","Distal Intergenic"
))

all_genome_or$Element <- factor(all_genome_or$Element, levels = state_order)
plot_or <- function(df, color_use, title_text){
  ggplot(df,
         aes(x = logOR,
             y = Element)) +
    geom_point(
      color = color_use,
      size = 3
    ) +  
    geom_errorbarh(
      aes(xmin = log_CI_low, xmax = log_CI_high),
      height = 0.2,
      color = color_use
    ) +    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +    
    labs(
      title = paste0(title_text," annotation"),
      x = "Odds Ratio",
      y = ""
    ) +
    theme_classic() + 
    theme(
      axis.text = element_text(size = 13, colour = "black"),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 15, hjust = 0.5)
    )
}
pdf(file = "E:/swy/EnsembleTR/精细定位四种分子表型STR对应基因组元件结果图_CI_log.pdf",width = 4,height = 3)
plot_or(all_genome_or[all_genome_or$phenotype %in% "aSTR",], "#E58B87","aSTR")
plot_or(all_genome_or[all_genome_or$phenotype %in% "eSTR",], "#89A3CD","eSTR")
plot_or(all_genome_or[all_genome_or$phenotype %in% "sSTR",], "#D28AC4","sSTR")
plot_or(all_genome_or[all_genome_or$phenotype %in% "mSTR",], "#E79F57","mSTR")
dev.off()


## B
library(dplyr)
library(ggplot2)
all_motif_or <- read.csv(file = "./figS7_b_精细定位之后motif_length画图数据.csv")
all_motif_or$phenotype <- factor(
  all_motif_or$phenotype,
  levels = c("aSTR","eSTR","sSTR","mSTR")
)

lm_results <- all_motif_or %>%
  group_by(phenotype) %>%
  do({
    fit <- lm(OR ~ motif_length, data = .)
    data.frame(
      R2 = summary(fit)$r.squared,
      p = summary(fit)$coefficients[2,4]
    )
  })


lm_results <- lm_results %>%
  mutate(
    label = paste0("R = ", round(R2,3), 
                   ",P = ", signif(p,3))
  )
pdf("E:/swy/EnsembleTR/fine_mapped_四种分子表型STR的motif长度分布特征.pdf",width = 7,height = 5)
ggplot(all_motif_or, aes(x = motif_length, y = OR, color = phenotype)) +
  
  geom_point(size = 3) +
  
  annotate(geom = 'rect',xmin = -Inf, xmax = Inf,ymin = 1, ymax = 2,fill = "#F1C89A",alpha = 0.15 ) +
  
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    color = "#F1C89A",
    size = 0.8
  ) +
  
  geom_smooth(
    method = lm,
    se = FALSE
  ) +
  
  scale_color_manual(values = c(
    "#f1c0be",  # aSTR
    "#acbedc",  # eSTR
    "#DEAAD4",  # mSTR
    "#E79F57"   # sSTR
  )) +
  
  scale_y_continuous(
    breaks = seq(0,2,0.5),
    expand = c(0,0.01),
    limits = c(0,2)
  ) +
  
  labs(
    x = "Motif length",
    y = "Odds Ratio",
    title = "Motif size features of fine-mapped STRs"
  ) +
  theme(
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 15, colour = "black"),
    axis.title.y = element_text(size = 15, colour = "black"),
    axis.line = element_line(size = 1),
    plot.title = element_text(hjust = 0.5, size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )
dev.off()


## C
library(dplyr)
library(ggplot2)
fine_mapped_aSTR_motif_or_results_top10 <- read.csv(file = "./figS7_c_精细定位后aSTR_motif_unit的数据.csv")
fine_mapped_eSTR_motif_or_results_top10 <- read.csv(file = "./figS7_c_精细定位后eSTR_motif_unit的数据.csv")
fine_mapped_sSTR_motif_or_results_top10 <- read.csv(file = "./figS7_c_精细定位后sSTR_motif_unit的数据.csv")
fine_mapped_mSTR_motif_or_results_top10 <- read.csv(file = "./figS7_c_精细定位后mSTR_motif_unit的数据.csv")
plot_or <- function(df, color_use, title_text){
  df <- df[order(df$OR, decreasing = T), ]
  ggplot(df, aes(x = Motif, y = OR, fill = sig)) + 
    geom_bar(stat = "identity") + 
    annotate(geom = 'rect', xmin = -Inf, xmax = Inf, ymin = 1, ymax = 1,
             fill = "#e66c2a", alpha = 0.2, color = '#ffd2ba') +  
    geom_hline(yintercept = 1, linetype = "dashed", color = "#e66c2a", size = 0.8) +
    # 自定义显著性颜色
    scale_fill_manual(values = c("Significant" = color_use, "Not significant" = "grey80")) +
    
    labs(title = title_text, x = "Motif unit", y = "Odds Ratio") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(hjust = 0.8, size = 14, angle = 30),
      axis.text.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15),
      axis.line = element_line(size = 1),
      # 加上刻度线：
      axis.ticks = element_line(color = "black", size = 1),
      axis.ticks.length = unit(3, "pt"),
      plot.title = element_text(hjust = 0.5, size = 15),
      legend.title = element_blank(),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
}
pdf(file = "E:/swy/EnsembleTR/精细定位四种分子表型STR对应motif_unit.pdf",width = 4,height = 3)
plot_or(fine_mapped_aSTR_motif_or_results_top10, "#E58B87","aSTR Motif unit")
plot_or(fine_mapped_eSTR_motif_or_results_top10, "#89A3CD","eSTR Motif unit")
plot_or(fine_mapped_sSTR_motif_or_results_top10, "#D28AC4","sSTR Motif unit")
plot_or(fine_mapped_mSTR_motif_or_results_top10, "#E79F57","mSTR Motif unit")
dev.off()


## D-G
library(dplyr)
library(ggplot2)
all_hmm_or <- read.csv(file = "./figS7_d_to_g精细定位后染色质状态富集CI取log的数据.csv")
all_hmm_or$Element <- factor(all_hmm_or$Element, levels = state_order)
all_hmm_or$phenotype <- factor(all_hmm_or$phenotype, levels = c("aSTR","eSTR","sSTR","mSTR"))

# 每个 phenotype 的颜色
phenotype_colors <- c(
  aSTR = "#E58B87",
  eSTR = "#89A3CD",
  sSTR = "#D28AC4",
  mSTR = "#E79F57"
)

# 循环分别作图
for (ph in levels(all_hmm_or$phenotype)) {
  
  plot_data <- all_hmm_or %>% filter(phenotype == ph)
  
  p <- ggplot(plot_data, aes(x = logOR, y = Element)) +
    geom_point(
      color = phenotype_colors[ph],
      size = 3
    ) +
    geom_errorbarh(
      aes(xmin = log_CI_low, xmax = log_CI_high),
      color = phenotype_colors[ph],
      height = 0.2
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    labs(
      title = ph,
      x = "Odds Ratio",
      y = ""
    ) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 13, colour = "black"),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
    )
  
  ggsave(
    filename = paste0("E:/swy/EnsembleTR/susieR/", ph, "_chromatin_state_enrichment_CI.pdf"),
    plot = p,
    width = 6,
    height = 7
  )
}


