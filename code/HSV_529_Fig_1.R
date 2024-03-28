
# code/HSV_529_Fig_1.R

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())

# INPUTS
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript'
} else {
  stop("set repo loc and repo manually")
}

filename_f1bc = file.path(repo_loc, 'data/tcellcounts.csv')
filename_f1ghi = file.path(repo_loc, 'data/fig1_ghi.csv')

person_color <-  c("1" = "#8d03a9", "2" = "#FFAD4A", "3" = "#0263be", 
                   "4" = "#FFE900", "5" = "#f93f61", "6" = "#6a4ee3", 
                   "7" = "#2e9a04", "8" = "#c8c89e", "9" = "#94ee78")

## OUTPUTS
figure_1b_filename = file.path(repo_loc, 'figures/fig_1b.pdf')
figure_1c_filename = file.path(repo_loc, 'figures/fig_1c.pdf')
figure_1g_filename = file.path(repo_loc, 'figures/fig_1g.pdf')
figure_1h_filename = file.path(repo_loc, 'figures/fig_1h.pdf')
figure_1i_filename = file.path(repo_loc, 'figures/fig_1i.pdf')

## Figure 1b, c
## prep 
tcells = read.csv(filename_f1bc, header = T)
colnames(tcells)
tcells = select(tcells, id, cd, control0, lesion, d0, d10, d30, d40, d180, d190, control10, control40, control190)
graphable = pivot_longer(tcells, cols = c(3:13), names_to = "day")
colnames(graphable) = c('id', 'cd', 'day', 'cells')
graphable$t = rep(c("arm", "lesion", 0, 10, 30, 40, 180, 190, 10, 40, 190), 18)
#graphable$t_plot = rep(c("a", "b", "c", "d", "e", "f", "g", "h",
#                         "i", "j", "k"), 18)
graphable$t_plot = rep(c(1:11), 18)
graph.cd4 = filter(graphable, cd == 4)
graph.cd8 = filter(graphable, cd == 8)


my_comparisons_paired <- list(c(1, 3), c(1, 9), c(1, 10), 
                              c(3, 4), c(5, 6), c(4, 9), c(6, 10))
my_comparisons_unpaired <- list(c(1, 11), c(7, 8), c(8, 11))
table(is.na(graph.cd4$cells), graph.cd4$day)

## CD4
plot <- graph.cd4 %>%
  select(t_plot, cells, id) %>%
  drop_na(cells) %>%
  ggplot(aes(x=t_plot, y=cells, group = t_plot), theme_bw()) + 
  geom_boxplot(color="black", fill="#bababa", outlier.shape=NA, width = 0.75, na.rm=TRUE,
               coef = 1.5) + 
  geom_point(aes(fill=factor(id)), size=3, shape=21, color='black') +
  geom_vline(xintercept=2.5, linetype = "dotted") + 
  geom_vline(xintercept=8.5, linetype = "dotted") +
  scale_x_continuous(breaks = c(1:11), labels=c("Arm", "Lesion", '0', '10', '30\n    Lesion Area Skin', '40', '180', '190', "10", "40 \n Arm", "190")) +
  scale_color_manual(values=person_color, guide = 'none') +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle(CD4^{"+"}~T~cells) + labs(y=Cells~per~mm^{2}, x="Biopsy time points (day)") + 
  theme_bw() + 
  theme(axis.text=element_text(size=10, hjust = 0.5),
        axis.title=element_text(size=10, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) + 
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, label.y = c(340, 360, 380, 300, 300, 315, 330),
                            comparisons = my_comparisons_paired, paired = TRUE, method = "wilcox", 
                            size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE) + 
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, label.y = c(400, 300, 345),
                             comparisons = my_comparisons_unpaired, paired = FALSE, method = "wilcox", 
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE)

pdf(figure_1b_filename, width = 6, height = 4)
plot
dev.off()
## note that you have to tweak the alignment of the text a little bit for visibility  

## CD8 
plot <- graph.cd8 %>%
  select(t_plot, cells, id) %>%
  drop_na(cells) %>%
  ggplot(aes(x=t_plot, y=cells, group = t_plot), theme_bw()) + 
  geom_boxplot(color="black", fill="#bababa", outlier.shape=NA, width = 0.75, na.rm=TRUE,
               coef = 1.5) + 
  geom_point(aes(fill=factor(id)), size=3, shape=21, color='black') +
  geom_vline(xintercept=2.5, linetype = "dotted") + 
  geom_vline(xintercept=8.5, linetype = "dotted") +
  scale_x_continuous(breaks = c(0.5,1.5,3:11), labels=c("Arm", "Lesion", '0', '10', '30\n    Lesion Area Skin', '40', '180', '190', "10", "40 \n Arm", "190")) +
  scale_color_manual(values=person_color, guide = 'none') +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle(CD8^{"+"}~T~cells) + labs(y=Cells~per~mm^{2}, x="Biopsy time points (day)") + 
  theme_bw() + 
  theme(axis.text=element_text(size=10, hjust = 0.5),
        axis.title=element_text(size=10, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, label.y = c(190, 205, 220, 130, 130, 150, 170),
                             comparisons = my_comparisons_paired, paired = TRUE, method = "wilcox", 
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE) + 
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, label.y = c(230, 130, 185),
                             comparisons = my_comparisons_unpaired, paired = FALSE, method = "wilcox", 
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE)

pdf(figure_1c_filename, width = 6, height = 4)
plot
dev.off()

## Stats for comparisons - 
cd4 <- filter(tcells, cd == "4")
cd8 <- filter(tcells, cd == "8")

wilcox.test(tcells$control0, tcells$control10, paired = T) #p = 0.005
wilcox.test(tcells$control0, tcells$control40, paired = T) #p = 0.003
wilcox.test(tcells$control0, tcells$control190, paired = T) #p = 0.01
wilcox.test(tcells$d0, tcells$d10, paired = T) #p = 0.98
wilcox.test(tcells$d0, tcells$d30, paired = T) #p = 0.23
wilcox.test(tcells$d0, tcells$d40, paired = T) 
wilcox.test(tcells$d0, tcells$d180, paired = T) 
wilcox.test(tcells$d0, tcells$d190, paired = T) # p = 0.02

wilcox.test(cd4$control0, cd4$control10, paired = T) #p = 0.018
wilcox.test(cd4$control0, cd4$control40, paired = T) #p = 0.012
wilcox.test(cd4$control0, cd4$control190, paired = T) #p = 0.1

wilcox.test(cd4$control0, cd4$d0, paired = T) #p = 0.01
wilcox.test(cd4$control0, cd4$d10, paired = T) #p = 0.007
wilcox.test(cd4$control0, cd4$d30, paired = T) #p = 0.004
wilcox.test(cd4$control0, cd4$d40, paired = T) #p = 0.01
wilcox.test(cd4$control0, cd4$d180, paired = T) #p = 0.007
wilcox.test(cd4$control0, cd4$d190, paired = T) #p = 0.03
mean(cd4$control0) #14.3
mean(cd4$control10) #40.0
mean(cd4$control40) #49.8

wilcox.test(cd8$control0, cd8$control10, paired = T) #p = 0.26
wilcox.test(cd8$control0, cd8$control40, paired = T) #p = 0.11
wilcox.test(cd8$control0, cd8$control190, paired = T) #p = 0.06

wilcox.test(cd8$control0, cd8$d0, paired = T) #p = 0.01
wilcox.test(cd8$control0, cd8$d10, paired = T) #p = 0.004
wilcox.test(cd8$control0, cd8$d30, paired = T) #p = 0.009
wilcox.test(cd8$control0, cd8$d40, paired = T) #p = 0.01
wilcox.test(cd8$control0, cd8$d180, paired = T) #p = 0.03
wilcox.test(cd8$control0, cd8$d190, paired = T) #p = 0.03

## Figure 1ghi 
clonal <- read.csv(filename_f1ghi, header = T, stringsAsFactors = FALSE)

## prep 
filename = file.path(repo_loc, '/data/vax_nt_all_dei.csv')
vax <- read_csv(filename)
nunique <- vax %>%
  group_by(person) %>%
  summarise(g0 = sum(g0>0), g10 = sum(g10>0), 
            a0 = sum(a0>0), a10 = sum(a10>0), 
            a30 = sum(a30>0), a40 = sum(a40>0), 
            a80 = sum(a80>0), a90 = sum(a90>0), 
            c0 = sum(c0>0), les = sum(les>0),
            c10 = sum(c10>0), c40 = sum(c40>0), 
            c90 = sum(c90>0))
wilcox.test(nunique$a10, nunique$c10, paired = T, exact = T)
nabove4 <- vax %>%
  group_by(person) %>%
  summarise(g0 = sum(g0>4), g10 = sum(g10>4), 
            a0 = sum(a0>4), a10 = sum(a10>4), 
            a30 = sum(a30>4), a40 = sum(a40>4), 
            a80 = sum(a80>4), a90 = sum(a90>4), 
            c0 = sum(c0>4), les = sum(les>4),
            c10 = sum(c10>4), c40 = sum(c40>4), 
            c90 = sum(c90>4))
wilcox.test(nabove4$a10, nabove4$c10, paired = T, exact = T) #.02
wilcox.test(nabove4$a40, nabove4$c40, paired = T, exact = T) #.12
wilcox.test(nabove4$a90, nabove4$c90, paired = T, exact = T) #.03

#wide versions for stats
nunique_wide <- pivot_wider(clonal[,c(1:3)], names_from = time, values_from = n_unique)
nabove4_wide <- pivot_wider(clonal[,c(1,2,5)], names_from = time, values_from = n_above4)
clonality_wide <- pivot_wider(clonal[,c(1,2,4)], names_from = time, values_from = clonality)

#long version for graphing 
clonal$time <- factor(clonal$time, levels = c("l0", "l10", "l30", "l40", "l80", "l90", "c10", "c40", "c90"),
                         labels = c(1,2,3,4,5,6,7,8,9))
colnames(clonal)

my_comparisons_paired <- list(c(1, 2), c(3, 4), c(7, 8), c(2, 7), c(4, 8), c(6, 9))
my_comparisons_unpaired <- list(c(5, 6), c(8, 9))

# Fig 1g (all unique)
p <- ggplot(data = clonal, aes(x=time, y=n_unique, group = time)) + 
  geom_boxplot(color="black", fill="#bababa", width = 0.75, na.rm=TRUE, coef = 1.5) + 
  geom_point(aes(fill=factor(person)), size=3, shape=21, color='black') +
  geom_vline(xintercept=6.5, linetype = "dotted") + 
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7,8,9), 
                   labels=c('0', '10', '30\n     Lesion Area Skin', '40', '180', '190', '10', '40\nArm', '190')) +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle("All unique TCR clonotypes") + 
  labs(y="Unique productive TCR sequences", x="Biopsy time points (day)") +
  theme_bw() + ylim(0,8400) + 
  theme(axis.text=element_text(size=10, hjust = 0.5),
        axis.title=element_text(size=10, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggpubr::stat_compare_means(mapping = aes(group=t_plot),  vjust = 0, 
                             comparisons = my_comparisons_paired, paired = TRUE, method = "wilcox.test", methods.args = "exact = TRUE",
                             label.y = c(6600, 6600, 6600, 7000, 7400, 7800),
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE) + 
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, 
                           comparisons = my_comparisons_unpaired, paired = FALSE, method = "wilcox.test", methods.args = "exact = TRUE",
                           size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE, label.y = c(6600, 6800))

pdf(figure_1g_filename, width = 6, height = 4)
p 
dev.off()

### stats (confirm exact p value in base R with stat_compare_means)
colnames(nunique_wide)
wilcox.test(nunique_wide$c10, nunique_wide$l10, paired = T) #p = 0.0744
wilcox.test(nunique_wide$c40, nunique_wide$l40, paired = T) #p = 0.35
wilcox.test(nunique_wide$c90, nunique_wide$l90, paired = T) #p = 0.0156
wilcox.test(nunique_wide$l0, nunique_wide$l10, paired = T) #p = 1
wilcox.test(nunique_wide$l30, nunique_wide$l40, paired = T) #p = 0.9
wilcox.test(nunique_wide$l80, nunique_wide$l90, paired = T) #p = 0.8
wilcox.test(nunique_wide$c10, nunique_wide$c40, paired = T) #p = 0.4
wilcox.test(nunique_wide$c40, nunique_wide$c90, paired = T) #p = 0.0156

# Fig 1h (number of unique clonotypes observed at greater than 4 copies )
p <- ggplot(data = clonal, aes(x=time, y=n_above4, group = time)) + 
  geom_boxplot(color="black", fill="#bababa", width = 0.75, na.rm=TRUE, coef = 1.5) + 
  geom_point(aes(fill=factor(person)), size=3, shape=21, color='black') +
  geom_vline(xintercept=6.5, linetype = "dotted") + 
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7,8,9), 
                   labels=c('0', '10', '30\n     Lesion Area Skin', '40', '180', '190', '10', '40\nArm', '190')) +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle("TCR clonotypes detected at >4 copies") + 
  labs(y="Unique productive TCR sequences", x="Biopsy time points (day)") + 
  theme_bw() + 
  theme(axis.text=element_text(size=10, hjust = 0.5),
        axis.title=element_text(size=10, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggpubr::stat_compare_means(mapping = aes(group=t_plot),  vjust = 0, 
                             comparisons = my_comparisons_paired, paired = TRUE, method = "wilcox", label.y = c(690, 690, 690, 730, 780, 820),
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE) + 
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, 
                             comparisons = my_comparisons_unpaired, paired = FALSE, method = "wilcox", 
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE, label.y = c(690, 700))

pdf(figure_1h_filename, width = 6, height = 4)
p 
dev.off()

### stats (confirm exact p value in base R with stat_compare_means)
colnames(nabove4_wide)
wilcox.test(nabove4_wide$c10, nabove4_wide$l10, paired = T) #p = 0.028
wilcox.test(nabove4_wide$c40, nabove4_wide$l40, paired = T) #p = 0.129
wilcox.test(nabove4_wide$c90, nabove4_wide$l90, paired = T) #p = 0.031
wilcox.test(nabove4_wide$l0, nabove4_wide$l10, paired = T) #p = 0.73
wilcox.test(nabove4_wide$l30, nabove4_wide$l40, paired = T) #p = 0.16
wilcox.test(nabove4_wide$l80, nabove4_wide$l90, paired = T) #p = 0.498
wilcox.test(nabove4_wide$c10, nabove4_wide$c40, paired = T) #p = 0.476
wilcox.test(nabove4_wide$c40, nabove4_wide$c90, paired = T) #p = 0.042

# Fig 1i (clonality)
p <- ggplot(data = clonal, aes(x=time, y=clonality, group = time)) + 
  geom_boxplot(color="black", fill="#bababa", width = 0.75, na.rm=TRUE, coef = 1.5) + 
  geom_point(aes(fill=factor(person)), size=3, shape=21, color='black') +
  geom_vline(xintercept=6.5, linetype = "dotted") + 
  scale_x_discrete(breaks = c(1,2,3,4,5,6,7,8,9), 
                   labels=c('0', '10', '30\n     Lesion Area Skin', '40', '180', '190', '10', '40\nArm', '190')) +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle("Repertoire clonality") + 
  labs(y="Clonality (1-normalized Shannon entropy)", x="Biopsy time points (day)") +
  theme_bw() + ylim(0,0.32) + 
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=10,face="bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5)) +
  ggpubr::stat_compare_means(mapping = aes(group=t_plot),  vjust = 0, 
                             comparisons = my_comparisons_paired, paired = TRUE, method = "wilcox", label.y = c(.24, .24, .24, .265, .285, .3),
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE) + 
  ggpubr::stat_compare_means(mapping = aes(group=t_plot), vjust = 0, 
                             comparisons = my_comparisons_unpaired, paired = FALSE, method = "wilcox", 
                             size = 3, label = "p.signif", label.y.npc = 0, tip.length = 0, hide.ns = FALSE, label.y = c(.24, .25))

pdf(figure_1i_filename, width = 6, height = 4)
p 
dev.off()

### stats (confirm exact p value in base R with stat_compare_means)
colnames(clonality_wide)
wilcox.test(clonality_wide$c10, clonality_wide$l10, paired = T) #p = 0.039
wilcox.test(clonality_wide$c40, clonality_wide$l40, paired = T) #p = 0.129
wilcox.test(clonality_wide$c90, clonality_wide$l90, paired = T) #p = 0.156
wilcox.test(clonality_wide$l0, clonality_wide$l10, paired = T)  #p = 1
wilcox.test(clonality_wide$l30, clonality_wide$l40, paired = T) #p = 0.074
wilcox.test(clonality_wide$l80, clonality_wide$l90, paired = T) #p = 0.22
wilcox.test(clonality_wide$c10, clonality_wide$c40, paired = T) #p = 0.359
wilcox.test(clonality_wide$c40, clonality_wide$c90, paired = T) #p = 0.219
