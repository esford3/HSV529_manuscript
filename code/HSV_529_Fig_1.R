
# code/HSV_529_Fig_1.R

library(ggplot2)
library(dplyr)
library(tidyr)

rm(list=ls())

# INPUTS
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript/'
} else {
  stop("set repo loc and repo manually")
}

filename_f1bc = file.path(repo_loc, '/data/tcellcounts.csv')
filename_f1ghi = file.path(repo_loc, '/data/fig1_ghi.csv')

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
graphable$t_plot = rep(c(0.5,1.5,3:11), 18)
graph.cd4 = filter(graphable, cd == 4)
graph.cd8 = filter(graphable, cd == 8)

## CD4
p <- ggplot(data = graph.cd4, aes(x=t_plot, y=cells, group = t_plot), theme_bw()) + 
  geom_boxplot(color="black", fill="#bababa", outlier.shape=NA, width = 0.75, na.rm=TRUE,
               coef = 1.5) + 
  geom_point(aes(fill=factor(id)), size=3, shape=21, color='black') +
  geom_vline(xintercept=2.5, linetype = "dotted") + 
  geom_vline(xintercept=8.5, linetype = "dotted") +
  scale_x_continuous(breaks = c(0.5,1.5,3:11), labels=c("Arm   ", "  Lesion", '0', '10', '30\n    Lesion Area Skin', '40', '180', '190', "10", "40 \n Arm", "190")) +
  scale_color_manual(values=person_color, guide = 'none') +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle(CD4^{"+"}~T~cells) + labs(y=Cells~per~mm^{2}, x="Biopsy time points (day)") + 
  theme_bw() + ylim(0,350) +
  theme(axis.text=element_text(size=12, hjust = 0.5),
        axis.title=element_text(size=12, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

pdf(figure_1b_filename, width = 6, height = 4)
p 
dev.off()

## CD8 
p <- ggplot(data = graph.cd8, aes(x=t_plot, y=cells, group = t_plot), theme_bw()) + 
  geom_boxplot(color="black", fill="#bababa", outlier.shape=NA, width = 0.75, na.rm=TRUE,
               coef = 1.5) + 
  geom_point(aes(fill=factor(id)), size=3, shape=21, color='black') +
  geom_vline(xintercept=2.5, linetype = "dotted") + 
  geom_vline(xintercept=8.5, linetype = "dotted") +
  scale_x_continuous(breaks = c(0.5,1.5,3:11), labels=c("Arm   ", "  Lesion", '0', '10', '30\n    Lesion Area Skin', '40', '180', '190', "10", "40 \n Arm", "190")) +
  scale_color_manual(values=person_color, guide = 'none') +
  scale_fill_manual(values=person_color, guide = guide_legend(title = 'Participant')) +
  ggtitle(CD8^{"+"}~T~cells) + labs(y=Cells~per~mm^{2}, x="Biopsy time points (day)") + 
  theme_bw() + ylim(0,200) +
  theme(axis.text=element_text(size=12, hjust = 0.5),
        axis.title=element_text(size=12, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

pdf(figure_1c_filename, width = 6, height = 4)
p 
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
clonal$time <- factor(clonal$time, levels = c("l0", "l10", "l30", "l40", "l80", "l90", "c10", "c40", "c90"),
                         labels = c(1,2,3,4,5,6,7,8,9))
colnames(clonal)

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
  theme_bw() + 
  theme(axis.text=element_text(size=12, hjust = 0.5),
        axis.title=element_text(size=12, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

pdf(figure_1g_filename, width = 6, height = 4)
p 
dev.off()

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
  theme(axis.text=element_text(size=12, hjust = 0.5),
        axis.title=element_text(size=12, face="bold", hjust = 0.5),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

pdf(figure_1h_filename, width = 6, height = 4)
p 
dev.off()


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
  theme_bw() + ylim(0,0.25) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

pdf(figure_1i_filename, width = 6, height = 4)
p 
dev.off()

### stats 
wilcox.test(clonality ~ time, clonal, subset = time == "a0" | time == "a90")
wilcox.test(clonality ~ time, clonal, subset = time == "a0" | time == "a10", paired = T)
wilcox.test(clonality ~ time, clonal, subset = time == "a10" | time == "c10", paired = T)
wilcox.test(clonality ~ time, clonal, subset = time == "a40" | time == "c40", paired = T)
wilcox.test(clonality ~ time, clonal, subset = time == "a90" | time == "c90", paired = T)
wilcox.test(clonality ~ time, clonal, subset = time == "a0" | time == "c10", paired = T)