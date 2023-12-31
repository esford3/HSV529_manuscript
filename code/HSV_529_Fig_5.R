
#code/HSV_529_Fig_5.R  
## HLA heatmap, logoplot 

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(viridis)
library(cowplot)
library(ggpubr)

rm(list=ls())

## INPUTS 
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript'
} else {
  stop("set repo loc and repo manually")
}

filename_f5c= file.path(repo_loc, 'data/clustered_trb_blood.csv')
filename_f5d = file.path(repo_loc, 'data/hla_long.csv')

## OUTPUTS
figure_5c_filename = file.path(repo_loc, 'figures/fig_5c.pdf')
figure_5d_filename = file.path(repo_loc, 'figures/fig_5d.pdf')

## Figure 5c (logoplots)

## network diagrams 
require(ggplot2)
require(ggseqlogo)

chemistry2 <- make_col_scheme(
  chars = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'V', "Z"),
  groups = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 7), rep('Insert', 1)),
  cols = c(rep('#109648', 5), rep('#5E239D', 2), rep('#255C99', 3), rep('#D62839', 2), rep('#221E22', 7), rep('#333333', 1)),
  name = "X")

hsv <- read.csv(filename_f5c)
table(hsv$cluster_id)
hsv_0 <- filter(hsv, cluster_id ==0)
hsv_2 <- filter(hsv, cluster_id ==2)
hsv_4 <- filter(hsv, cluster_id ==15)

a <- ggseqlogo(hsv_0$cdr3_b_aa, method = "probability", col_scheme = "chemistry2")
b <- ggseqlogo(hsv_2$cdr3_b_aa, method = "probability", col_scheme = "chemistry2") + theme(legend.position = 'none', axis.text.x = element_blank())
d <- ggseqlogo(hsv_4$cdr3_b_aa, method = "probability", col_scheme = "chemistry2") + theme(legend.position = 'none', axis.text.x = element_blank()) 

grobs <- ggplotGrob(a)$grobs
legend <- get_legend(a + theme(legend.box.margin = margin(0, 0, 0, 0)))
pgrid <- plot_grid(a + theme(legend.position = 'none', axis.text.x = element_blank()),
                   b,d, ncol = 1)
# add legend
p <- plot_grid(pgrid, legend, nrow = 2, rel_heights = c(10, 1), rel_widths = c(3, 1))

pdf(figure_5c_filename, width = 5, height = 4)
p 
dev.off()

## Figure 5d (HLA heatmap)
hla <- read_csv(filename_f5d)
colnames(hla) 

hla <- hla[with(hla, order(ptid)),]

sharedObs <- function(i) {
  p <- do.call(paste, subset(dfs[[i[1]]], select = DRB1))
  q <- do.call(paste, subset(dfs[[i[2]]], select = DRB1))
  t <- do.call(paste, subset(dfs[[i[1]]], select = DQB1))
  u <- do.call(paste, subset(dfs[[i[2]]], select = DQB1))
  v <- do.call(paste, subset(dfs[[i[1]]], select = DPB1))
  w <- do.call(paste, subset(dfs[[i[2]]], select = DPB1))
  length(intersect(p,q))  + length(intersect(t,u)) + length(intersect(v,w))
}

dfs <- split(hla, hla$ptid)
n <- length(dfs)
mat <- `dimnames<-`(matrix(0,n,n),list(names(dfs),names(dfs))) #makes an empty matrix w/ rownames
mat[lower.tri(mat, diag = F)] <- combn(n,2,sharedObs) #calculates the sharing in the lower half of the matrix
res <- t(mat) + mat # replicates in mirror 
res

pdf(figure_5d_filename, width = 4, height = 4)
heatmap(res, Rowv = NA, Colv = NA, scale="none", ##the Colv NA makes it compatible with the other fig
        labRow = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
        labCol = c("1", "2", "3", "4", "5", "6", "7", "8", "9"),
        las = 2, 
        col = viridis(n))
dev.off()
range(res)
#legend(x= "bottomright", legend=c("0 (min)", "1 (mean)", "3 (max)"), title = "Shared HLA-II \nAlleles", 
#       fill=viridis(3), xjust = 1)

