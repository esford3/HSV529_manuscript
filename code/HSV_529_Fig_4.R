

#code/HSV_529_Fig_4.R  

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(readxl)
library(ggrepel)

rm(list=ls())

# INPUTS
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript/'
} else {
  stop("set repo loc and repo manually")
}

filename_f4 = file.path(repo_loc, '/data/vax_nt_all_dei.csv') ## too big for github, available online
color_path = file.path(repo_loc, 'data/color_code.xlsx')

## OUTPUTS
figure_4a_filename = file.path(repo_loc, 'figures/fig_4a.pdf')
figure_s3b_filename = file.path(repo_loc, 'figures/fig_s3b.pdf')

## pivot to long first by time, then filter 
vax <- read_csv(filename_f4)
colnames(vax)
map <- vax[,c("v", "j", "aminoAcid", "g0", "g10", "a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10",
                     "c40", "c90", "les", "person")]
map <- mutate(map, person = factor(person))

## condense v down to just family 
map$v <- sapply(strsplit(map$v,"-"), `[`, 1)
map$v <- gsub("TCRBV", "TRBV", map$v)
map$j <- gsub("TCRBJ", "TRBJ", map$j)

## make two summary charts for v family and j gene usage by time
map <- filter(map, !is.na(v))
map_v <- map %>%
  group_by(person, v) %>%
  summarise(g0_unique = NROW(g0[g0>0 & !is.na(g0)]),
            g0_sum = sum(g0[g0>0 & !is.na(g0)]),
            g10_unique = NROW(g10[g10>0 & !is.na(g10)]),
            g10_sum = sum(g10[g10>0 & !is.na(g10)]),
            a0_unique = NROW(a0[a0>4 & !is.na(a0)]),
            a0_sum = sum(a0[a0>4 & !is.na(a0)]),
            a10_unique = NROW(a10[a10>4 & !is.na(a10)]),
            a10_sum = sum(a10[a10>4 & !is.na(a10)]),
            a30_unique = NROW(a30[a30>4 & !is.na(a30)]),
            a30_sum = sum(a30[a30>4 & !is.na(a30)]),
            a40_unique = NROW(a40[a40>4 & !is.na(a40)]),
            a40_sum = sum(a40[a40>4 & !is.na(a40)]),
            a80_unique = NROW(a80[a80>4 & !is.na(a80)]),
            a80_sum = sum(a80[a80>4 & !is.na(a80)]),
            a90_unique = NROW(a90[a90>4 & !is.na(a90)]),
            a90_sum = sum(a90[a90>4 & !is.na(a90)]),
            c0_unique = NROW(c0[c0>4 & !is.na(c0)]),
            c0_sum = sum(c0[c0>4 & !is.na(c0)]),
            c10_unique = NROW(c10[c10>4 & !is.na(c10)]),
            c10_sum = sum(c10[c10>4 & !is.na(c10)]),
            c40_unique = NROW(c40[c40>4 & !is.na(c40)]),
            c40_sum = sum(c40[c40>4 & !is.na(c40)]),
            c90_unique = NROW(c90[c90>4 & !is.na(c90)]),
            c90_sum = sum(c90[c90>4 & !is.na(c90)]),
            les_unique = NROW(les[les>4 & !is.na(les)]),
            les_sum = sum(les[les>4 & !is.na(les)]))
map_v <- map_v %>%
  group_by(person) %>%
  mutate(g0_percent = g0_sum/sum(g0_sum),
         g10_percent = g10_sum/sum(g10_sum),
         a0_percent = a0_sum/sum(a0_sum),
         a10_percent = a10_sum/sum(a10_sum),
         a30_percent = a30_sum/sum(a30_sum),
         a40_percent = a40_sum/sum(a40_sum),
         a80_percent = a80_sum/sum(a80_sum),
         a90_percent = a90_sum/sum(a90_sum),
         c0_percent = c0_sum/sum(c0_sum),
         c10_percent = c10_sum/sum(c10_sum),
         c40_percent = c40_sum/sum(c40_sum),
         c90_percent = c90_sum/sum(c90_sum),
         les_percent = les_sum/sum(les_sum))

map <- filter(map, !is.na(j))
map_j <- map %>%
  group_by(person, j) %>%
  summarise(g0_unique = NROW(g0[g0>0 & !is.na(g0)]),
            g0_sum = sum(g0[g0>0 & !is.na(g0)]),
            g10_unique = NROW(g10[g10>0 & !is.na(g10)]),
            g10_sum = sum(g10[g10>0 & !is.na(g10)]),
            a0_unique = NROW(a0[a0>4 & !is.na(a0)]),
            a0_sum = sum(a0[a0>4 & !is.na(a0)]),
            a10_unique = NROW(a10[a10>4 & !is.na(a10)]),
            a10_sum = sum(a10[a10>4 & !is.na(a10)]),
            a30_unique = NROW(a30[a30>4 & !is.na(a30)]),
            a30_sum = sum(a30[a30>4 & !is.na(a30)]),
            a40_unique = NROW(a40[a40>4 & !is.na(a40)]),
            a40_sum = sum(a40[a40>4 & !is.na(a40)]),
            a80_unique = NROW(a80[a80>4 & !is.na(a80)]),
            a80_sum = sum(a80[a80>4 & !is.na(a80)]),
            a90_unique = NROW(a90[a90>4 & !is.na(a90)]),
            a90_sum = sum(a90[a90>4 & !is.na(a90)]),
            c0_unique = NROW(c0[c0>4 & !is.na(c0)]),
            c0_sum = sum(c0[c0>4 & !is.na(c0)]),
            c10_unique = NROW(c10[c10>4 & !is.na(c10)]),
            c10_sum = sum(c10[c10>4 & !is.na(c10)]),
            c40_unique = NROW(c40[c40>4 & !is.na(c40)]),
            c40_sum = sum(c40[c40>4 & !is.na(c40)]),
            c90_unique = NROW(c90[c90>4 & !is.na(c90)]),
            c90_sum = sum(c90[c90>4 & !is.na(c90)]),
            les_unique = NROW(les[les>4 & !is.na(les)]),
            les_sum = sum(les[les>4 & !is.na(les)]))
map_j <- map_j %>%
  group_by(person) %>%
  mutate(g0_percent = g0_sum/sum(g0_sum),
         g10_percent = g10_sum/sum(g10_sum),
         a0_percent = a0_sum/sum(a0_sum),
         a10_percent = a10_sum/sum(a10_sum),
         a30_percent = a30_sum/sum(a30_sum),
         a40_percent = a40_sum/sum(a40_sum),
         a80_percent = a80_sum/sum(a80_sum),
         a90_percent = a90_sum/sum(a90_sum),
         c0_percent = c0_sum/sum(c0_sum),
         c10_percent = c10_sum/sum(c10_sum),
         c40_percent = c40_sum/sum(c40_sum),
         c90_percent = c90_sum/sum(c90_sum),
         les_percent = les_sum/sum(les_sum))

## get color scale
color <- read_xlsx(color_path)
color$color <- gsub("beige", "navy", color$color)
color$color <- gsub("olive", "olivedrab", color$color)
color$color <- gsub("indigo", "mediumpurple4", color$color)
color$color <- gsub("lime", "olivedrab1", color$color)
color$color <- gsub("teal", "darkcyan", color$color)
color <- unique(color$color)

require(data.table)
map_v_long <- pivot_longer(map_v, cols  = -c(person, v), names_to = c("time", "type"), 
                           names_sep = "\\_")
map_v_wide <- pivot_wider(map_v_long, id_cols = c(person, v, time), names_from = type)
map_v_wide$time = factor(map_v_wide$time, levels = c("g0", "g10", "c0", "c10", "c40", "c90", "les", "a0", "a10", "a30", "a40", "a80", "a90"))

map_j_long <- pivot_longer(map_j, cols  = -c(person, j), names_to = c("time", "type"), 
                           names_sep = "\\_")
map_j_wide <- pivot_wider(map_j_long, id_cols = c(person, j, time), names_from = type)
map_j_wide$time = factor(map_j_wide$time, levels = c("g0", "g10", "c0", "c10", "c40", "c90", "les", "a0", "a10", "a30", "a40", "a80", "a90"))

## check percentages
sum(map_v$a10_percent[map_v$person == "1"]) #(these all sum to 1)

## start with v families
colnames(map_v_wide)
pie = ggplot(data = map_v_wide, aes(x="", y=percent, fill=v)) + geom_bar(stat="sum", width=1) + 
  coord_polar("y", start=90) + 
  #geom_text_repel(aes(label = ifelse(percent>0.1, paste0(round(percent*100), "%"), "")), position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values = color) + 
  guides(fill = guide_legend(ncol=1)) +
  labs(x = NULL, y = NULL, fill = "V gene usage", title = "V gene usage by person \n(>4 copies in tissue)") + 
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666")) +
  theme(legend.key.size = unit(0.2, 'in')) + 
  facet_grid(cols = vars(time), rows = vars(person), labeller = 
               labeller(.multi_line = FALSE), 
  )
pdf(figure_4a_filename, height = 8, width = 12)
pie
dev.off()

pie = ggplot(data = map_j_wide, aes(x="", y=percent, fill=j)) + geom_bar(stat="sum", width=1) + 
  coord_polar("y", start=90) + 
  #geom_text_repel(aes(label = ifelse(percent>0.1, paste0(round(percent*100), "%"), "")), position = position_stack(vjust = 0.5)) + 
  scale_fill_manual(values = color) + 
  guides(fill = guide_legend(ncol=1)) +
  labs(x = NULL, y = NULL, fill = "J gene usage", title = "J gene usage by person \n(>4 copies in tissue)") + 
  theme_classic() + 
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5, color = "#666666")) +
  theme(legend.key.size = unit(0.2, 'in')) + 
  facet_grid(cols = vars(time), rows = vars(person), labeller = 
               labeller(.multi_line = FALSE), 
  )
pdf(figure_s3b_filename, height = 8, width = 12)
pie
dev.off()



