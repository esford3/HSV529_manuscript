
#code/HSV_529_Fig_3.R  

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(data.table)
library(stringr)

rm(list=ls())

# INPUTS
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript/'
} else {
  stop("set repo loc and repo manually")
}

filename_f3 = file.path(repo_loc, '/data/vax_nt_all_dei.csv') ## too big for github, available online

## OUTPUTS
figure_3b_filename = file.path(repo_loc, 'figures/fig_3b1.pdf')
figure_3b_subset_filename = file.path(repo_loc, 'figures/fig_3b2.pdf')
figure_3c_filename = file.path(repo_loc, 'figures/fig_3c1.pdf')
figure_3c_subset_filename = file.path(repo_loc, 'figures/fig_3c2.pdf')
figure_3d_filename = file.path(repo_loc, 'figures/fig_3d.pdf')
figure_3e_filename = file.path(repo_loc, 'figures/fig_3e.pdf')

## PREP
vax <- read_csv(filename_f3)
map <- filter(vax, a0!=0 & a10/a0 >=6)
map$blood <- ifelse(map$g0>0 | map$g10>0, "TRUE", "FALSE")
blood <- c("TRUE" = "red", "FALSE" = "black", 'NA' = "black")
map$blood[(is.na(map$blood))] <- "FALSE"

## Figure 3b
plot <- map %>% 
  select(person, a0, a10, a30, a40, a80, a90, c0, c10, c40, c90, blood) %>%
  mutate(uid = seq_along(.$person)) %>%
  mutate(c0 = NA) %>%
  tidyr::gather(time, count, -person, -uid, -blood) %>%
  mutate(xpos = factor(time, levels = c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90"))) %>%
  mutate(person = factor(person)) %>%
  mutate(value = ifelse(count == 0, 0.1, count)) %>%
  mutate(value = ifelse(value >1000, 1000, value)) %>%
  ggplot(aes(x = xpos, y = value, group = uid, col = blood)) +
  geom_line(size = .3, alpha = 0.5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                limits = c(0.1, 1000),
                labels = c("ND", 1, 10, 100, 1000)) +
  scale_x_discrete(labels=c('0', '10', '30\n        Study day', '40', '180', '190', '',  "10", "40 \n Study day", "190"),
                   breaks=c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90")) +
  theme_classic() +
  facet_wrap(~person, ncol = 5, labeller = label_parsed) +
  annotation_logticks(side = "l" , size = .1) +
  scale_color_manual(values = blood) +
  theme(legend.position = "none") +
  geom_text(aes(label = paste("participant", person), x = 8, y = 500), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  xlab("") + 
  ylab("Prevalent clonotypes exp >=6 over dose 1 \n(clonotype copies)")
pdf(figure_3b_filename, width = 10, height = 2.5)
plot 
dev.off()

##need dummy set as subset only has person 2 and 4 
df1<-data.frame(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df2<-data.frame(2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df3<-data.frame(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df4<-data.frame(4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df5<-data.frame(5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df6<-data.frame(6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df7<-data.frame(7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df8<-data.frame(8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
df9<-data.frame(9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, "FALSE")
list <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9)
df <- rbindlist(list, use.names = FALSE)
colnames(df)<- c("person", "a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90", "blood")

plot <- map %>% 
  subset(g0!=0 | g10!=0) %>%
  select(person, a0, a10, a30, a40, a80, a90, c0, c10, c40, c90, blood) %>%
  mutate(uid = seq_along(.$person)) %>%
  mutate(c0 = NA) %>%
  bind_rows(df) %>%
  tidyr::gather(time, count, -person, -uid, -blood) %>%
  mutate(xpos = factor(time, levels = c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90"))) %>%
  mutate(person = factor(person)) %>%
  mutate(value = ifelse(count == 0, 0.1, count)) %>%
  mutate(value = ifelse(value >1000, 1000, value)) %>%
  ggplot(aes(x = xpos, y = value, group = uid, col = blood)) +
  geom_line(size = .3, alpha = 0.5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                limits = c(0.1, 1000),
                labels = c("ND", 1, 10, 100, 1000)) +
  scale_x_discrete(labels=c('0', '10', '30\n        Study day', '40', '180', '190', '',  "10", "40 \n Study day", "190"),
                   breaks=c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90")) +
  theme_classic() +
  facet_wrap(~person, ncol = 5, labeller = label_parsed) +
  annotation_logticks(side = "l" , size = .1) +
  scale_color_manual(values = blood) +
  theme(legend.position = "none") +
  geom_text(aes(label = paste("participant", person), x = 8, y = 500), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  xlab("") + 
  ylab("Prevalent clonotypes exp >=6 over dose 1 \n(clonotype copies)")
pdf(figure_3b_subset_filename, width = 10, height = 2.5)
plot 
dev.off()

## Figure 3c  
map <- filter(vax, a0==0 & a10 >=6)
map$blood <- ifelse(map$g0>0 | map$g10>0, "TRUE", "FALSE")
blood <- c("TRUE" = "red", "FALSE" = "black")
map$blood[(is.na(map$blood))] <- "FALSE"
map <- select(map, person, a0, a10, a30, a40, a80, a90, c0, c10, c40, c90, blood)

## elicits are missing person 1 and person 5 - still need dummy rows 
plot <- map %>% 
  select(person, a0, a10, a30, a40, a80, a90, c0, c10, c40, c90, blood) %>%
  mutate(uid = seq_along(.$person)) %>%
  mutate(c0 = NA) %>%
  bind_rows(df) %>% 
  tidyr::gather(time, count, -person, -uid, -blood) %>%
  mutate(xpos = factor(time, levels = c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90"))) %>%
  mutate(person = factor(person)) %>%
  mutate(value = ifelse(count == 0, 0.1, count)) %>%
  mutate(value = ifelse(value >1000, 1000, value)) %>%
  ggplot(aes(x = xpos, y = value, group = uid, col = blood)) +
  geom_line(size = .3, alpha = 0.5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                limits = c(0.1, 1000),
                labels = c("ND", 1, 10, 100, 1000)) +
  scale_x_discrete(labels=c('0', '10', '30\n        Study day', '40', '180', '190', '',  "10", "40 \n Study day", "190"),
                   breaks=c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90")) +
  theme_classic() +
  facet_wrap(~person, ncol = 5, labeller = label_parsed) +
  annotation_logticks(side = "l" , size = .1) +
  scale_color_manual(values = blood) +
  theme(legend.position = "none") +
  geom_text(aes(label = paste("participant", person), x = 8, y = 500), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  xlab("") + 
  ylab("Elicited clonotypes >=6 over dose 1 \n(clonotype copies)")
pdf(figure_3c_filename, width = 10, height = 2.5)
plot 
dev.off()

plot <- map %>% 
  subset(blood == TRUE) %>%
  select(person, a0, a10, a30, a40, a80, a90, c0, c10, c40, c90, blood) %>%
  mutate(uid = seq_along(.$person)) %>%
  mutate(c0 = NA) %>%
  bind_rows(df) %>%
  tidyr::gather(time, count, -person, -uid, -blood) %>%
  mutate(xpos = factor(time, levels = c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90"))) %>%
  mutate(person = factor(person)) %>%
  mutate(value = ifelse(count == 0, 0.1, count)) %>%
  mutate(value = ifelse(value >1000, 1000, value)) %>%
  ggplot(aes(x = xpos, y = value, group = uid, col = blood)) +
  geom_line(size = .3) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                limits = c(0.1, 1000),
                labels = c("ND", 1, 10, 100, 1000)) +
  scale_x_discrete(labels=c('0', '10', '30\n        Study day', '40', '180', '190', '',  "10", "40 \n Study day", "190"),
                   breaks=c("a0", "a10", "a30", "a40", "a80", "a90", "c0", "c10", "c40", "c90")) +
  theme_classic() +
  facet_wrap(~person, ncol = 5, labeller = label_parsed) +
  annotation_logticks(side = "l" , size = .1) +
  scale_color_manual(values = blood) +
  theme(legend.position = "none") +
  geom_text(aes(label = paste("participant", person), x = 8, y = 500), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.title = element_text(size = 8)) +
  theme(axis.text = element_text(size = 8)) +
  xlab("") + 
  ylab("Elicited clonotypes >=6 over dose 1 \n(clonotype copies)")
plot
pdf(figure_3c_subset_filename, width = 10, height = 2.5)
plot 
dev.off()


## Figure 3d
person_color <-  c("1" = "#8d03a9", "2" = "#FFAD4A", "3" = "#0263be", 
                   "4" = "#FFE900", "5" = "#f93f61", "6" = "#6a4ee3", 
                   "7" = "#2e9a04", "8" = "#c8c89e", "9" = "#94ee78")
str(vax)
vax$prev <- ifelse(vax$a0>0, "TRUE", "FALSE")
vax$dose1 <- ifelse(vax$a0!=0 & vax$a10!=0, vax$a10/vax$a0, 
                      ifelse(vax$a0==0 & vax$a10!=0, vax$a10,
                             ifelse(vax$a10==0, 0, NA)))
table(vax$dose1)
vax$dose2 <- ifelse(vax$a30!=0 & vax$a40!=0, vax$a40/vax$a30, 
                      ifelse(vax$a30==0 & vax$a40!=0, vax$a40,
                             ifelse(vax$a40==0, 0, NA)))
table(vax$dose2)

vax$dose3 <- ifelse(vax$a80>0 & vax$a90>0, vax$a90/vax$a80, 
                      ifelse((vax$a80==0 | is.na(vax$a80)) & vax$a90>0, vax$a90,
                             ifelse(vax$a90==0, 0, NA)))
table(vax$dose3)

matrix <- vax %>% 
  group_by(person) %>%
  summarise(dose1 = NROW(dose1[dose1>=6]),
            dose2 = NROW(dose2[dose2>=6]),
            dose3 = NROW(dose3[dose3>=6 & !is.na(dose3)]))

m_long <- matrix %>% select(person, dose1, dose2, dose3) %>%
  pivot_longer(cols = c(2:4)) %>%
  mutate(xpos = name) %>%
  mutate(xpos = factor(xpos, levels = c("dose1", "dose2", "dose3")))

my_comparisons <- list(c("dose1", "dose2"), c("dose2", "dose3"), c("dose1", "dose3"))

a <- ggplot(m_long, aes(x = xpos, y = value, fill = factor(person), group = person)) +
  geom_bar(stat = "identity", show.legend = T, na.rm = T) +
  theme_bw() + ggtitle("Expanding clonotypes by dose") +
  scale_x_discrete(labels = c("dose1" = "Dose 1",
                              "dose2" = "Dose 2",
                              "dose3" = "Dose 3")) + 
  labs(y="Unique clonotypes expanding >6", x= "") +
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  theme_classic() + 
  scale_fill_manual(values = person_color, name = "Person", 
                    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9")) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, paired = TRUE, na.rm = FALSE, 
                             method = "wilcox", size = 4, label = "p.signif", label.y = c(550, 500, 590)) +
  ggpubr::stat_compare_means(method = "wilcox", size = 4, label = "p.signif", 
                             col = "black", hide.ns = TRUE)

pdf(figure_3d_filename, width = 3.5, height = 5)
a 
dev.off()

## Figure 3e
matrix <- vax %>% 
  group_by(person, prev) %>%
  summarise(dose1 = NROW(dose1[dose1>=6]),
            dose2 = NROW(dose2[dose2>=6]),
            dose3 = NROW(dose3[dose3>=6 & !is.na(dose3)]))

m_long <- matrix %>%
  pivot_longer(cols = c(3:5)) %>%
  mutate(prev = str_replace_all(prev, c("TRUE" = "prev", "FALSE" = "elic"))) %>%
  mutate(xpos = paste(prev, name, sep = "_")) %>%
  mutate(xpos = factor(xpos, levels = c("prev_dose1", "elic_dose1", 
                                        "prev_dose2", "elic_dose2", 
                                        "prev_dose3", "elic_dose3")))
my_comparisons <- list(c("prev_dose1", "elic_dose1"), c("prev_dose2", "elic_dose2"), c("prev_dose3", "elic_dose3"))

a <- ggplot(m_long, aes(x = xpos, y = value, fill = factor(person), group = person)) +
  geom_bar(stat = "identity", show.legend = T, na.rm = T) +
  theme_bw() + ggtitle("Expanding clonotypes by dose and detection at day 0") +
  labs(x="Dose 1                Dose 2                Dose 3", 
       y="Unique clonotypes expanding >6") +
  theme(plot.title = element_text(size = 12, hjust = 0.5)) +
  scale_x_discrete(labels = c("Prevalent", "Elicited","Prevalent", "Elicited", "Prevalent", "Elicited"), 
                   guide = guide_axis(angle = 45)) + 
  theme_classic() + 
  scale_fill_manual(values = person_color, name = "Person", 
                    labels = c("1", "2", "3", "4", "5", "6", "7", "8", "9")) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons, paired = FALSE, na.rm = FALSE, 
                             method = "wilcox", size = 4, label = "p.signif", label.y = c(275, 275, 275))+
  ggpubr::stat_compare_means(method = "wilcox", size = 4, label = "p.signif", hide.ns = TRUE, col = "black")
  
pdf(figure_3e_filename, width = 5, height = 5)
a 
dev.off()
