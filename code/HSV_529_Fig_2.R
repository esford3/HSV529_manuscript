
#code/HSV_529_Fig_2.R  

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

rm(list=ls())

# INPUTS
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript/'
} else {
  stop("set repo loc and repo manually")
}

filename_f2 = file.path(repo_loc, '/data/vax_nt_all_dei.csv') ## too big for github, available online

## OUTPUTS
figure_2a_filename = file.path(repo_loc, 'figures/fig_2a.pdf')
figure_2b_filename = file.path(repo_loc, 'figures/fig_2b.pdf')
figure_2c_filename = file.path(repo_loc, 'figures/fig_2c.pdf')
figure_2d_filename = file.path(repo_loc, 'figures/fig_2d.pdf')

## Figure 2a - day 0 histogram 
prev_color <- c("TRUE" = "purple4", "FALSE" = "red")
hist <- read.csv(filename_f2)
hist$a80[is.na(hist$a80)] <- 0
hist$a90[is.na(hist$a90)] <- 0
hist$g0[is.na(hist$g0)] <- 0
hist$blood = ifelse((hist$g0!=0 & hist$a0==0 & hist$a10==0 & hist$a30==0 & hist$a40==0 & hist$a80==0 & hist$a90==0), TRUE, FALSE)
hist$both = ifelse((hist$g0!=0 & hist$blood == FALSE), TRUE, FALSE)
hist$skin = ifelse((hist$a0!=0 & hist$blood == FALSE & hist$both == FALSE), TRUE, FALSE)
hist$prev = ifelse(hist$a0>0, TRUE, FALSE)
table(hist$blood, hist$both)
table(hist$blood, hist$skin)

hist_g <- hist %>%
  select(person, g0, a0, a10, a30, a40, a80, a90, blood, both, skin, prev) %>%
  group_by(person, prev) %>% 
  summarise(blood = NROW(prev[blood == TRUE]),
            both = NROW(prev[both == TRUE]),
            skin = NROW(prev[skin == TRUE])) %>%
  pivot_longer(cols = c(blood, both, skin), names_to = "site") %>% 
  mutate(person = factor(person))

site_shape <- c("blood" = 24, "skin" = 25, "both" = 23)
hist_g0 <- ggplot(hist_g, aes(fill = prev, y = value, x = site)) +
  geom_point(aes(shape = site, fill = prev), size = 2, show.legend = T) + 
  scale_fill_manual(values = prev_color) + 
  scale_shape_manual(values = site_shape) + 
  facet_wrap(~person, nrow = 1) + 
  scale_y_log10(breaks = c(0, 1, 10, 100, 1000, 10000),
                limits = c(0.9, 10000), 
                labels = c("ND", "1", "10", "100", "1000", "10000")) + 
  scale_x_discrete(labels = c("Blood \nonly", "Both", "Skin \nonly")) + 
  annotation_logticks(sides = "l") + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  geom_text(aes(label = paste("Participant", person), x = 2, y = 10000), color = "black", 
            size = 3, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(axis.title = element_text(size = 6)) +
  theme(axis.text = element_text(size = 6)) +
  xlab("") +
  ylab("Number of unique blood-detected clones at day 0 \nby detection in skin")
hist_g0

pdf(figure_2a_filename, width = 11, height = 2)
hist_g0 
dev.off()

## FIGURE 2b,d
#select overlapping clonotypes 
vax <- read.csv(filename_f2)
map <- filter(vax, g0 !=0 | g10 != 0)
map$g0[is.na(map$g0)] <- 0
map$a80[is.na(map$a80)] <- 0
map$a90[is.na(map$a90)] <- 0
map$prev <- ifelse(map$a0>0, "TRUE", "FALSE")
table(is.na(map$a0))
table(map$prev, map$person)
exp <- filter(map, (prev == TRUE & a10/a0 >=6) | (prev == FALSE & a10 >=6))
label_color <-  c("1" = "#8d03a9", "2" = "#FFAD4A", "3" = "#0263be", 
                  "4" = "#FFE900", "5" = "#f93f61", "6" = "#6a4ee3", 
                  "7" = "#2e9a04", "8" = "#c8c89e", "9" = "#94ee78")

## FIGURE 2B
prev_color <- c("FALSE" = "red", "TRUE" = "purple4")
plot <- map %>% 
  filter(g0 > 0) %>%
  select(person, a0, a10, a30, a40, a80, a90, prev) %>%
  mutate(uid = seq_along(.$person)) %>%
  tidyr::gather(time, count, -person, -uid, -prev) %>%
  mutate(xpos = factor(time, levels = c("a0", "a10", "a30", "a40", "a80", "a90"))) %>%
  mutate(person = factor(person)) %>%
  mutate(value = ifelse(count == 0, 0.1, count)) %>%
  ggplot(aes(x = xpos, y = value, group = uid, col = prev)) +
  geom_line(size = .3, alpha = .5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                limits = c(0.1, 1000),
                labels = c("ND", 1, 10, 100, 1000)) +
  scale_x_discrete(labels=c('0', '10', '30\n        Study day', '40', '180', '190'),
                   breaks=c("a0", "a10", "a30", "a40", "a80", "a90")) +
  theme_classic() +
  facet_wrap(~factor(prev, levels = c("TRUE", "FALSE"))+person, nrow = 2, labeller = label_parsed, drop = FALSE) + #if want wide version
  annotation_logticks(side = "l" , size = .1) +
  scale_color_manual(values = prev_color) +
  theme(legend.position = "none") +
  geom_text(aes(label = paste("Participant", person), x = 3, y = 500), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.title = element_text(size = 6)) +
  theme(axis.text = element_text(size = 6)) +
  xlab("") + 
  ylab("Blood-detected clones at day 0 \n(clonotype copies)")
pdf(figure_2b_filename, width = 8.5, height = 3)
plot 
dev.off()

## Figure 2C
prev_color <- c("pre-vax" = "purple4", "post-vax" = "red")
hist <- vax
hist$a80[is.na(hist$a80)] = 0
hist$a90[is.na(hist$a90)] = 0
hist$g10[is.na(hist$g10)] = 0
hist$g0[is.na(hist$g0)] = 0

hist_g <- hist %>%
  select(person, g10, a0, a10, a30, a40, a80, a90) %>%
  mutate(prev = ifelse(a0>0, "pre-vax", "post-vax")) %>%
  group_by(person, prev) %>% 
  summarise(blood = NROW(g10[g10>0 & a10==0]),
            both = NROW(g10[g10>0 & a10>0]),
            skin = NROW(a10[a10>0])) %>%
  pivot_longer(cols = c(blood, both, skin), names_to = "site", values_to = "value") %>% 
  mutate(person = factor(person))

site_shape <- c("blood" = 24, "skin" = 25, "both" = 23)
hist_g10 <- ggplot(hist_g, aes(fill = prev, y = value, x = site)) +
  geom_point(aes(shape = site, fill = prev), size = 2) + 
  scale_fill_manual(values = prev_color) + 
  scale_shape_manual(values = site_shape) + 
  facet_wrap(~person, nrow = 1) + 
  scale_y_log10(breaks = c(0, 1, 10, 100, 1000, 10000),
                labels = c("ND", "1", "10", "100", "1000", "10000")) + 
  scale_x_discrete(labels = c("Blood \nonly", "Both", "Skin \nonly")) + 
  annotation_logticks(sides = "l") + 
  theme_classic() + 
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  geom_text(aes(label = paste("Participant", person), x = 1.75, y = 1E4), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(axis.title = element_text(size = 6)) +
  theme(axis.text = element_text(size = 6)) +
  ylab("Number of unique clones at day 10 \nby site and detection pre-vaccine")
pdf(figure_2c_filename, width = 11, height = 2)
hist_g10 
dev.off()


## Figure 2d
prev_color <- c("TRUE" = "purple4", "FALSE" = "red")
plot <- map %>% filter(g10 > 0) %>%
  select(person, a0, a10, a30, a40, a80, a90, prev) %>%
  mutate(uid = seq_along(.$person)) %>%
  tidyr::gather(time, count, -person, -uid, -prev) %>%
  mutate(xpos = factor(time, levels = c("a0", "a10", "a30", "a40", "a80", "a90"))) %>%
  mutate(person = factor(person)) %>%
  mutate(value = ifelse(count == 0, 0.1, count)) %>%
  ggplot(aes(x = xpos, y = value, group = uid, col = prev)) +
  geom_line(size = .3, alpha = .5) +
  scale_y_log10(breaks = c(0.1, 1, 10, 100, 1000),
                limits = c(0.1, 1000),
                labels = c("ND", 1, 10, 100, 1000)) +
  scale_x_discrete(labels=c('0', '10', '30\n        Study day', '40', '180', '190'),
                   breaks=c("a0", "a10", "a30", "a40", "a80", "a90")) +
  theme_classic() +
  facet_wrap(~factor(prev, levels = c("TRUE", "FALSE"))+person, nrow = 2, labeller = label_parsed, drop = FALSE) + #if want wide version
  annotation_logticks(side = "l" , size = .1) +
  scale_color_manual(values = prev_color) +
  theme(legend.position = "none") +
  geom_text(aes(label = paste("Participant", person), x = 3, y = 500), color = "black", 
            size = 2, family = "sans", fontface = "plain", check_overlap = TRUE) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) + 
  theme(axis.title = element_text(size = 6)) +
  theme(axis.text = element_text(size = 6)) +
  xlab("") + 
  ylab("Blood-detected clones at day 10 \n(clonotype copies)")
pdf(figure_2d_filename, width = 8.5, height = 3)
plot 
dev.off()

