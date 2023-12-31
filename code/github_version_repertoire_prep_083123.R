rm(list=ls())

library(readr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)

#load and process
#day 0
ap1 <- read_tsv("file1.tsv")
#day 10 
bp1 <- read_tsv("file2.tsv")
#day 30 
cp1 <- read_tsv("file3.tsv")
#day 40 
dp1 <- read_tsv("file4.tsv")
#day 180 
ep1 <- read_tsv("file5.tsv")
#day 190 
fp1 <- read_tsv("file6.tsv")

##controls
#day 10 
gp1 <- read_tsv("file7.tsv")
#day 40 
hp1 <- read_tsv("file8.tsv")
#day 190 
ip1 <- read_tsv("file9.tsv")

##blood
#day 0 - IFNg
jp1 <- read_tsv("file10.tsv")
#day 10 - IFNg 
kp1 <- read_tsv("file11.tsv")

# background
lp1 <- read_tsv("file12.tsv")
mp1 <- read_tsv("file13.tsv")

##add in location and time codes 
ap1$loc <- "a0"
bp1$loc <- "a10"
cp1$loc <- "a30"
dp1$loc <- "a40"
ep1$loc <- "a80"
fp1$loc <- "a90"
gp1$loc <-  "c10"
hp1$loc <-  "c40"
ip1$loc <-  "c90"
jp1$loc <-  "g0"
kp1$loc <-  "g10"
lp1$loc <-  "c0"
mp1$loc <-  "les" 

##rbind 
p1 <- rbind(ap1, bp1, cp1, dp1, ep1, fp1, gp1, hp1, ip1, jp1, kp1, lp1, mp1)
remove(ap1, bp1, cp1, dp1, ep1, fp1, gp1, hp1, ip1, jp1, kp1, lp1, mp1)

## p1
p1 <- drop_na(p1, aminoAcid)
p1 <- p1[!grepl('[*]', p1$aminoAcid),]
p1 <- p1[, -c(6,9:13,16:20,23:52)]

p1[ p1 == "unresolved" ] <- NA
p1$v <- ifelse(is.na(p1$vGeneName), p1$vFamilyName, p1$vGeneName)
p1$d <- ifelse(is.na(p1$dGeneName), p1$dFamilyName, p1$dGeneName)
p1$j <- ifelse(is.na(p1$jGeneName), p1$jFamilyName, p1$jGeneName)
p1 <- p1[, -c(6:11)]

p1c <- p1[, -4] #make a count spreadsheet
p1f <- p1[, -3] #make a frequency spreadsheet
colnames(p1c) <- c("nucleotide", "aminoAcid", "count", "id", "v", "j") 
colnames(p1f) <- c("nucleotide", "aminoAcid", "freq", "id", "v", "j") 

p1c <- dplyr::select(p1c, nucleotide:j)
p1c <- tidyr::spread(p1c, id, count)

p1f <- dplyr::select(p1f, nucleotide:j)
p1f <- tidyr::spread(p1f, id, freq)

# consolidate - nucleotide and amino acid versions, by count (c) and by frequency (f)
AAc <- p1c %>% 
  group_by(aminoAcid, v, j) %>%
  summarise(g0 = sum(g0[!is.na(g0)]), g10 = sum(g10[!is.na(g10)]), 
            a0 = sum(a0[!is.na(a0)]), a10 = sum(a10[!is.na(a10)]), 
            a30 = sum(a30[!is.na(a30)]), a40 = sum(a40[!is.na(a40)]), 
            a80 = sum(a80[!is.na(a80)]), a90 = sum(a90[!is.na(a90)]), 
            c0 = sum(c0[!is.na(c0)]), les = sum(les[!is.na(les)]),
            c10 = sum(c10[!is.na(c10)]), c40 = sum(c40[!is.na(c40)]), 
            c90 = sum(c90[!is.na(c90)]))

NCc <- p1c %>% 
  group_by(nucleotide, aminoAcid, v, j) %>%
  summarise(g0 = sum(g0[!is.na(g0)]), g10 = sum(g10[!is.na(g10)]), 
            a0 = sum(a0[!is.na(a0)]), a10 = sum(a10[!is.na(a10)]), 
            a30 = sum(a30[!is.na(a30)]), a40 = sum(a40[!is.na(a40)]), 
            a80 = sum(a80[!is.na(a80)]), a90 = sum(a90[!is.na(a90)]), 
            c0 = sum(c0[!is.na(c0)]), les = sum(les[!is.na(les)]),
            c10 = sum(c10[!is.na(c10)]), c40 = sum(c40[!is.na(c40)]), 
            c90 = sum(c90[!is.na(c90)]))

AAf <- p1f %>% 
  group_by(aminoAcid, v, j) %>%
  summarise(g0 = sum(g0[!is.na(g0)]), g10 = sum(g10[!is.na(g10)]), 
            a0 = sum(a0[!is.na(a0)]), a10 = sum(a10[!is.na(a10)]), 
            a30 = sum(a30[!is.na(a30)]), a40 = sum(a40[!is.na(a40)]), 
            a80 = sum(a80[!is.na(a80)]), a90 = sum(a90[!is.na(a90)]), 
            c0 = sum(c0[!is.na(c0)]), les = sum(les[!is.na(les)]),
            c10 = sum(c10[!is.na(c10)]), c40 = sum(c40[!is.na(c40)]), 
            c90 = sum(c90[!is.na(c90)]))

NCf <- p1f %>% 
  group_by(nucleotide, aminoAcid, v, j) %>%
  summarise(g0 = sum(g0[!is.na(g0)]), g10 = sum(g10[!is.na(g10)]), 
            a0 = sum(a0[!is.na(a0)]), a10 = sum(a10[!is.na(a10)]), 
            a30 = sum(a30[!is.na(a30)]), a40 = sum(a40[!is.na(a40)]), 
            a80 = sum(a80[!is.na(a80)]), a90 = sum(a90[!is.na(a90)]), 
            c0 = sum(c0[!is.na(c0)]), les = sum(les[!is.na(les)]),
            c10 = sum(c10[!is.na(c10)]), c40 = sum(c40[!is.na(c40)]), 
            c90 = sum(c90[!is.na(c90)]))

AAc$person <- "p1"
NCc$person <- "p1"
AAf$person <- "p1"
NCf$person <- "p1"
remove(p1, p1c, p1f)