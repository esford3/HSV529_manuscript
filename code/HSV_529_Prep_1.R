

# code/HSV_529_Prep_1.R

library(ggplot2)
library(dplyr)
library(tidyr)


# INPUTS
user = 'esf'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software/HSV529_manuscript/'
} else {
  stop("set repo loc and repo manually")
}

filename = file.path(repo_loc, '/data/vax_nt_all_dei.csv')
vax <- read_csv(filename)

g <- filter(vax, g0>0 & !is.na(g0) | g10>0 & !is.na(g10))
g <- select(g, c("nucleotide", "aminoAcid", "v", "j", "g0", "g10", "les", "a0",  "a10", "a30", "a40", "a80", "a90",
                 "c0", "c10", "c40", "c90", "person"))
g <- mutate(g, person = factor(person))
colnames(g) <- c("cdr3_b_nt", "cdr3_b_aa", "v_b_gene", "j_b_gene", "blood_d0", "blood_d10", "lesion", "skin_d0",  "skin_d10", "skin_d30", "skin_d40", "skin_d180", "skin_d190",
                 "arm_d0", "arm_d10", "arm_d40", "arm_d190", "person")
#dd[with(dd, order(-z, b)), ] #https://stackoverflow.com/questions/1296646/sort-order-data-frame-rows-by-multiple-columns 
g <- g[with(g, order(person, -blood_d0, -blood_d10, -lesion, -skin_d0, -skin_d10)),]
write.csv(g, "blood-identified clonotypes 9.13.23.csv", row.names = F)
prev <- filter(g, !is.na(skin_d0) & skin_d0>0)
write.csv(prev, "sup tab 4a prev blood-identified clonotypes 9.13.23.csv", row.names = F)
elic <- filter(g, skin_d0==0 & (skin_d10>0 | skin_d30>0 | skin_d40>0))
write.csv(elic, "sup tab 4b elic blood-identified clonotypes 9.13.23.csv", row.names = F)
lesg <- filter(g, !is.na(lesion) & lesion>0)
write.csv(lesg, "sup tab 4c lesion blood-identified clonotypes 9.13.23.csv", row.names = F)


## Table prep - supplemental table 6
d <- select(vax, c("nucleotide", "aminoAcid", "v", "j", "g0", "g10", "les", "a0",  "a10", "a30", "a40", "a80", "a90",
                          "c0", "c10", "c40", "c90", "person"))
d <- mutate(d, person = factor(person))
colnames(d) <- c("cdr3_b_nt", "cdr3_b_aa", "v_b_gene", "j_b_gene", "blood_d0", "blood_d10", "lesion", "skin_d0",  "skin_d10", "skin_d30", "skin_d40", "skin_d180", "skin_d190",
                 "arm_d0", "arm_d10", "arm_d40", "arm_d190", "person")
d <- d[with(d, order(person, -blood_d0, -blood_d10, -lesion, -skin_d0, -skin_d10)),]
prev <- filter(d, skin_d0>0)
elic <- filter(d, skin_d0==0 & (skin_d10>0 | skin_d30>0 | skin_d40>0 | (!is.na(skin_d180) & skin_d180>0) | (!is.na(skin_d190) & skin_d190>0)))
prev$fold_d1 <- prev$skin_d10/prev$skin_d0
prev_d1 <- filter(prev, fold_d1 >=6)
prev_d1 <- prev_d1[with(prev_d1, order(person, -fold_d1)),]
write.csv(prev_d1, "sup table 6a prev_exp clonotypes 9.13.23.csv", row.names = F)

prev$fold_d2 <- prev$skin_d40/prev$skin_d30
prev_d2 <- filter(prev, fold_d2 >=6 & fold_d2!="Inf")
prev_d2 <- prev_d2[with(prev_d2, order(person, -fold_d2)),]
write.csv(prev_d2, "sup table 6c prev_exp_d2 clonotypes 9.13.23.csv", row.names = F)

elic$fold_d1 <- elic$skin_d10
elic_d1 <- filter(elic, fold_d1 >=6)
elic_d1 <- elic_d1[with(elic_d1, order(person, -fold_d1)),]
write.csv(elic_d1, "sup table 6b elic_exp clonotypes 9.13.23.csv", row.names = F)

elic$fold_d2 <- elic$skin_d40/elic$skin_d30
elic_d2 <- filter(elic, fold_d2 >=6 & fold_d2!="Inf")
elic_d2 <- elic_d2[with(elic_d2, order(person, -fold_d2)),]
write.csv(elic_d2, "sup table 6c elic_exp_d2 clonotypes 9.13.23.csv", row.names = F)
