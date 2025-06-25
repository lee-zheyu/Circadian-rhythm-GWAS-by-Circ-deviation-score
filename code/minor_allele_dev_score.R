library(ggplot2)
library(biomaRt)

set.seed(23)

load('Desktop/circ_gene/data/R_weighted_complete2.RData')
load('Desktop/circ_gene/data/genotype_group2.RData')

alt_allele_is_minor_allele <- vector('character', length = nrow(R_complete))
names(alt_allele_is_minor_allele) <- R_complete$RS_id
for (i in 1:length(alt_allele_is_minor_allele)) {
  rsid <- R_complete$RS_id[i]
  genotypes <- genotype_group[[rsid]]
  major_allele_freq <- length(genotypes[['0']])
  minor_allele_freq <- length(genotypes[['1']]) + length(genotypes[['2']])
  if (minor_allele_freq < major_allele_freq) {
    alt_allele_is_minor_allele[rsid] <- 'Y'
  } else {
    alt_allele_is_minor_allele[rsid] <- 'N'
  }
}

minor_allele_has_greater_dev_score <- vector('character', length = nrow(R_complete))
names(minor_allele_has_greater_dev_score) <- names(alt_allele_is_minor_allele)
for (i in 1:length(alt_allele_is_minor_allele)) {
  rsid <- R_complete$RS_id[i]
  regression_coef <- R_complete$coef[R_complete$RS_id == rsid]
  if (alt_allele_is_minor_allele[rsid] == 'Y' & regression_coef > 0) {
    minor_allele_has_greater_dev_score[rsid] <- 'Y'
  } else if (alt_allele_is_minor_allele[rsid] == 'N' & regression_coef < 0) {
    minor_allele_has_greater_dev_score[rsid] <- 'Y'
  } else {
    minor_allele_has_greater_dev_score[rsid] <- 'N'
  }
}

minor_allele_has_greater_dev_score_circ <- minor_allele_has_greater_dev_score

library(dplyr)
library(ggpubr)
library(nortest)

load('Desktop/circ_gene/data/random_SNP_weighted_complete.RData')
load('Desktop/circ_gene/data/deviation_by_tissue_weighted_complete.RData')

tissues <- unique(R$tissue)
tissue_corr <- vector('list', length = length(tissues))
names(tissue_corr) <- tissues
# for (i in 1:length(tissue_corr)) {
#   tissue <- tissues[i]
#   print(tissue)
#   deviation <- deviation_by_tissue[[tissue]]
#   donors <- rownames(deviation)
#   corr <- matrix(nrow = length(donors), ncol = length(donors))
#   rownames(corr) <- donors
#   colnames(corr) <- donors
#   for (j in 1:length(donors)) {
#     d1 <- donors[j]
#     d1_deviation <- deviation[rownames(deviation)==d1,]
#     for (k in j:length(donors)) {
#       d2 <- donors[k]
#       d2_deviation <- deviation[rownames(deviation)==d2,]
#       correlation <- cor(d1_deviation, d2_deviation)
#       corr[j,k] <- correlation
#       corr[k,j] <- correlation
#     }
#   }
#   tissue_corr[[i]] <- corr
# }
# 
# save('tissue_corr', file = 'Desktop/circ_gene/data/random_tissue_corr.RData')
load('Desktop/circ_gene/data2/random_tissue_corr.RData')

# gt_directory <- '~/cmb0/circadian_gene/result/indiv_gt'
gt_directory <- 'Desktop/circ_gene/data/indiv_gt2'
genotype_group <- vector('list', length = nrow(R))
names(genotype_group) <- R$variant_id
# for (i in 1:length(genotype_group)) {
#   print(i)
#   tissue <- R$tissue[i]
#   donor_id <- rownames(deviation_by_tissue[[tissue]])
#   ind_parts <- unlist(strsplit(R$ind[i],':'))
#   gt_file <- file.path(gt_directory, paste0('indiv_gt_frag_', ind_parts[1], '.RData'))
#   load(gt_file)
#   colnames(gt) <- unlist(lapply(colnames(gt), function(x) {
#     gsub('GTEX-','',x)
#   }))
#   donor_id_gt <- colnames(gt)
#   gt_sub <- gt[,donor_id_gt %in% donor_id]
#   gt_snp <- t(gt_sub[as.numeric(ind_parts[2]),])
# 
#   group0_donors <- rownames(gt_snp)[gt_snp==0]
#   group1_donors <- rownames(gt_snp)[gt_snp==1]
#   group2_donors <- rownames(gt_snp)[gt_snp==2]
# 
#   genotype_group[[i]] <- list('0'=group0_donors, '1'=group1_donors, '2'=group2_donors)
# }
# 
# save('genotype_group', file = 'Desktop/circ_gene/data/random_genotype_group.RData')
load('Desktop/circ_gene/data/random_genotype_group.RData')

alt_allele_is_minor_allele <- vector('character', length = nrow(R))
names(alt_allele_is_minor_allele) <- R$variant_id
for (i in 1:length(alt_allele_is_minor_allele)) {
  variant_id <- R$variant_id[i]
  genotypes <- genotype_group[[variant_id]]
  major_allele_freq <- length(genotypes[['0']])
  minor_allele_freq <- length(genotypes[['1']]) + length(genotypes[['2']])
  if (minor_allele_freq < major_allele_freq) {
    alt_allele_is_minor_allele[variant_id] <- 'Y'
  } else {
    alt_allele_is_minor_allele[variant_id] <- 'N'
  }
}

minor_allele_has_greater_dev_score <- vector('character', length = nrow(R))
names(minor_allele_has_greater_dev_score) <- names(alt_allele_is_minor_allele)
for (i in 1:length(alt_allele_is_minor_allele)) {
  variant_id <- R$variant_id[i]
  regression_coef <- R$coef[R$variant_id == variant_id]
  if (alt_allele_is_minor_allele[variant_id] == 'Y' & regression_coef > 0) {
    minor_allele_has_greater_dev_score[variant_id] <- 'Y'
  } else if (alt_allele_is_minor_allele[variant_id] == 'N' & regression_coef < 0) {
    minor_allele_has_greater_dev_score[variant_id] <- 'Y'
  } else {
    minor_allele_has_greater_dev_score[variant_id] <- 'N'
  }
}

minor_allele_has_greater_dev_score_random <- minor_allele_has_greater_dev_score
minor_allele_has_greater_dev_score_random <- minor_allele_has_greater_dev_score_random[
  sample(x = seq(1, length(minor_allele_has_greater_dev_score_random)), size = 2000)]

circ_tot <- length(minor_allele_has_greater_dev_score_circ)
Y_freq_circ <- sum(minor_allele_has_greater_dev_score_circ == 'Y')
N_freq_circ <- 0
Y_perc_circ <- Y_freq_circ / circ_tot
N_perc_circ <- N_freq_circ / circ_tot

perc_df_circ <- data.frame(frequency = c(Y_freq_circ, N_freq_circ),
                           percentage = c(Y_perc_circ, N_perc_circ),
                           result = c('Minor allele has greater deviation score', 
                                      'Minor allele has smaller deviation score'),
                           group = 'Circ-SNPs')

random_tot <- length(minor_allele_has_greater_dev_score_random)
Y_freq_random <- sum(minor_allele_has_greater_dev_score_random == 'Y')
N_freq_random <- sum(minor_allele_has_greater_dev_score_random == 'N')
Y_perc_random <- Y_freq_random / random_tot
N_perc_random <- N_freq_random / random_tot

perc_df_random <- data.frame(frequency = c(Y_freq_random, N_freq_random),
                             percentage = c(Y_perc_random, N_perc_random),
                             result = c('Minor allele has greater deviation score', 
                                        'Minor allele has smaller deviation score'),
                             group = 'random SNPs')

perc_df <- rbind(perc_df_circ, perc_df_random)
perc_df$result <- factor(perc_df$result, levels = c('Minor allele has greater deviation score', 
                                                        'Minor allele has smaller deviation score'))

p <- ggplot(perc_df, aes(x=group, y=percentage, fill=result)) +
  geom_bar(stat = "identity", position=position_dodge(.9), color='black') +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        legend.text = element_text(color='black', size = 8),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', size = 8),
        axis.text.y = element_text(color='black', size = 8)) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("#00AFBB", "#FC4E07")) +
  geom_text(aes(y=percentage + 0.05, group=result, label=frequency),
            position=position_dodge(.9), size = 3) +
  guides(fill = guide_legend(nrow = 2)) 
p

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'f2_minor_allele_dev_score.png',
  plot = last_plot(),
  device = "png",  
  width = 3,
  height = 2.5,
  units = 'in',
  dpi = 300,
  bg = NULL,
)







