library(dplyr)
library(ggpubr)
library(nortest)

load('Desktop/circ_gene/data/R_weighted_complete2.RData')
load('Desktop/circ_gene/data/normalized_tpm_by_tissue.RData')

tissues <- unique(R_complete$tissue)
tissue_corr <- vector('list', length = length(tissues))
names(tissue_corr) <- tissues
for (i in 1:length(tissue_corr)) {
  tissue <- tissues[i]
  print(tissue)
  normalized_tpm <- normalized_tpm_by_tissue[[tissue]]
  donors <- rownames(normalized_tpm)
  corr <- matrix(nrow = length(donors), ncol = length(donors))
  rownames(corr) <- donors
  colnames(corr) <- donors
  for (j in 1:length(donors)) {
    d1 <- donors[j]
    d1_normalized_tpm <- normalized_tpm[rownames(normalized_tpm)==d1,]
    for (k in j:length(donors)) {
      d2 <- donors[k]
      d2_normalized_tpm <- normalized_tpm[rownames(normalized_tpm)==d2,]
      correlation <- cor(d1_normalized_tpm, d2_normalized_tpm, method = 'spearman')
      corr[j,k] <- correlation
      corr[k,j] <- correlation
    }
  }
  tissue_corr[[i]] <- corr
}

save('tissue_corr', file = 'Desktop/circ_gene/data/tissue_corr_normalized_transformed_tpm_spearman.RData')
load('Desktop/circ_gene/data/tissue_corr_transformed_spearman.RData')

gt_directory <- 'Desktop/circ_gene/data/indiv_gt2'
genotype_group <- vector('list', length = nrow(R_complete))
names(genotype_group) <- R_complete$RS_id
for (i in 1:length(genotype_group)) {
  print(i)
  tissue <- R_complete$tissue[i]
  donor_id <- rownames(normalized_tpm_by_tissue[[tissue]])
  ind_parts <- unlist(strsplit(R_complete$ind[i],':'))
  gt_file <- file.path(gt_directory, paste0('indiv_gt_frag_', ind_parts[1], '.RData'))
  load(gt_file)
  colnames(gt) <- unlist(lapply(colnames(gt), function(x) {
    gsub('GTEX-','',x)
  }))
  donor_id_gt <- colnames(gt)
  gt_sub <- gt[,donor_id_gt %in% donor_id]
  gt_snp <- t(gt_sub[as.numeric(ind_parts[2]),])

  group0_donors <- rownames(gt_snp)[gt_snp==0]
  group1_donors <- rownames(gt_snp)[gt_snp==1]
  group2_donors <- rownames(gt_snp)[gt_snp==2]

  genotype_group[[i]] <- list('0'=group0_donors, '1'=group1_donors, '2'=group2_donors)
}

load('Desktop/circ_gene/data/genotype_group2.RData')

corr_by_genotype <- vector('list', length = length(genotype_group))
names(corr_by_genotype) <- names(genotype_group)
for (i in 1:length(corr_by_genotype)) {
  print(i)
  tissue <- R_complete$tissue[i]
  corr <- tissue_corr[[tissue]]
  genotype <- genotype_group[[i]]

  group0_donors <- genotype[['0']]
  if (length(group0_donors) > 1) {
    group0_combination_num <- choose(length(group0_donors), 2)
    group0_donor_combo <- vector('character', length = group0_combination_num)
    group0_corr <- vector('numeric', length = group0_combination_num)
    group0_element_ind <- 1
    for (j in 1:length(group0_donors)) {
      d1 <- group0_donors[j]
      row_ind <- seq(1,nrow(corr))[rownames(corr)==d1]
      for (k in j:length(group0_donors)) {
        if (j==k) {
          next
        }
        d2 <- group0_donors[k]
        col_ind <- seq(1,ncol(corr))[colnames(corr)==d2]
        group0_donor_combo[group0_element_ind] <- paste(d1,'-', d2)
        group0_corr[group0_element_ind] <- corr[row_ind, col_ind]
        group0_element_ind <- group0_element_ind + 1
      }
    }
    names(group0_corr) <- group0_donor_combo
  } else {
    group0_corr <- c()
  }

  group1_donors <- genotype[['1']]
  if (length(group1_donors) > 1) {
    group1_combination_num <- choose(length(group1_donors), 2)
    group1_donor_combo <- vector('character', length = group1_combination_num)
    group1_corr <- vector('numeric', length = group1_combination_num)
    group1_element_ind <- 1
    for (j in 1:length(group1_donors)) {
      d1 <- group1_donors[j]
      row_ind <- seq(1,nrow(corr))[rownames(corr)==d1]
      for (k in j:length(group1_donors)) {
        if (j==k) {
          next
        }
        d2 <- group1_donors[k]
        col_ind <- seq(1,ncol(corr))[colnames(corr)==d2]
        group1_donor_combo[group1_element_ind] <- paste(d1,'-', d2)
        group1_corr[group1_element_ind] <- corr[row_ind, col_ind]
        group1_element_ind <- group1_element_ind + 1
      }
    }
    names(group1_corr) <- group1_donor_combo
  } else {
    group1_corr <- c()
  }

  group2_donors <- genotype[['2']]
  if (length(group2_donors) > 1) {
    group2_combination_num <- choose(length(group2_donors), 2)
    group2_donor_combo <- vector('character', length = group2_combination_num)
    group2_corr <- vector('numeric', length = group2_combination_num)
    group2_element_ind <- 1
    for (j in 1:length(group2_donors)) {
      d1 <- group2_donors[j]
      row_ind <- seq(1,nrow(corr))[rownames(corr)==d1]
      for (k in j:length(group2_donors)) {
        if (j==k) {
          next
        }
        d2 <- group2_donors[k]
        col_ind <- seq(1,ncol(corr))[colnames(corr)==d2]
        group2_donor_combo[group2_element_ind] <- paste(d1,'-', d2)
        group2_corr[group2_element_ind] <- corr[row_ind, col_ind]
        group2_element_ind <- group2_element_ind + 1
      }
    }
    names(group2_corr) <- group2_donor_combo
  } else {
    group2_corr <- c()
  }

  if (length(group0_donors) > 0 & length(group1_donors) > 0) {
    group01_combination_num <- length(group0_donors) * length(group1_donors)
    group01_donor_combo <- vector('character', length = group01_combination_num)
    group01_corr <- vector('numeric', length = group01_combination_num)
    group01_element_ind <- 1
    for (j in 1:length(group0_donors)) {
      d1 <- group0_donors[j]
      row_ind <- seq(1,nrow(corr))[rownames(corr)==d1]
      for (k in 1:length(group1_donors)) {
        d2 <- group1_donors[k]
        col_ind <- seq(1,ncol(corr))[colnames(corr)==d2]
        group01_donor_combo[group01_element_ind] <- paste(d1,'-', d2)
        group01_corr[group01_element_ind] <- corr[row_ind, col_ind]
        group01_element_ind <- group01_element_ind + 1
      }
    }
    names(group01_corr) <- group01_donor_combo
  } else {
    group01_corr <- c()
  }

  if (length(group0_donors) > 0 & length(group2_donors) > 0) {
    group02_combination_num <- length(group0_donors) * length(group2_donors)
    group02_donor_combo <- vector('character', length = group02_combination_num)
    group02_corr <- vector('numeric', length = group02_combination_num)
    group02_element_ind <- 1
    for (j in 1:length(group0_donors)) {
      d1 <- group0_donors[j]
      row_ind <- seq(1,nrow(corr))[rownames(corr)==d1]
      for (k in 1:length(group2_donors)) {
        d2 <- group2_donors[k]
        col_ind <- seq(1,ncol(corr))[colnames(corr)==d2]
        group02_donor_combo[group02_element_ind] <- paste(d1,'-', d2)
        group02_corr[group02_element_ind] <- corr[row_ind, col_ind]
        group02_element_ind <- group02_element_ind + 1
      }
    }
    names(group02_corr) <- group02_donor_combo
  } else {
    group02_corr <- c()
  }

  if (length(group1_donors) > 0 & length(group2_donors) > 0) {
    group12_combination_num <- length(group1_donors) * length(group2_donors)
    group12_donor_combo <- vector('character', length = group12_combination_num)
    group12_corr <- vector('numeric', length = group12_combination_num)
    group12_element_ind <- 1
    for (j in 1:length(group1_donors)) {
      d1 <- group1_donors[j]
      row_ind <- seq(1,nrow(corr))[rownames(corr)==d1]
      for (k in 1:length(group2_donors)) {
        d2 <- group2_donors[k]
        col_ind <- seq(1,ncol(corr))[colnames(corr)==d2]
        group12_donor_combo[group12_element_ind] <- paste(d1,'-', d2)
        group12_corr[group12_element_ind] <- corr[row_ind, col_ind]
        group12_element_ind <- group12_element_ind + 1
      }
    }
    names(group12_corr) <- group12_donor_combo
  } else {
    group12_corr <- c()
  }

  corr_by_genotype[[i]] <- list('A' = c(group0_corr, group1_corr, group2_corr),
                                'B' = c(group01_corr, group02_corr, group12_corr))
}

save('corr_by_genotype', file = 'Desktop/circ_gene/data/corr_by_genotype_normalized_tpm_spearman.RData')
load('Desktop/circ_gene/data/corr_by_genotype_normalized_tpm_spearman.RData')

wilcoxon_result <- data.frame(matrix(nrow = length(corr_by_genotype), ncol = 5))
colnames(wilcoxon_result) <- c('rs_id', 'tissue', 'sample_size','p_val', 'result')
p_val_threshold <- .05
bad_ind <- c()
for (i in 1:nrow(wilcoxon_result)) {
  print(i)
  snp_corr <- corr_by_genotype[[i]]
  
  wilcoxon_result$rs_id[i] <- R_complete$RS_id[i]
  tissue <- R_complete$tissue[i]
  wilcoxon_result$tissue[i] <- tissue
  wilcoxon_result$sample_size[i] <- nrow(tissue_corr[[tissue]])
  
  group_A_corr <- snp_corr[['A']]
  group_A_exists <- !is.null(group_A_corr)
  
  group_B_corr <- snp_corr[['B']]
  group_B_exists <- !is.null(group_B_corr)
  
  if (sum(group_A_exists, group_B_exists) < 2) {
    bad_ind <- c(bad_ind, i)
    next
  }
  
  group_AB_result <- wilcox.test(group_A_corr, group_B_corr, alternative = 'greater', conf.level = .95)
  p_val <- group_AB_result$p.value
  wilcoxon_result$p_val[i] <- p_val
  if (p_val <= p_val_threshold) {
    wilcoxon_result$result[i] <- 'bad'
  } else {
    wilcoxon_result$result[i] <- 'good'
  }
}
good_rsid <- wilcoxon_result$rs_id[wilcoxon_result$result=='good']

unsynced_rsid <- good_rsid
length(unsynced_rsid)/nrow(wilcoxon_result)

