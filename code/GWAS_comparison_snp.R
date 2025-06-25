# GWAS comparison for complete SNP list

library(rsnps)
library(ggplot2)
library(ggtext)
library(tidyverse)
library("writexl")
source('Desktop/circ_gene/data/diseaseSNPsFunction.R')
load('Desktop/circ_gene/data/R_weighted_complete2.RData')


# cleanup GWAS catalog
GWAS_catalog <- read.delim('Desktop/circ_gene/data/gwas-catalog-associations_ontology-annotated.tsv',
                           quote="")

GWAS_clean_1 <- GWAS_catalog[GWAS_catalog$SNP_ID_CURRENT!='',]
GWAS_clean_2 <- GWAS_clean_1[GWAS_clean_1$MAPPED_TRAIT!='',]
GWAS_clean_3 <- GWAS_clean_2[!duplicated(GWAS_clean_2),]

record_ind_w_missing_info <- seq(1,nrow(GWAS_clean_3))[GWAS_clean_3$CHR_ID=='' |
                                                             GWAS_clean_3$CHR_POS=='']

GWAS_missing_info <- GWAS_clean_3[record_ind_w_missing_info,]
dbSNP_status <- data.frame(matrix(nrow = length(record_ind_w_missing_info), ncol = 3))
colnames(dbSNP_status) <- c('ind','rsid','status')
for (i in 1:length(record_ind_w_missing_info)) {
  print(i)
  ind <- record_ind_w_missing_info[i]
  dbSNP_status$ind[i] <- ind
  rsid <- paste0('rs',GWAS_clean_3$SNP_ID_CURRENT[ind])
  dbSNP_status$rsid[i] <- rsid
  tryCatch(
    expr = {
      snp_record <- ncbi_snp_query(rsid)
      if (snp_record$class=='snv') {
        if (nrow(snp_record)!=0) {
          GWAS_clean_3$CHR_ID[ind] <- snp_record$chromosome
          GWAS_clean_3$CHR_POS[ind] <- snp_record$bp
          dbSNP_status$status[i] <- 'found'
        }
      } else {
        dbSNP_status$status[i] <- 'not snv'
      }
    },
    error = function(e) {
      dbSNP_status$status[i] <<- e$message
    },
    warning = function(w) {
      # W <<- w
      dbSNP_status$status[i] <<- w$message
    }
  )
}

save(dbSNP_status, file = 'Desktop/circ_gene/data/dbSNP_status_complete_v2.RData')

GWAS_clean_4 <- GWAS_clean_3[-dbSNP_status$ind[dbSNP_status$status!='found'],]
GWAS_clean_4$CHR_POS <- as.numeric(GWAS_clean_4$CHR_POS)
GWAS_clean_4$SNP_ID_CURRENT <- paste0('rs',GWAS_clean_4$SNP_ID_CURRENT)

save(GWAS_clean_4, file = 'Desktop/circ_gene/data/gwas_snp_clean_4_complete_v2.RData')

load('Desktop/circ_gene/data/gwas_snp_clean_4_complete_v2.RData')

# tag GWAS SNPs by significant SNPs
GWAS_chr <- paste0('chr',GWAS_clean_4$CHR_ID)
GWAS_range_up <- GWAS_clean_4$CHR_POS
GWAS_range_down <- GWAS_clean_4$CHR_POS
GWAS_snp_range <- data.frame(chr = GWAS_chr,
                             start = GWAS_range_up,
                             end = GWAS_range_down)

GWAS_Granges <- range2GRanges(GWAS_snp_range)
names(GWAS_Granges) <- GWAS_clean_4$SNP_ID_CURRENT

hit_GWAS_snp <- vector('list',length = nrow(R_complete))
names(hit_GWAS_snp) <- R_complete$RS_id
missing_ind <- c()
for (i in 1:length(hit_GWAS_snp)) {
  print(i)

  # convert snp coor to GRange
  snpCoor <- R_complete$pos[i]
  snpsRanges <- string2range(snpCoor, delim=":", region=FALSE)
  snpsGranges <- range2GRanges(snpsRanges)
  names(snpsGranges) <- snpCoor

# find overlap of snps and genes
  overlap <- GenomicRanges::findOverlaps(snpsGranges, GWAS_Granges, maxgap = 1000000, type = 'any')

  # make vector of SNPs to genes
  hit_snp <- unique(names(GWAS_Granges[overlap@to]))

  if (length(hit_snp)==0) {
    missing_ind <- c(missing_ind,i)
    next
  }

  traits <- vector('list', length = length(hit_snp))
  names(traits) <- hit_snp
  for (j in 1:length(hit_snp)) {
    rsid <- hit_snp[j]
    traits[[j]] <- unique(GWAS_clean_4$MAPPED_TRAIT[GWAS_clean_4$SNP_ID_CURRENT %in% rsid])
  }

  hit_GWAS_snp[[i]] <- list(rsids=hit_snp,traits=traits)
}
save(list=c('hit_GWAS_snp', 'missing_ind'),
     file = 'Desktop/circ_gene/data/hit_GWAS_SNP_weighted_complete_v2.RData')

load('Desktop/circ_gene/data/hit_GWAS_SNP_weighted_complete_v2.RData')

if (length(missing_ind)==0) {
  hit_GWAS_snp_clean <- hit_GWAS_snp
} else {
  hit_GWAS_snp_clean <- hit_GWAS_snp[-missing_ind]
}

target_traits <- c('sleep latency', 'sleep depth', 'sleep apnea', 'Sleep efficiency', 'sleep quality',
                   'sleep measurement', 'sleep time', 'sleep duration', 'sleep apnea measurement during REM sleep',
                   'sleep apnea measurement during non-REM sleep', 'short sleep', 'obstructive sleep apnea',
                   'sleep disorder', 'Wake After Sleep Onset', 'REM sleep behavior disorder', 'sleep apnea measurement',
                   'nighttime rest measurement', 'bruxism', 'daytime rest measurement', 'narcolepsy',
                   'insomnia measurement', 'Somnambulism', 'snoring measurement', 'periodic limb movement disorder',
                   'narcolepsy without cataplexy', 'narcolepsy-cataplexy syndrome', 'hypersomnia',
                   'excessive daytime sleepiness measurement', 'Kleine-Levin syndrome', 'insomnia',
                   'circadian rhythm', 'sleepiness', 'chronotype measurement', 'manic episode measurement',
                   'irritability measurement', 'hypertension', 'low density lipoprotein cholesterol measurement',
                   'depressive symptom measurement', 'triglyceride measurement', 'high density lipoprotein cholesterol measurement',
                   'coronary artery disease', 'emotional symptom measurement')

target_traits <- target_traits[target_traits %in% GWAS_clean_4$MAPPED_TRAIT]

# comprare tagged GWAS SNPs to all GWAS SNPs
hit_GWAS_snp_unlisted <- c()
for (i in 1:length(hit_GWAS_snp_clean)) {
  hit_GWAS_snp_unlisted <- union(hit_GWAS_snp_unlisted,hit_GWAS_snp_clean[[i]]$rsids)
}
length(hit_GWAS_snp_unlisted)

# all target traits
GWAS_target_trait_ind <- seq(1,nrow(GWAS_clean_4))[GWAS_clean_4$MAPPED_TRAIT %in% target_traits]
all_snp_target_trait <- unique(GWAS_clean_4$SNP_ID_CURRENT[GWAS_target_trait_ind])
all_snp_non_target_trait <- unique(GWAS_clean_4$SNP_ID_CURRENT[
  !GWAS_clean_4$SNP_ID_CURRENT %in% all_snp_target_trait])
all_gwas_snp <- union(all_snp_target_trait,all_snp_non_target_trait)

hit_snp_target_trait <- hit_GWAS_snp_unlisted[hit_GWAS_snp_unlisted %in% all_snp_target_trait]
hit_snp_non_target_trait <- hit_GWAS_snp_unlisted[hit_GWAS_snp_unlisted %in% all_snp_non_target_trait]

h_pval_by_snp <- phyper(length(hit_snp_target_trait) - 1, # drawn white ball
                        length(all_snp_target_trait), # total white ball
                        length(all_snp_non_target_trait), # total black ball
                        length(hit_GWAS_snp_unlisted), # num of drawn
                        lower.tail = FALSE)

hit_snp_prop <- prop.test(length(hit_snp_target_trait),length(hit_GWAS_snp_unlisted))
all_snp_prop <- prop.test(length(all_snp_target_trait),length(all_gwas_snp))
hit_snp_prop$conf.int[1] > all_snp_prop$conf.int[2]

prop_test_by_snp <- vector('list', length = length(target_traits))
names(prop_test_by_snp) <- target_traits
enriched_trait_ind <- c()
conf_int_diff <- c()
h_test_by_snp <- vector('numeric', length = length(target_traits))
names(h_test_by_snp) <- target_traits
for (i in 1:length(target_traits)) {
  GWAS_target_trait_ind <- seq(1,nrow(GWAS_clean_4))[GWAS_clean_4$MAPPED_TRAIT %in% target_traits[i]]
  all_snp_target_trait <- unique(GWAS_clean_4$SNP_ID_CURRENT[GWAS_target_trait_ind])
  all_snp_non_target_trait <- unique(GWAS_clean_4$SNP_ID_CURRENT[
    !GWAS_clean_4$SNP_ID_CURRENT %in% all_snp_target_trait])
  all_gwas_snp <- union(all_snp_target_trait,all_snp_non_target_trait)

  hit_snp_target_trait <- hit_GWAS_snp_unlisted[hit_GWAS_snp_unlisted %in% all_snp_target_trait]
  hit_snp_non_target_trait <- hit_GWAS_snp_unlisted[hit_GWAS_snp_unlisted %in% all_snp_non_target_trait]

  h_test_by_snp[i] <- phyper(length(hit_snp_target_trait) - 1, # drawn white ball
                             length(all_snp_target_trait), # total white ball
                             length(all_snp_non_target_trait), # total black ball
                             length(hit_GWAS_snp_unlisted), # num of drawn
                             lower.tail = FALSE)

  hit_GWAS_prop <- prop.test(length(hit_snp_target_trait),length(hit_GWAS_snp_unlisted))
  all_GWAS_prop <- prop.test(length(all_snp_target_trait),length(all_gwas_snp))
  prop_test_by_snp[[i]] <- list(hit_prop=hit_GWAS_prop,all_prop=all_GWAS_prop)
  if (hit_GWAS_prop$conf.int[1] > all_GWAS_prop$conf.int[2]) {
    enriched_trait_ind <- c(enriched_trait_ind,i)
    conf_int_diff <- c(conf_int_diff,hit_GWAS_prop$conf.int[1] - all_GWAS_prop$conf.int[2])
  }
}

enriched_traits_by_snp_hyper <- h_test_by_snp[h_test_by_snp < .05]
enriched_traits_by_snp_proportion <- target_traits[enriched_trait_ind]

snp_prop_test_result_df <- c()
for (i in 1:length(prop_test_by_snp)) {
  trait_df <- data.frame(matrix(nrow = 2, ncol = 5))
  colnames(trait_df) <- c('trait', 'percentage', 'lower_lim', 'upper_lim', 'group')
  trait <- target_traits[i]
  trait_df$trait <- trait
  test_result <- prop_test_by_snp[[trait]]

  trait_df$group[1] <- 'hit'
  trait_df$percentage[1] <- test_result$hit_prop$estimate
  trait_df$lower_lim[1] <- test_result$hit_prop$conf.int[1]
  trait_df$upper_lim[1] <- test_result$hit_prop$conf.int[2]

  trait_df$group[2] <- 'overall'
  trait_df$percentage[2] <- test_result$all_prop$estimate
  trait_df$lower_lim[2] <- test_result$all_prop$conf.int[1]
  trait_df$upper_lim[2] <- test_result$all_prop$conf.int[2]

  snp_prop_test_result_df <- rbind(snp_prop_test_result_df,trait_df)
}

hit_df <- snp_prop_test_result_df[snp_prop_test_result_df$group=='hit',]
trait_ordered <- hit_df$trait[order(hit_df$percentage, decreasing = F)]
snp_prop_test_result_df$trait <- factor(snp_prop_test_result_df$trait, levels = trait_ordered)

hit_GWAS_trait_unlisted <- c()
for (i in 1:length(hit_GWAS_snp_clean)) {
  hit_GWAS_trait_unlisted <- union(hit_GWAS_trait_unlisted,unlist(hit_GWAS_snp_clean[[i]]$traits))
}

all_GWAS_traits <- unique(GWAS_clean_4$MAPPED_TRAIT)
length(all_GWAS_traits)
hit_target_traits <- hit_GWAS_trait_unlisted[hit_GWAS_trait_unlisted %in% target_traits]
non_target_traits <- all_GWAS_traits[!all_GWAS_traits %in% target_traits]

hyper_p_val_trait <- phyper(length(hit_target_traits) - 1, # drawn white ball
                            length(target_traits), # total white ball
                            length(non_target_traits), # total black ball
                            length(hit_GWAS_trait_unlisted), # num of drawn
                            lower.tail = FALSE)

target_traits_df <- data.frame(trait = target_traits, 
                               contain_proxy_SNPs = target_traits %in% hit_target_traits)
write_xlsx(target_traits_df,"Desktop/circ_gene/data/hit_target_traits_v2.xlsx")

hit_gene_id <- c()
for (i in 1:nrow(R_complete)) {
  genes <- unlist(strsplit(R_complete$Gene[i], ';'))
  hit_gene_id <- union(hit_gene_id,genes)
}

all_GWAS_genes_list <- vector('list', length = nrow(GWAS_clean_4))
for (i in 1:length(all_GWAS_genes_list)) {
  genes <- unlist(strsplit(GWAS_clean_4$SNP_GENE_IDS[i],','))
  if (length(genes)==1) {
    all_GWAS_genes_list[[i]] <- genes
  } else if (length(genes)>1) {
    all_GWAS_genes_list[[i]] <- unlist(lapply(genes, function(x) {gsub(' ','',x)}))
  } else {
    next
  }
}

gtf <- read.table('Desktop/circ_gene/data/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz', 
                  header = F, stringsAsFactors = F, sep = '\t')
gtf_gene <- gtf[gtf[,3]=="gene", ]
gene_type <- unlist(lapply(gtf_gene$V9, function(x) {
  y <- strsplit(x, ';')[[1]][2]
  gsub(' gene_type ','',y)
}))

gene_names <- unlist(lapply(gtf[gtf[,3]=="gene", 9], function(x) {
  y <- strsplit(x, ';')[[1]][1]
  z <- gsub('gene_id ','',y)
  gsub('\\.\\d+','',z)
}))

GWAS_target_trait_ind <- seq(1,nrow(GWAS_clean_4))[GWAS_clean_4$MAPPED_TRAIT %in% target_traits]
all_gene_target_trait <- unique(unlist(all_GWAS_genes_list[GWAS_target_trait_ind]))
all_gene_non_target_trait <- unique(gene_names[!gene_names %in% all_gene_target_trait])

hit_gene_target_trait <- hit_gene_id[hit_gene_id %in% all_gene_target_trait]
hit_gene_non_target_trait <- hit_gene_id[hit_gene_id %in% all_gene_non_target_trait]

hyper_p_val_gene_overall <- phyper(length(hit_gene_target_trait) - 1, # drawn white ball
                                   length(all_gene_target_trait), # total white ball
                                   length(all_gene_non_target_trait), # total black ball
                                   length(hit_gene_id), # num of drawn
                                   lower.tail = FALSE)

circ_master_genes_df <- data.frame(gene = hit_gene_id, 
                                   overlap_with_GWAS_circ_genes = hit_gene_id %in% all_gene_target_trait)
write_xlsx(circ_master_genes_df,"Desktop/circ_gene/data/circ_master_genes_df_v2.xlsx")


hyper_p_val_gene <- vector('numeric', length = length(target_traits))
names(hyper_p_val_gene) <- target_traits
prop_test_by_gene <- vector('list', length = length(target_traits)) 
names(prop_test_by_gene) <- target_traits
enriched_trait_ind_by_gene <- c()
conf_int_diff_by_gene <- c()
for (i in 1:length(target_traits)) {
  GWAS_target_trait_ind <- seq(1,nrow(GWAS_clean_4))[GWAS_clean_4$MAPPED_TRAIT %in% target_traits[i]]
  if (length(GWAS_target_trait_ind)==0) {
    next
  }
  all_gene_target_trait <- unique(unlist(all_GWAS_genes_list[GWAS_target_trait_ind]))
  all_gene_non_target_trait <- unique(gene_names[!gene_names %in% all_gene_target_trait])
  
  hit_gene_target_trait <- hit_gene_id[hit_gene_id %in% all_gene_target_trait]
  hit_gene_non_target_trait <- hit_gene_id[hit_gene_id %in% all_gene_non_target_trait]
  
  hyper_p_val_gene[i] <- phyper(length(hit_gene_target_trait) - 1, # drawn white ball
                                length(all_gene_target_trait), # total white ball
                                length(all_gene_non_target_trait), # total black ball
                                length(hit_gene_id), # num of drawn
                                lower.tail = FALSE)
  
  hit_GWAS_prop <- prop.test(length(hit_gene_target_trait),length(hit_gene_id))
  all_gene_target_trait_clean <- all_gene_target_trait[!all_gene_target_trait %in% hit_gene_target_trait]
  gene_names_clean <- gene_names[!gene_names %in% hit_gene_id]
  all_GWAS_prop <- prop.test(length(all_gene_target_trait_clean),length(gene_names_clean))
  prop_test_by_gene[[i]] <- list(hit_prop=hit_GWAS_prop,all_prop=all_GWAS_prop)
  if (hit_GWAS_prop$conf.int[1] > all_GWAS_prop$conf.int[2]) {
    enriched_trait_ind_by_gene <- c(enriched_trait_ind_by_gene,i)
    conf_int_diff_by_gene <- c(conf_int_diff_by_gene,hit_GWAS_prop$conf.int[1] - all_GWAS_prop$conf.int[2])
  }
}

enriched_traits_by_gene_hyper <- hyper_p_val_gene[hyper_p_val_gene < .05]
enriched_traits_by_gene_proportion <- target_traits[enriched_trait_ind_by_gene]

gene_prop_test_result_df <- c()
for (i in 1:length(prop_test_by_gene)) {
  trait_df <- data.frame(matrix(nrow = 2, ncol = 5))
  colnames(trait_df) <- c('trait', 'percentage', 'lower_lim', 'upper_lim', 'group')
  trait <- target_traits[i]
  trait_df$trait <- trait
  test_result <- prop_test_by_gene[[trait]]
  
  trait_df$group[1] <- 'Circ-SNPs'
  trait_df$percentage[1] <- test_result$hit_prop$estimate
  trait_df$lower_lim[1] <- test_result$hit_prop$conf.int[1]
  trait_df$upper_lim[1] <- test_result$hit_prop$conf.int[2]
  
  trait_df$group[2] <- 'GWAS SNPs'
  trait_df$percentage[2] <- test_result$all_prop$estimate
  trait_df$lower_lim[2] <- test_result$all_prop$conf.int[1]
  trait_df$upper_lim[2] <- test_result$all_prop$conf.int[2]
  
  gene_prop_test_result_df <- rbind(gene_prop_test_result_df,trait_df)
}

discovered_df <- gene_prop_test_result_df[gene_prop_test_result_df$group=='Circ-SNPs',]
trait_ordered <- discovered_df$trait[order(discovered_df$percentage, decreasing = F)]
gene_prop_test_result_df$trait <- factor(gene_prop_test_result_df$trait, levels = trait_ordered)

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "\\1e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "x10<sup>", l)
  # convert exponent
  p1 <- strsplit(l,'-')[[1]][1]
  p2 <- paste0(' ',strsplit(l,'-')[[1]][2])
  p2 <- gsub(' 0', '', p2)
  p2 <- paste0('-',gsub(' ','',p2))
  # return this as a string
  return(paste0(p1,p2,'</sup>'))
}

pValue_mod <- vector('character', length = length(enriched_traits_by_gene_hyper))
for (i in 1:length(pValue_mod)) {
  pValue_mod[i] <- fancy_scientific(signif(enriched_traits_by_gene_hyper[i],2))
}

pValue_y <- vector('numeric', length = length(enriched_traits_by_gene_hyper))
for (i in 1:length(pValue_y)) {
  test_result <- prop_test_by_gene[[names(enriched_traits_by_gene_hyper)[i]]]
  pValue_y[i] <- test_result$hit_prop$conf.int[2] + .022
}

pValue_DF <- data.frame(label = pValue_mod,
                        x = names(enriched_traits_by_gene_hyper),
                        y = pValue_y)

discovered_df_ind <- seq(1,nrow(gene_prop_test_result_df))[gene_prop_test_result_df$group==
                                                             'Circ-SNPs']
non_zero_traits <- as.character(gene_prop_test_result_df$trait[discovered_df_ind][
  gene_prop_test_result_df$percentage[discovered_df_ind] > 0])
gene_prop_test_result_df_clean <- gene_prop_test_result_df[gene_prop_test_result_df$trait %in% non_zero_traits,]

p <- ggplot(gene_prop_test_result_df_clean, aes(x=trait, y=percentage, fill=group)) +
  geom_bar(stat = "identity", position=position_dodge(-.9), color='black') +
  geom_errorbar(aes(ymin=lower_lim, ymax=upper_lim), width=.2, position=position_dodge(-.9)) +
  theme_minimal() +  coord_flip() + 
  theme(legend.position = c(.7,.1), legend.title = element_blank(),
        legend.text = element_text(color='black', size = 8),
        axis.title = element_blank(),
        axis.text.x = element_text(color='black', size = 8),
        axis.text.y = element_text(color='black', size = 8)) +
  geom_richtext(data = pValue_DF, aes(x = x, y = y, label = label), 
                fill = NA, label.color = NA,
                label.padding = grid::unit(rep(0, nrow(pValue_DF)), "pt"),
                size = 3) + ylim(0,.175)
p

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'comparison_by_gene.png',
  plot = last_plot(),
  device = "png",  
  width = 5.25,
  height = 3,
  units = 'in',
  dpi = 300,
  bg = NULL,
)


























