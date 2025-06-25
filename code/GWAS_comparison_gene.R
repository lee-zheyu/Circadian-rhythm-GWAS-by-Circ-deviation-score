library(ggplot2)
library(ggtext)

load('Desktop/circ_gene/data/R_weighted_complete2.RData')
load('Desktop/circ_gene/data/gwas_snp_clean_4_complete_v2.RData')

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

gene_names <- unlist(lapply(gtf_gene$V9, function(x) {
  y <- strsplit(x, ';')[[1]][1]
  z <- gsub('gene_id ','',y)
  gsub('\\.\\d+','',z)
}))

gencode_protein_coding_genes <- gene_names[gene_type=='protein_coding']
length(gene_names)

GWAS_target_trait_ind <- seq(1,nrow(GWAS_clean_4))[GWAS_clean_4$MAPPED_TRAIT %in% target_traits]
all_gene_target_trait <- unique(unlist(all_GWAS_genes_list[GWAS_target_trait_ind]))

gencode_target_trait <- all_gene_target_trait[all_gene_target_trait %in% gene_names]
length(gencode_target_trait)
gencode_non_target_triat <- unique(gene_names[!gene_names %in% gencode_target_trait])

hit_gene_id <- c()
for (i in 1:nrow(R_complete)) {
  genes <- unlist(strsplit(R_complete$Gene[i], ';'))
  hit_gene_id <- union(hit_gene_id,genes)
}
hit_gencode_id <- hit_gene_id[hit_gene_id %in% gene_names]
hit_gene_target_trait <- hit_gencode_id[hit_gencode_id %in% gencode_target_trait]
hit_gene_non_target_trait <- hit_gencode_id[hit_gencode_id %in% gencode_non_target_triat]

hyper_p_val_gene_overall <- phyper(length(hit_gene_target_trait) - 1, # drawn white ball
                                   length(gencode_target_trait), # total white ball
                                   length(gencode_non_target_triat), # total black ball
                                   length(hit_gencode_id), # num of drawn
                                   lower.tail = FALSE)

circ_master_genes_df <- data.frame(gene = hit_gencode_id, 
                                   overlap_with_GWAS_circ_genes = hit_gencode_id %in% gencode_target_trait)
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
  gencode_target_trait <- all_gene_target_trait[all_gene_target_trait %in% gene_names]
  gencode_non_target_triat <- unique(gene_names[!gene_names %in% gencode_target_trait])

  hit_gene_target_trait <- hit_gencode_id[hit_gencode_id %in% gencode_target_trait]
 
  hyper_p_val_gene[i] <- phyper(length(hit_gene_target_trait) - 1, # drawn white ball
                                length(gencode_target_trait), # total white ball
                                length(gencode_non_target_triat), # total black ball
                                length(hit_gencode_id), # num of drawn
                                lower.tail = FALSE)
  
  hit_GWAS_prop <- prop.test(length(hit_gene_target_trait),length(hit_gencode_id))
  all_GWAS_prop <- prop.test(length(gencode_target_trait),length(gene_names))
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
  
  trait_df$group[1] <- 'GWAS Circ-SNPs'
  trait_df$percentage[1] <- test_result$hit_prop$estimate
  trait_df$lower_lim[1] <- test_result$hit_prop$conf.int[1]
  trait_df$upper_lim[1] <- test_result$hit_prop$conf.int[2]
  
  trait_df$group[2] <- 'All GWAS SNPs'
  trait_df$percentage[2] <- test_result$all_prop$estimate
  trait_df$lower_lim[2] <- test_result$all_prop$conf.int[1]
  trait_df$upper_lim[2] <- test_result$all_prop$conf.int[2]
  
  gene_prop_test_result_df <- rbind(gene_prop_test_result_df,trait_df)
}

discovered_df <- gene_prop_test_result_df[gene_prop_test_result_df$group=='GWAS Circ-SNPs',]
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
                                                             'GWAS Circ-SNPs']
non_zero_traits <- as.character(gene_prop_test_result_df$trait[discovered_df_ind][
  gene_prop_test_result_df$percentage[discovered_df_ind] > 0])
gene_prop_test_result_df_clean <- gene_prop_test_result_df[gene_prop_test_result_df$trait %in% non_zero_traits,]
gene_prop_test_result_df_clean$group <- factor(gene_prop_test_result_df_clean$group, 
                                               levels = c('GWAS Circ-SNPs', 'All GWAS SNPs'))

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
                size = 3) + ylim(0,.17)
p

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'comparison_by_gene.png',
  plot = last_plot(),
  device = "png",  
  width = 5.5,
  height = 3,
  units = 'in',
  dpi = 300,
  bg = NULL,
)

load('Desktop/circ_gene/data/circadian_genes_v3.RData')

all_circ_genes <- c()
for (i in 1:length(circadian_gene_by_tissue)) {
  tissue_data <- circadian_gene_by_tissue[[i]]
  all_circ_genes <- union(all_circ_genes, tissue_data$gene_id)
}

sum(gencode_target_trait %in% all_circ_genes) / length(gencode_target_trait)

sum(hit_gencode_id %in% all_circ_genes) / length(hit_gencode_id)

