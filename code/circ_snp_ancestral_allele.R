library(writexl)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(biomaRt)

# load('Desktop/circ_gene/data/R_weighted_complete2.RData')
# 
# mart <- useMart(biomart = 'ENSEMBL_MART_SNP', dataset = 'hsapiens_snp', verbose = T)
# 
# snp_id <- R_complete$RS_id
# 
# sublist <- getBM(attributes = c('refsnp_id', 'allele_1', 'ensembl_gene_stable_id'), 
#                  filters = 'snp_filter', values = snp_id, 
#                  mart = mart, verbose = T)
# 
# ancestral_alleles_sub <- sublist[,c(1,2)]
# ancestral_alleles_clean <- ancestral_alleles_sub[ancestral_alleles_sub$allele_1!='',]
# 
# save(ancestral_alleles_clean, 
#      file = 'Desktop/circ_gene/data/circ_snp_ancestral_allele_lookup.RData')
# load('Desktop/circ_gene/data/circ_snp_ancestral_allele_lookup.RData')
# 
# circ_snp_ancestral_allele_df <- R_complete[,c(1,2,4,5)]
# ancestral_allele <- vector('character', length = nrow(circ_snp_ancestral_allele_df))
# missing_ind <- c()
# for (i in 1:length(ancestral_allele)) {
#   rs_id <- circ_snp_ancestral_allele_df$RS_id[i]
#   row_ind <- match(rs_id, ancestral_alleles_clean$refsnp_id)
#   if (!is.na(row_ind)) {
#     ancestral_allele[i] <- ancestral_alleles_clean$allele_1[row_ind]
#   } else {
#     missing_ind <- c(missing_ind, i)
#   }
# }
# circ_snp_ancestral_allele_df2 <- cbind(circ_snp_ancestral_allele_df, ancestral_allele)
# circ_snp_ancestral_allele_df3 <- circ_snp_ancestral_allele_df2[-missing_ind,]
# 
# save(circ_snp_ancestral_allele_df3,
#      file = 'Desktop/circ_gene/data/circ_snp_ancestral_allele_df.RData')
# 
# load('Desktop/circ_gene/data/random_SNP_weighted_complete_2.RData')
# 
# snp_pos <- unlist(lapply(R$pos, function(x){
#   gsub('chr','',x)
# }))
# rs_id <- vector('list', length = nrow(R))
# for (i in 1:length(rs_id)) {
#   print(i)
#   snp_grange <- GRanges(snp_pos[i])
#   overlap <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_grange)
#   if (length(overlap)==0) {
#     next
#   } else {
#     rs_id[[i]] <- overlap@elementMetadata@listData$RefSNP_id
#   }
# }
# save(rs_id, file = 'Desktop/circ_gene/data/random_SNP_rsid.RData')
# load('Desktop/circ_gene/data/random_SNP_rsid.RData')
# 
# rs_id_len <- vector('numeric', length = length(rs_id))
# for (i in 1:length(rs_id_len)) {
#   rs_id_len[i] <- length(rs_id[[i]])
# }
# 
# rs_id_clean <- rs_id[rs_id_len == 1]
# R_clean <- R[rs_id_len == 1,]
# R_2 <- cbind(rs_id = unlist(rs_id_clean), R_clean[,c(1,3,4)])
# 
# snp_id <- R_2$rs_id
# start_ind <- seq(1,length(snp_id),100)
# end_ind <- c(seq(100,length(snp_id),100),length(snp_id))
# 
# mart <- useMart(biomart = 'ENSEMBL_MART_SNP', dataset = 'hsapiens_snp', verbose = T)
# ancestral_alleles_dir <- 'Desktop/circ_gene/data/random_snp_ancestral_allele'
# if (!dir.exists(ancestral_alleles_dir)) {
#   dir.create(ancestral_alleles_dir)
# }
# 
# keep_checking <- 1
# while(keep_checking == 1) {
#   tryCatch({
#     for (i in 1973:length(start_ind)) {
#       print(i)
#       sublist <- getBM(attributes = c('refsnp_id', 'allele_1', 'ensembl_gene_stable_id'), 
#                        filters = 'snp_filter', values = snp_id[start_ind[i]:end_ind[i]], 
#                        mart = mart, verbose = F)
#       
#       ancestral_alleles_sub <- sublist[,c(1,2)]
#       ancestral_alleles_clean <- ancestral_alleles_sub[ancestral_alleles_sub$allele_1!='',]
#       save('ancestral_alleles_clean', 
#            file = file.path(ancestral_alleles_dir, paste0('ancestral_allele_',i,'.RData')))
#     }
#   },
#   error = function(cond) {
#     message("Here's the original error message:")
#     message(conditionMessage(cond))
#     # Choose a return value in case of error
#     NA
#   },
#   warning = function(cond) {
#     message("Here's the original warning message:")
#     message(conditionMessage(cond))
#     # Choose a return value in case of warning
#     NULL
#   },
#   finally = {
#     # NOTE:
#     # Here goes everything that should be executed at the end,
#     # regardless of success or error.
#     # If you want more than one expression to be executed, then you
#     # need to wrap them in curly brackets ({...}); otherwise you could
#     # just have written 'finally = <expression>' 
#     message(paste("Processed chunk:", i))
#   }
#   )
#   keep_checking = 0
# }
# 
# ancestral_alleles_dir <- 'Desktop/circ_gene/data/random_snp_ancestral_allele'
# fileList <- list.files(path = ancestral_alleles_dir)
# filePath <- file.path(ancestral_alleles_dir,fileList)
# 
# ancestral_allele <- c()
# for (i in 1:length(filePath)) {
#   load(filePath[i])
#   ancestral_allele <- rbind(ancestral_allele, ancestral_alleles_clean)
# }
# ancestral_allele_clean <- ancestral_allele[!duplicated(ancestral_allele),]
# 
# save(ancestral_allele_clean, 
#      file = 'Desktop/circ_gene/data/random_ancestral_allele.RData')
# load('Desktop/circ_gene/data/random_ancestral_allele.RData')
# 
# rsid_tally <- table(R_2$rs_id)
# rsid_unique <- names(rsid_tally)[rsid_tally == 1]
# R_2_unique <- R_2[R_2$rs_id %in% rsid_unique,]
# ancestral_allele_random <- vector('character', length = nrow(R_2_unique))
# missing_ind_random <- c()
# for (i in 1:length(ancestral_allele_random)) {
#   rs_id <- R_2_unique$rs_id[i]
#   row_ind <- match(rs_id, ancestral_allele_clean$refsnp_id)
#   if (!is.na(row_ind)) {
#     ancestral_allele_random[i] <- ancestral_allele_clean$allele_1[row_ind]
#   } else {
#     missing_ind_random <- c(missing_ind_random, i)
#   }
# }
# random_snp_ancestral_allele_df <- cbind(R_2_unique, ancestral_allele_random)
# random_snp_ancestral_allele_df_clean <- random_snp_ancestral_allele_df[-missing_ind_random,]
# save(random_snp_ancestral_allele_df_clean,
#      file = 'Desktop/circ_gene/data/random_ancestral_allele_df.RData')
# load('Desktop/circ_gene/data/random_ancestral_allele_df.RData')
# 
# sum(random_snp_ancestral_allele_df_clean$A1 == random_snp_ancestral_allele_df_clean$ancestral_allele) /
#   nrow(random_snp_ancestral_allele_df_clean)
# 
# nrow(random_snp_ancestral_allele_df_clean)

load('Desktop/circ_gene/data/circ_snp_ancestral_allele_df.RData')
load('Desktop/circ_gene/data/random_ancestral_allele_df.RData')

random_snp_ancestral_allele_df2 <- random_snp_ancestral_allele_df_clean[
  !random_snp_ancestral_allele_df_clean$variant_id %in% circ_snp_ancestral_allele_df3$variant_id,]
varaint_ids <- sort(c(circ_snp_ancestral_allele_df3$variant_id,
                      random_snp_ancestral_allele_df2$variant_id))

gt_directory <- 'Desktop/circ_gene/data/indiv_gt2'
gt_files <- list.files(gt_directory)
gt_file_ind <- 1
genotype_group <- data.frame(matrix(nrow = length(varaint_ids), ncol = 5))
colnames(genotype_group) <- c('variant_id', 'frag_name', 'variant_ind', 'major_allele', 'minor_allele')
processed_v_ids <- c()
genotype_group_row_ind <- 1
while (length(processed_v_ids) < length(varaint_ids)) {
  print(gt_file_ind)
  load(file.path(gt_directory,gt_files[gt_file_ind]))
  hit_v_ids <- varaint_ids[varaint_ids %in% variant_id]
  if (length(hit_v_ids) > 0) {
    for (i in 1:length(hit_v_ids)) {
      v_id <- hit_v_ids[i]
      v_id_parts <- unlist(strsplit(v_id, '_'))
      ref_allele <- v_id_parts[3]
      alt_allele <- v_id_parts[4]
      snp_ind <- match(v_id, variant_id)
      gt_snp <- t(gt[snp_ind,])
      ref_allele_freq <- sum(gt_snp == 0)
      alt_allele_freq <- length(gt_snp) - ref_allele_freq
      if (ref_allele_freq > alt_allele_freq) {
        major_allele <- ref_allele
        minor_allele <- alt_allele
      } else if (ref_allele_freq < alt_allele_freq) {
        major_allele <- alt_allele
        minor_allele <- ref_allele
      }
      genotype_group$variant_id[genotype_group_row_ind] <- v_id
      genotype_group$frag_name[genotype_group_row_ind] <- gt_files[gt_file_ind]
      genotype_group$variant_ind[genotype_group_row_ind] <- snp_ind
      genotype_group$major_allele[genotype_group_row_ind] <- major_allele
      genotype_group$minor_allele[genotype_group_row_ind] <- minor_allele
      genotype_group_row_ind <- genotype_group_row_ind + 1
    }
    processed_v_ids <- c(processed_v_ids, hit_v_ids)
  }
  gt_file_ind <- gt_file_ind + 1
} 

save(genotype_group, file = 'Desktop/circ_gene/data/ancestral_allele_lookup.RData')
load('Desktop/circ_gene/data/ancestral_allele_lookup.RData')

mm_allele_circ <- data.frame(matrix(nrow = nrow(circ_snp_ancestral_allele_df3),
                                    ncol = 2))
colnames(mm_allele_circ) <- c('major_allele', 'minor_allele')
for (i in 1:nrow(mm_allele_circ)) {
  v_id <- circ_snp_ancestral_allele_df3$variant_id[i]
  lookup_ind <- match(v_id, genotype_group$variant_id)
  record <- genotype_group[lookup_ind,]
  mm_allele_circ[i,] <- record[,c(4,5)]
  genotype_group[-lookup_ind,]
}
circ_snp_ancestral_allele_df4 <- cbind(circ_snp_ancestral_allele_df3, mm_allele_circ)
sum(circ_snp_ancestral_allele_df4$ancestral_allele == circ_snp_ancestral_allele_df4$major_allele) /
  nrow(circ_snp_ancestral_allele_df4)

write_xlsx(circ_snp_ancestral_allele_df4,
           "Desktop/circ_gene/data/circ_snp_ancestral_allel.xlsx")

mm_allele_random <- data.frame(matrix(nrow = nrow(random_snp_ancestral_allele_df2),
                                      ncol = 2))
colnames(mm_allele_random) <- c('major_allele', 'minor_allele')
for (i in 1:nrow(mm_allele_random)) {
  v_id <- random_snp_ancestral_allele_df2$variant_id[i]
  lookup_ind <- match(v_id, genotype_group$variant_id)
  record <- genotype_group[lookup_ind,]
  mm_allele_random[i,] <- record[,c(4,5)]
  genotype_group <- genotype_group[-lookup_ind,]
}
random_snp_ancestral_allele_df3 <- cbind(random_snp_ancestral_allele_df2, mm_allele_random)
save(random_snp_ancestral_allele_df3,
     file = 'Desktop/circ_gene/data/circ_snp_ancestral_allele_df.RData')
sum(random_snp_ancestral_allele_df3$ancestral_allele == random_snp_ancestral_allele_df3$major_allele) /
  nrow(random_snp_ancestral_allele_df3)

x_freq <- c(sum(circ_snp_ancestral_allele_df4$ancestral_allele == circ_snp_ancestral_allele_df4$major_allele),
            sum(random_snp_ancestral_allele_df3$ancestral_allele == random_snp_ancestral_allele_df3$major_allele))
x_tot <- c(nrow(circ_snp_ancestral_allele_df4), nrow(random_snp_ancestral_allele_df3))
prop.test(x = x_freq , n = x_tot, alternative = 'greater', correct = TRUE)
