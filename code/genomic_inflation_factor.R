tissue_sumstats_dir <- '~/cmb0/circadian_gene/result/sumstats_final'
tissue_sumstats_files <- list.files(tissue_sumstats_dir)

genomic_inflation_factor <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(genomic_inflation_factor) <- c('tissue','lambda')
write.table(genomic_inflation_factor, file = '~/cmb0/circadian_gene/result/genomic_inflation_factor.txt',
            append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)

for (i in 1:length(tissue_sumstats_files)) {
  tissue_file <- tissue_sumstats_files[i]
  print(tissue_file)
  tissue <- gsub('\\.txt','',tissue_file)
  tissue_sumstats <- read.delim(file.path(tissue_sumstats_dir, tissue_file))
  
  pvals <- tissue_sumstats$p_val
  
  chisq_stats <- qchisq(1 - pvals, df = 1)
  
  median_chisq_obs <- median(chisq_stats, na.rm = TRUE)
  
  lambda_gc <- median_chisq_obs / .456
  
  lambda_tissue <- data.frame(tissue = tissue, lambda = lambda_gc)
  write.table(lambda_tissue, file = '~/cmb0/circadian_gene/result/genomic_inflation_factor.txt',
              append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}




