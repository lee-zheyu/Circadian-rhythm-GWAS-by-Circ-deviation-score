tissue_ind <-1

load('~/cmb0/circadian_gene/data/deviation_by_tissue_weighted_complete.RData')

input_dir <- '~/cmb0/circadian_gene/result/regression_weighted_complete2'

output_dir <- '~/cmb0/circadian_gene/result/sumstats_final'
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

tissues <- names(deviation_by_tissue)
for (i in tissue_ind) {
  tissue <- tissues[i]
  print(tissue)
  sumstat_tissue <- data.frame(matrix(nrow = 0, ncol = 6))
  colnames(sumstat_tissue) <- c('Predictor','A1','A2','n','Direction','p_val')
  file_name <- paste0(tissue,'.txt')
  write.table(sumstat_tissue, file = file.path(output_dir,file_name),
              append = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE)
  for (j in 0:2165) {
    load(file.path(input_dir,paste0('regression_frag_',j,'.Rdata')))
    print(paste0('Processing frag #',j,'.'))
    ss_tissue <- Regression_p_val[Regression_p_val$tissue==tissue,]
    ss_tissue_clean <- ss_tissue[!is.na(ss_tissue$pos),]
    
    predictor <- unlist(lapply(ss_tissue_clean$pos, function(x) {
      gsub('chr','',x)
    }))
    
    sumstat_frag <- cbind(predictor,ss_tissue_clean$A1,ss_tissue_clean$A1,ss_tissue_clean$n,
                          ss_tissue_clean$slope_se,ss_tissue_clean$p_val)
    colnames(sumstat_frag) <- c('Predictor','A1','A2','n','Direction','p_val')
    write.table(sumstat_frag, file = file.path(output_dir,file_name),
                append = TRUE, sep = "\t", dec = ".", row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}
