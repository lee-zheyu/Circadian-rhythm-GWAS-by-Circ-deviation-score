setwd('~/cmb0/circadian_gene/result/regression_weighted_complete2')

file_names <- list.files()
R <- c()
record_num_significant_F <- vector('numeric', length = length(file_names))
record_num_all_F <- record_num_significant_F
for (i in 1:length(file_names)) {
  print(i)
  load(file_names[i])
  Regression_p_val_clean <- Regression_p_val[!is.na(Regression_p_val$p_val),]
  significant_record <- Regression_p_val[Regression_p_val$p_val<5E-8,]
  R <- rbind(R,significant_record)
  record_num_significant_F[i] <- nrow(significant_record)
  record_num_all_F[i] <- nrow(Regression_p_val)
}

save('R',
     file = '~/cmb0/circadian_gene/result/significant_SNP_weighted_complete2.RData')

save(list = c('record_num_significant_F', 'record_num_all_F'),
     file = '~/cmb0/circadian_gene/result/snp_num.RData')





