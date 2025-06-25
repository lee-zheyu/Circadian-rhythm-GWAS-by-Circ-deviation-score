gt_directory <- '~/cmb0/circadian_gene/result/indiv_gt/chrX'
output_dir <- "~/cmb0/circadian_gene/result/indiv_gt_sex"

frag_range <- seq(336,348)

for (l in frag_range) {
  frag_ind <- l
  print(paste0('Frag #',frag_ind))
  gt_file_name <- paste0('chrX_gt_',frag_ind,'.Rdata')
  load(file.path(gt_directory,gt_file_name))
  bad_ind <- union(gt_criteria_ind,impute_criteria_ind)
  gt <- gt[-bad_ind,]
  variant_id <- variant_id[-bad_ind]
  record_ind <- record_ind[-bad_ind]
  pos <- pos[-bad_ind]
  
  # save results
  obj_list <- c('variant_id','record_ind','pos','gt')
  fileName <- file.path(output_dir,paste0('indiv_gt_frag_',l,'.Rdata'))
  save(list = obj_list, file = fileName)
}
