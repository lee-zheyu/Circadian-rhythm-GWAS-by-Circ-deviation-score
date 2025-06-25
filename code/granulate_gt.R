gt_directory <- '~/cmb0/circadian_gene/result/indiv_gt3'
new_gt_dir <- '~/cmb0/circadian_gene/result/indiv_gt_final'

frag_range <- seq(0,335)
frag_ind <- 1
head_ind <- seq(1,5000)
new_frag_ind <- 0

GT <- c()
POS <- c()
Record_ind <- c()
Variant_id <- c()
while (frag_ind <= frag_range[length(frag_range)]+1) {
  print(paste0('Frag #',frag_ind))
  gt_file_name <- paste0('indiv_gt_frag_',frag_range[frag_ind],'.Rdata')
  load(file.path(gt_directory,gt_file_name))
  GT <- rbind(GT,gt)
  POS <- c(POS,pos)
  # Record_ind <- c(Record_ind,record_ind)
  Variant_id <- c(Variant_id,variant_id)
  while (nrow(GT)>5000) {
    gt <- GT[head_ind,]
    pos <- POS[head_ind]
    record_ind <- paste0(new_frag_ind,':',seq(1,5000))
    variant_id <- Variant_id[head_ind]
    
    obj_list <- c('variant_id','record_ind','pos','gt')
    fileName <- file.path(new_gt_dir,paste0('indiv_gt_frag_',new_frag_ind,'.Rdata'))
    save(list = obj_list, file = fileName)
    
    GT <- GT[-head_ind,]
    POS <- POS[-head_ind]
    Variant_id <- Variant_id[-head_ind]
    new_frag_ind <- new_frag_ind+1
  }
  frag_ind <- frag_ind+1
}

gt <- GT
pos <- POS
record_ind <- paste0(new_frag_ind,':',seq(1,nrow(gt)))
variant_id <- Variant_id
obj_list <- c('variant_id','record_ind','pos','gt')
fileName <- file.path(new_gt_dir,paste0('indiv_gt_frag_',new_frag_ind,'.Rdata'))
save(list = obj_list, file = fileName)

gt_directory <- '~/cmb0/circadian_gene/result/indiv_gt_sex'
new_gt_dir <- '~/cmb0/circadian_gene/result/indiv_gt_final'

frag_range <- seq(336,348)
frag_ind <- 1
head_ind <- seq(1,5000)
new_frag_ind <- 2095

GT <- c()
POS <- c()
Record_ind <- c()
Variant_id <- c()
while (frag_ind <= length(frag_range)) {
  print(paste0('Frag #',frag_ind))
  gt_file_name <- paste0('indiv_gt_frag_',frag_range[frag_ind],'.Rdata')
  load(file.path(gt_directory,gt_file_name))
  GT <- rbind(GT,gt)
  POS <- c(POS,pos)
  Variant_id <- c(Variant_id,variant_id)
  while (nrow(GT)>5000) {
    gt <- GT[head_ind,]
    pos <- POS[head_ind]
    record_ind <- paste0(new_frag_ind,':',seq(1,5000))
    variant_id <- Variant_id[head_ind]
    
    obj_list <- c('variant_id','record_ind','pos','gt')
    fileName <- file.path(new_gt_dir,paste0('indiv_gt_frag_',new_frag_ind,'.Rdata'))
    save(list = obj_list, file = fileName)
    
    GT <- GT[-head_ind,]
    POS <- POS[-head_ind]
    # Record_ind <- Record_ind[-head_ind]
    Variant_id <- Variant_id[-head_ind]
    new_frag_ind <- new_frag_ind+1
  }
  frag_ind <- frag_ind+1
}

gt <- GT
pos <- POS
record_ind <- paste0(new_frag_ind,':',seq(1,nrow(gt)))
variant_id <- Variant_id
obj_list <- c('variant_id','record_ind','pos','gt')
fileName <- file.path(new_gt_dir,paste0('indiv_gt_frag_',new_frag_ind,'.Rdata'))
save(list = obj_list, file = fileName)
