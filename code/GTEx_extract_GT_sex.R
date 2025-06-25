library(vcfR)

load('~/cmb0/circadian_gene/result/indiv_sex.RData')

donor_id_male <- sex_df_clean$sample_id[sex_df_clean$sex==1] 
donor_id_female <- sex_df_clean$sample_id[sex_df_clean$sex==2] 

setwd('/scratch1/zheyuli/indiv_gt/')

segment_ind <- 1

# partition work
start_ind <- seq(0,348,25)
end_ind <- c(seq(24,348,25),348)
file_ind <- seq(start_ind[segment_ind],end_ind[segment_ind])

outputDir <-  "~/cmb0/circadian_gene/result/indiv_gt/chrX"

set.seed(100)
for (i in file_ind) {
  print(paste0('Processing file #',i))
  
  # load data
  file_name <- paste0(sprintf("p%03d", i),'.vcf')
  indiv_gt <- read.vcfR(file_name)
  fix_region <- getFIX(indiv_gt,getINFO = T)
  
  # check filter & chr
  pass_ind <- seq(1,nrow(fix_region))[fix_region[,7]=='PASS' & 
                                        fix_region[,1]=='chrX']
  
  if (length(pass_ind)==0) {
    next
  }

  # extract SNP record
  SNP_ind <- c()
  for (j in 1:length(pass_ind)) {
    ref_len <- nchar(fix_region[pass_ind[j],4])
    alt_len <- nchar(fix_region[pass_ind[j],5])
    
    if (ref_len==1 & alt_len==1) {
      SNP_ind <- c(SNP_ind,pass_ind[j])
    }
  }
  
  # extract common SNPs
  common_SNP_ind <- c()
  for (j in 1:length(SNP_ind)) {
    info <- fix_region[SNP_ind[j],8]
    allele_freq <- as.numeric(strsplit(strsplit(info,';')[[1]][2],'=')[[1]][2])
    if (allele_freq >= .01) {
      common_SNP_ind <- c(common_SNP_ind,SNP_ind[j])
    }
  }
  # Record_ind <- c(Record_ind,paste0(i,':',common_SNP_ind))
  record_ind <- paste0(i,':',common_SNP_ind)
  
  # record SNP position
  pos <- paste0(fix_region[common_SNP_ind,1],':',fix_region[common_SNP_ind,2])
  
  # variant id
  variant_id <- fix_region[common_SNP_ind,3]
  
  # extract GT
  gt_region <- indiv_gt@gt[common_SNP_ind,]
  sample_id <- colnames(gt_region)[-1]
  gt <- data.frame(matrix(nrow = nrow(gt_region), ncol = nrow(sex_df_clean)))
  colnames(gt) <- c(colnames(gt_region)[colnames(gt_region) %in% donor_id_male],
                    colnames(gt_region)[colnames(gt_region) %in% donor_id_female])
  impute_criteria_ind <- c()
  gt_criteria_ind <- c()
  for (j in 1:nrow(gt)) {
    if (j%%1000==0) {
      print(paste0(j,'/',nrow(gt),' records were processed.'))
    }
    
    # male
    r_og_m <- gt_region[j,colnames(gt_region) %in% donor_id_male]
    r_1_m <- vector('numeric',length = length(r_og_m))
    r_2_m <- r_1_m
    bad_r1_ind_m <- c()
    bad_r2_ind_m <- bad_r1_ind_m
    for (k in 1:length(r_og_m)) {
      if(is.na(r_og_m[k])) {
        r_1_m[k] <- NA
        bad_r1_ind_m <- c(bad_r1_ind_m,k)
        r_2_m[k] <- NA
        bad_r2_ind_m <- c(bad_r2_ind_m,k)
        next
      }
      
      n1 <- strsplit(strsplit(r_og_m[k],':')[[1]][1],'/')[[1]][1]
      if (n1 == '.') {
        r_1_m[k] <- NA
        bad_r1_ind_m <- c(bad_r1_ind_m,k)
      } else {
        r_1_m[k] <- as.numeric(n1)
      }
      
      n2 <- strsplit(strsplit(r_og_m[k],':')[[1]][1],'/')[[1]][2]
      if (n2 == '.') {
        r_2_m[k] <- NA
        bad_r2_ind_m <- c(bad_r2_ind_m,k)
      } else {
        r_2_m[k] <- as.numeric(n2)
      }
    }
    
    if (length(union(bad_r1_ind_m,bad_r2_ind_m))>(length(r_og_m)*.05)) {
      impute_criteria_ind <- c(impute_criteria_ind,j)
      next
    }
    
    gt_m <- r_1_m + r_2_m
    
    r_p_m <- sum(gt_m==2, na.rm = TRUE)/sum(!is.na(gt_m))
    
    for (l in union(bad_r1_ind_m,bad_r2_ind_m)) {
      gt_m[l] <- rbinom(1,1,r_p_m)*2
    }
    
    bad_gt_m_ind <- seq(1,length(gt_m))[gt_m==1]
    for (l in bad_gt_m_ind) {
      gt_m[l] <- rbinom(1,1,.5)*2
    }
    
    r_tally_m <- table(gt_m)
    if (sum(r_tally_m >= 5)>=2) {
      m_pass <- TRUE
    } else {
      gt_criteria_ind <- c(gt_criteria_ind,j)
      next
    }
    
    # female
    r_og_f <- gt_region[j,colnames(gt_region) %in% donor_id_female]
    r_1_f <- vector('numeric',length = length(r_og_f))
    r_2_f <- r_1_f
    bad_r1_ind_f <- c()
    bad_r2_ind_f <- bad_r1_ind_f
    for (k in 1:length(r_og_f)) {
      if(is.na(r_og_f[k])) {
        r_1_f[k] <- NA
        bad_r1_ind_f <- c(bad_r1_ind_f,k)
        r_2_f[k] <- NA
        bad_r2_ind_f <- c(bad_r2_ind_f,k)
        next
      }
      
      n1 <- strsplit(strsplit(r_og_f[k],':')[[1]][1],'/')[[1]][1]
      if (n1 == '.') {
        r_1_f[k] <- NA
        bad_r1_ind_f <- c(bad_r1_ind_f,k)
      } else {
        r_1_f[k] <- as.numeric(n1)
      }
      
      n2 <- strsplit(strsplit(r_og_f[k],':')[[1]][1],'/')[[1]][2]
      if (n2 == '.') {
        r_2_f[k] <- NA
        bad_r2_ind_f <- c(bad_r2_ind_f,k)
      } else {
        r_2_f[k] <- as.numeric(n2)
      }
    }
    
    if (length(union(bad_r1_ind_f,bad_r2_ind_f))>(length(r_og_f)*.05)) {
      impute_criteria_ind <- c(impute_criteria_ind,j)
      next
    }
    
    r1_p_f <- sum(r_1_f==1, na.rm = TRUE)/sum(!is.na(r_1_f))
    r2_p_f <- sum(r_2_f==1, na.rm = TRUE)/sum(!is.na(r_2_f))
    for (l in union(bad_r1_ind_f,bad_r2_ind_f)) {
      r_1_f[l] <- rbinom(1,1,r1_p_f)
      r_2_f[l] <- rbinom(1,1,r2_p_f)
    }
    gt_f <- r_1_f + r_2_f
    
    r_tally_f <- table(gt_f)
    if (sum(r_tally_f >= 5)>=2 | length(r_tally)==3) {
      female_pass <- TRUE
    } else {
      gt_criteria_ind <- c(gt_criteria_ind,j)
      next
    }
    
    gt[j,] <- c(gt_m,gt_f)
  }
    
  # save results
  obj_list <- c('variant_id','record_ind','pos','gt','impute_criteria_ind','gt_criteria_ind')
  fileName <- file.path(outputDir,paste0('chrX_gt_',i,'.Rdata'))
  save(list = obj_list, file = fileName)
}






