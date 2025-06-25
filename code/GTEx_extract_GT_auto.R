library(vcfR)
setwd('/scratch1/zheyuli/indiv_gt/')

segment_ind <- 1

# partition work
start_ind <- seq(0,348,25)
end_ind <- c(seq(24,348,25),348)
file_ind <- seq(start_ind[segment_ind],end_ind[segment_ind])

outputDir <-  "~/cmb0/circadian_gene/result/indiv_gt2"

set.seed(100)
for (i in file_ind) {
  print(paste0('Processing file #',i))
  
  # load data
  file_name <- paste0(sprintf("p%03d", i),'.vcf')
  indiv_gt <- read.vcfR(file_name)
  fix_region <- getFIX(indiv_gt,getINFO = T)
  
  # check filter
  pass_ind <- seq(1,nrow(fix_region))[fix_region[,7]=='PASS']

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
  record_ind <- paste0(i,':',common_SNP_ind)
  
  # record SNP position
  pos <- paste0(fix_region[common_SNP_ind,1],':',fix_region[common_SNP_ind,2])
  
  # variant id
  variant_id <- fix_region[common_SNP_ind,3]
  
  # extract GT
  gt_region <- indiv_gt@gt[common_SNP_ind,]
  sample_id <- colnames(gt_region)[-1]
  gt <- data.frame(matrix(nrow = nrow(gt_region), ncol = length(sample_id)))
  colnames(gt) <- sample_id
  impute_criteria_ind <- c()
  gt_criteria_ind <- c()
  for (j in 1:nrow(gt)) {
    if (j%%1000==0) {
      print(paste0(j,'/',nrow(gt),' records were processed.'))
    }
    r_og <- gt_region[j,-1]
    r_1 <- vector('numeric',length = length(r_og))
    r_2 <- r_1
    bad_r1_ind <- c()
    bad_r2_ind <- bad_r1_ind
    for (k in 1:length(r_og)) {
      if(is.na(r_og[k])) {
        bad_r1_ind <- c(bad_r1_ind,k)
        bad_r2_ind <- c(bad_r2_ind,k)
        r_1[k] <- NA
        r_2[k] <- NA
      } else {
        n1 <- strsplit(strsplit(r_og[k],':')[[1]][1],'/')[[1]][1]
        if (n1 == '.') {
          r_1[k] <- NA
          bad_r1_ind <- c(bad_r1_ind,k)
        } else {
          r_1[k] <- as.numeric(n1)
        }
        
        n2 <- strsplit(strsplit(r_og[k],':')[[1]][1],'/')[[1]][2]
        if (n2 == '.') {
          r_2[k] <- NA
          bad_r2_ind <- c(bad_r2_ind,k)
        } else {
          r_2[k] <- as.numeric(n2)
        }
      }
    }
    
    if (length(union(bad_r1_ind,bad_r2_ind))>(length(r_og)*.05)) {
      impute_criteria_ind <- c(impute_criteria_ind,j)
      next
    }
    
    r1_p <- sum(r_1==1, na.rm = TRUE)/sum(!is.na(r_1))
    for (k in bad_r1_ind) {
      r_1[k] <- rbinom(1,1,r1_p)
    }
    
    r2_p <- sum(r_2==1, na.rm = TRUE)/sum(!is.na(r_2))
    for (k in bad_r2_ind) {
      r_2[k] <- rbinom(1,1,r2_p)
    }
    
    r_gt <- r_1 + r_2
    
    r_tally <- table(r_gt)
    if (sum(r_tally >= 5)>=2 | length(r_tally)==3) {
      gt[j,] <- r_gt
    } else {
      gt_criteria_ind <- c(gt_criteria_ind,j)
      next
    }
  }
    
  # save results
  obj_list <- c('variant_id','record_ind','pos','gt','impute_criteria_ind','gt_criteria_ind')
  fileName <- file.path(outputDir,paste0('indiv_gt_frag_',i,'.Rdata'))
  save(list = obj_list, file = fileName)
}






