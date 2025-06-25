segment_ind <-1

# partition work
start_ind <- seq(0,2165,5)
end_ind <- c(seq(4,2165,5),2165)
frag_range <- seq(start_ind[segment_ind],end_ind[segment_ind])

gt_directory <- '~/cmb0/circadian_gene/result/indiv_gt_final'
load('~/cmb0/circadian_gene/data/deviation_by_tissue_weighted_complete.RData')
covariates_path <- '~/cmb0/circadian_gene/data/GTEx_Analysis_v8_eQTL_covariates'
output_dir <- "~/cmb0/circadian_gene/result/regression_weighted_complete2"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

for (l in frag_range) {
  frag_ind <- l
  print(paste0('Frag #',frag_ind))
  gt_file_name <- paste0('indiv_gt_frag_',frag_ind,'.Rdata')
  load(file.path(gt_directory,gt_file_name))
  
  Regression_p_val <- c()
  for (i in 1:length(deviation_by_tissue)) {
    print(paste0('Processing ',i,'/',length(deviation_by_tissue),' tissue.'))
    # load expression deviation
    tissue_name <- names(deviation_by_tissue)[i]
    tissue_data <- deviation_by_tissue[[i]]
    donor_id_deviation <- rownames(tissue_data)
    deviation_sum <- rowSums(tissue_data)
    
    # load covariates
    new_tissue_name <- gsub(' ','_',gsub(' - ','_',gsub(')', '',gsub('\\(','',tissue_name))))
    tissue_covariates_path <- file.path(covariates_path, paste0(new_tissue_name,'.v8.covariates.txt'))
    tissue_covariates <- read.delim(tissue_covariates_path, header = T)
    donor_id_covariates <- unlist(lapply(colnames(tissue_covariates), function(x) {
      gsub('GTEX.','',x)
    }))[-1]
    
    # find donors who have covariates
    donor_id <- intersect(donor_id_deviation, donor_id_covariates)
    
    # assemble covariate and response
    M <- data.frame(matrix(nrow = length(donor_id), ncol = nrow(tissue_covariates) + 1))
    rownames(M) <- donor_id
    colnames(M) <- c(tissue_covariates$ID,'deviation')
    for (j in 1:length(donor_id)) {
      M[j,] <- unlist(c(tissue_covariates[paste0('GTEX.',donor_id[j])],deviation_sum[donor_id[j]]))
    }
    
    # extract gt column
    donor_id_gt <- unlist(lapply(colnames(gt), function(x) {
      gsub('GTEX-','',x)
    }))
    gt_col_ind <- vector('numeric', length = length(donor_id))
    for (j in 1:length(donor_id)) {
      gt_col_ind[j] <- match(donor_id[j],donor_id_gt)
    }
    gt_sub <- gt[,gt_col_ind]
    
    # fit models
    regression_p_val <- data.frame(matrix(nrow = nrow(gt_sub), ncol = 10))
    colnames(regression_p_val) <- c('variant_id','pos','A1','A2','n','ind','coef',
                                    'slope_se','p_val','tissue')
    regression_p_val$tissue <- tissue_name
    bad_snp_ind <- c()
    for (k in 1:nrow(gt_sub)) {
      if (k%%1000==0) {
        print(paste0(round(k/nrow(gt_sub),4)*100,'% of variations processed.'))
      }
      
      gt_snp <- t(gt_sub[k,])
      
      r_tally <- table(gt_snp)
      if (sum(r_tally >= 5)>=2 | length(r_tally)==3) {
        v_id <- variant_id[k]
        regression_p_val$variant_id[k] <- v_id
        regression_p_val$pos[k] <- pos[k]
        v_id_parts <- unlist(strsplit(v_id,'_'))
        regression_p_val$A1[k] <- v_id_parts[3]
        regression_p_val$A2[k] <- v_id_parts[4]
        regression_p_val$n[k] <- nrow(M)
        regression_p_val$ind[k] <- record_ind[k]
        M_gt <- cbind(gt_snp,M)
        colnames(M_gt)[1] <- 'gt'
        model <- lm(deviation ~., data = M_gt)
        model_summary <- summary(model)
        regression_p_val$coef[k] <- model_summary$coefficients[2,1]
        regression_p_val$slope_se[k] <- model_summary$coefficients[2,2]
        regression_p_val$p_val[k] <- model_summary$coefficients[2,4]
      } else {
        bad_snp_ind <- c(bad_snp_ind,k)
        next
      }
    }
    if (length(bad_snp_ind)>0) {
      regression_p_val <- regression_p_val[-bad_snp_ind,]
    }
    Regression_p_val <- rbind(Regression_p_val,regression_p_val)
  }
  # save results
  fileName <- file.path(output_dir,
                        paste0('regression_frag_',frag_ind,'.Rdata'))
  save('Regression_p_val', file = fileName)
  
}

print('Done.')
