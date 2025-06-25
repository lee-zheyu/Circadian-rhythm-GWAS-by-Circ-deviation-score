library(CePa)
library(edgeR)
library(RNOmni)

# load GTEx gene count and sample attributes
gene_read_ct <- read.gct("~/cmb0/circadian_gene/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz")
gene_read_tpm <- read.gct("~/cmb0/circadian_gene/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
sample_attributes <- read.delim("~/cmb0/circadian_gene/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
                                header = TRUE)

# load tissue category mapping from baboon to GTEx
category_mapping <- read.delim("~/cmb0/circadian_gene/data/GTEx_tissues.txt",
                                  header = TRUE)
category_mapping_clean <- category_mapping[category_mapping$baboon!='',]

# load circadian gene ID
load('~/cmb0/circadian_gene/data/circadian_genes_v3.RData')
print('Finished loading files.')

# modify formats of sample ID and gene ID in read counts data
read_ct_sample_id <- colnames(gene_read_ct)
for (k in 1:length(read_ct_sample_id)) {
  read_ct_sample_id[k] <- paste0(unlist(strsplit(read_ct_sample_id[k],'[.]')), collapse = '-')
}
colnames(gene_read_ct) <- read_ct_sample_id

read_ct_gene_id <- rownames(gene_read_ct)
for (k in 1:length(read_ct_gene_id)) {
  read_ct_gene_id[k] <- unlist(strsplit(read_ct_gene_id[k],'[.]'))[1]
}
rownames(gene_read_ct) <- read_ct_gene_id

read_tpm_sample_id <- colnames(gene_read_tpm)
for (k in 1:length(read_tpm_sample_id)) {
  read_tpm_sample_id[k] <- paste0(unlist(strsplit(read_tpm_sample_id[k],'[.]')), collapse = '-')
}
colnames(gene_read_tpm) <- read_tpm_sample_id

read_tpm_gene_id <- rownames(gene_read_tpm)
for (k in 1:length(read_tpm_gene_id)) {
  read_tpm_gene_id[k] <- unlist(strsplit(read_tpm_gene_id[k],'[.]'))[1]
}
rownames(gene_read_tpm) <- read_tpm_gene_id

# extract read counts of circadian genes by tissue
deviation_by_tissue <- vector('list',length = nrow(category_mapping_clean))
tissue_names <- category_mapping_clean$GTEx
names(deviation_by_tissue) <- tissue_names
for (i in 1:length(deviation_by_tissue)) { 
  # extract expression values per tissue
  GTEx_tissue <- tissue_names[i]
  print(GTEx_tissue)
  GTEx_sample_id <- sample_attributes$SAMPID[sample_attributes$SMTSD==GTEx_tissue]
  tissue_ct <- t(gene_read_ct[,read_ct_sample_id%in%GTEx_sample_id])
  tissue_tpm <- t(gene_read_tpm[,read_tpm_sample_id%in%GTEx_sample_id])
  
  # find cols with not enough ct
  col_ind_valid_ct <- c()
  for (j in 1:ncol(tissue_ct)) {
    valid_ct_num <- sum(tissue_ct[,j]>=6)
    if (valid_ct_num > nrow(tissue_ct)*.2) {
      col_ind_valid_ct <- c(col_ind_valid_ct,j)
    }
  }
  
  # find cols with not enough tpm
  col_ind_valid_tpm <- c()
  for (j in 1:ncol(tissue_tpm)) {
    valid_ct_num <- sum(tissue_tpm[,j]>.1)
    if (valid_ct_num > nrow(tissue_tpm)*.2) {
      col_ind_valid_tpm <- c(col_ind_valid_tpm,j)
    }
  }
  
  # remove genes w/o enough data
  col_ind_valid <- union(col_ind_valid_ct,col_ind_valid_tpm)
  tissue_tpm_clean <- tissue_tpm[,col_ind_valid]
  
  # normalization
  # TMM
  d <- DGEList(counts=tissue_tpm_clean, group=colnames(tissue_tpm_clean))
  TMM <- calcNormFactors(d, method="TMM")
  tissue_tpm_TMM <- cpm(TMM)
  
  # collect circadian gene ID
  circ_gene_by_tissue <- c()
  baboon_tissue <- unlist(strsplit(category_mapping_clean$baboon[i],';'))
  for (j in 1:length(baboon_tissue)) {
    circ_gene_by_tissue <- rbind(circ_gene_by_tissue,circadian_gene_by_tissue[[baboon_tissue[j]]])
  }

  if (is.null(circ_gene_by_tissue)) {
    next
  } else {
    redundant_gene_id <- circ_gene_by_tissue$gene_id[duplicated(circ_gene_by_tissue$gene_id)]
    
    if (length(redundant_gene_id)>0) {
      Redundant_records <- c()
      for (k in 1:length(redundant_gene_id)) {
        redundant_records <- circ_gene_by_tissue[circ_gene_by_tissue$gene_id==redundant_gene_id[k],]
        redundant_rows <- order(redundant_records$log_p, decreasing = TRUE)[-1]
        Redundant_records <- rbind(Redundant_records,redundant_records[redundant_rows,])
      }
      
      redundant_rowname <- rownames(Redundant_records)
      redundant_row_ind <- seq(1,nrow(circ_gene_by_tissue))[rownames(circ_gene_by_tissue) %in% redundant_rowname]
      circ_gene_by_tissue_clean <- circ_gene_by_tissue[-redundant_row_ind,]
    } else {
      circ_gene_by_tissue_clean <- circ_gene_by_tissue
    }
  }
  
  gene_id_by_tissue <- circ_gene_by_tissue_clean$gene_id
  
  # extract read counts
  circadian_gene_tpm <- tissue_tpm_TMM[,colnames(tissue_tpm_TMM)%in%gene_id_by_tissue]
  # change column name to donor ID
  donor_id <- rownames(circadian_gene_tpm)
  for (k in 1:length(donor_id)) {
    donor_id[k] <- unlist(strsplit(donor_id[k],'-'))[2] 
  }
  rownames(circadian_gene_tpm) <- donor_id
   
  # INT & deviation
  circadian_gene_INT <- c()
  circ_gene_id <- colnames(circadian_gene_tpm)
  for (j in 1:ncol(circadian_gene_tpm)) {
    new_col <- RankNorm(circadian_gene_tpm[,j], k = 0.375)
    weight <- circ_gene_by_tissue_clean$log_p[circ_gene_by_tissue_clean$gene_id==circ_gene_id[j]]
    deviation_col <- abs(new_col - median(new_col)) * weight
    circadian_gene_INT <- cbind(circadian_gene_INT,deviation_col)
  }
  colnames(circadian_gene_INT) <- colnames(circadian_gene_tpm)
  
  deviation_by_tissue[[i]] <- circadian_gene_INT
}

print('Finished compiling deviations.')

valid_ind <- c()
for (i in 1:length(deviation_by_tissue)) {
  if (!is.null(deviation_by_tissue[[i]]) ) {
    if (nrow(deviation_by_tissue[[i]])>=70) {
      valid_ind <- c(valid_ind,i)
    }
  }
}
deviation_by_tissue <- deviation_by_tissue[valid_ind]
print('Finished cleaning deviation list.')

save('deviation_by_tissue',
     file = '~/cmb0/circadian_gene/data/deviation_by_tissue_weighted_complete.RData')
print('Data saved.')




