library(readxl)

# collect circadian genes in baboon by tissue, along with their p values

baboon_gene_db_1 <- read_excel('Desktop/circ_gene/data/Database S1-Baboon_metacycle_statistics.xls',
                              sheet = 1)

baboon_gene_db_2 <- read_excel('Desktop/circ_gene/data/Database S1-Baboon_metacycle_statistics.xls',
                              sheet = 2)

gene_id_baboon <- baboon_gene_db_1$...1[2:nrow(baboon_gene_db_1)]

# load baboon human conversion downloaded from Ensembl
baboon2human_gene_id <- read.delim("Desktop/circ_gene/data/baboon2human_gene_id_v5.txt")
baboon2human_gene_id_1to1 <- baboon2human_gene_id[baboon2human_gene_id$Human.homology.type 
                                                  =='ortholog_one2one',]
baboon2human_gene_DF <- data.frame(baboon = baboon2human_gene_id_1to1$Gene.stable.ID,
                                   human = baboon2human_gene_id_1to1$Human.gene.stable.ID)

gene_id_baboon_clean <- gene_id_baboon[gene_id_baboon %in% baboon2human_gene_DF$baboon]

# record circadian genes by tissue, sheet 1
baboon_gene_db_1_clean <- baboon_gene_db_1[baboon_gene_db_1$...1 %in% gene_id_baboon_clean,
                                           3:ncol(baboon_gene_db_1)]
circadian_gene_by_tissue_1 <- vector('list', length = ncol(baboon_gene_db_1_clean)/6)
tissues_1 <- vector('character',length = length(circadian_gene_by_tissue_1))
p_val_col_ind <- seq(1,ncol(baboon_gene_db_1_clean),6)
for (i in 1:length(p_val_col_ind)) {
  # extract circadian genes based on p-value
  p_val_tissue <- as.numeric(unlist(array(baboon_gene_db_1_clean[,p_val_col_ind[i]])))
  circ_gene_ind <- seq(1,nrow(baboon_gene_db_1_clean))[p_val_tissue<=.05]
  
  # convert baboon gene to human gene
  circ_gene_id <- gene_id_baboon_clean[circ_gene_ind]
  circ_gene_id_new <- vector('character',length = length(circ_gene_ind))
  for (j in 1:length(circ_gene_id)) {
    baboon_gene_id <- circ_gene_id[j]
    circ_gene_id_new[j] <- baboon2human_gene_DF$human[baboon2human_gene_DF$baboon==baboon_gene_id]
  }
  
  # convert p-value
  circ_gene_p_val <- p_val_tissue[circ_gene_ind]
  mu <- mean(circ_gene_p_val)
  std <- sd(circ_gene_p_val)
  circ_gene_z_score <- abs((circ_gene_p_val - mu) / std)
  circ_gene_log_p <- -log10(circ_gene_p_val)
  
  tissues_1[i] <- gsub('...\\d+','',colnames(baboon_gene_db_1_clean)[p_val_col_ind[i]])
  circ_gene_df <- data.frame(gene_id=circ_gene_id_new,z_score=circ_gene_z_score,
                             log_p=circ_gene_log_p)
  circadian_gene_by_tissue_1[[i]] <- circ_gene_df
}
names(circadian_gene_by_tissue_1) <- tissues_1

# record circadian genes by tissue, sheet 2
baboon_gene_db_2_clean <- baboon_gene_db_2[baboon_gene_db_2$...1 %in% gene_id_baboon_clean,
                                           3:ncol(baboon_gene_db_2)]
circadian_gene_by_tissue_2 <- vector('list', length = ncol(baboon_gene_db_2_clean)/6)
tissues_2 <- vector('character',length = length(circadian_gene_by_tissue_2))
for (i in 1:length(p_val_col_ind)) {
  # extract circadian genes based on p-value
  p_val_tissue <- as.numeric(unlist(array(baboon_gene_db_2_clean[,p_val_col_ind[i]])))
  circ_gene_ind <- seq(1,nrow(baboon_gene_db_2_clean))[p_val_tissue<=.05]
  
  # convert baboon gene to human gene
  circ_gene_id <- gene_id_baboon_clean[circ_gene_ind]
  circ_gene_id_new <- vector('character',length = length(circ_gene_ind))
  for (j in 1:length(circ_gene_id)) {
    baboon_gene_id <- circ_gene_id[j]
    circ_gene_id_new[j] <- baboon2human_gene_DF$human[baboon2human_gene_DF$baboon==baboon_gene_id]
  }
  
  # convert p-value
  circ_gene_p_val <- p_val_tissue[circ_gene_ind]
  mu <- mean(circ_gene_p_val)
  std <- sd(circ_gene_p_val)
  circ_gene_z_score <- abs((circ_gene_p_val - mu) / std)
  circ_gene_log_p <- -log10(circ_gene_p_val)
  
  tissues_2[i] <- gsub('...\\d+','',colnames(baboon_gene_db_2_clean)[p_val_col_ind[i]])
  circ_gene_df <- data.frame(gene_id=circ_gene_id_new,z_score=circ_gene_z_score,
                             log_p=circ_gene_log_p)
  circadian_gene_by_tissue_2[[i]] <- circ_gene_df
}
names(circadian_gene_by_tissue_2) <- tissues_2

circadian_gene_by_tissue <- c(circadian_gene_by_tissue_1,circadian_gene_by_tissue_2)

save(circadian_gene_by_tissue, 
     file = 'Desktop/circ_gene/data/circadian_genes_v3.RData')  
  
  







