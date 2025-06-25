load('Desktop/circ_gene/data/R_weighted_complete2.RData')

snp_clusters <- vector('list', length = 0)
cluster_names <- c()
window_size <- 500000

chromosomes <- unlist(lapply(R_complete$pos, function(x){
  gsub(':\\d+', '', x)
}))
unique_chromosomes <- sort(unique(chromosomes))
# extract snps by chr
for (j in 1:length(unique_chromosomes)) {
  snps_chr <- R_complete[chromosomes==unique_chromosomes[j],]
  while (nrow(snps_chr)>0) {
    p_val_order <- order(snps_chr$p_val)
    snps_chr <- snps_chr[p_val_order,]
    pos <- as.numeric(unlist(lapply(snps_chr$pos, function(x){
      gsub('\\w+:', '', x)
    })))
    peak_pos <- pos[1]
    cluster_start <- peak_pos-window_size
    if (cluster_start < 0) {
      cluster_start <- 0
    }
    cluster_end <- peak_pos+window_size
    cluster_range <- seq(cluster_start,cluster_end)
    cluster_df <- data.frame(rs_id = snps_chr$RS_id[1], 
                             variant_id = snps_chr$variant_id[1],
                             chr = unique_chromosomes[j], pos = peak_pos, 
                             p_val = snps_chr$p_val[1],
                             tissue = snps_chr$tissue[1], 
                             peak = 'Yes')
    member_ind <- seq(1,length(p_val_order))[pos %in% cluster_range]
    if (length(member_ind)==1) {
      snp_clusters <- c(snp_clusters, list(cluster_df)) 
      cluster_name <- paste0(unique_chromosomes[j],':',cluster_start,'-',cluster_end)
      cluster_names <- c(cluster_names,cluster_name)
      snps_chr <- snps_chr[-member_ind,]
    } else {
      for (l in 2:length(member_ind)) {
        member_df <- data.frame(rs_id = snps_chr$RS_id[member_ind[l]],
                                variant_id = snps_chr$variant_id[member_ind[l]],
                                chr = unique_chromosomes[j], 
                                pos = pos[member_ind[l]], 
                                p_val = snps_chr$p_val[member_ind[l]],
                                tissue = snps_chr$tissue[member_ind[l]], 
                                peak = 'No')
        cluster_df <- rbind(cluster_df,member_df)
      }
      snp_clusters <- c(snp_clusters, list(cluster_df))
      cluster_name <- paste0(unique_chromosomes[j],':',cluster_start,'-',cluster_end)
      cluster_names <- c(cluster_names,cluster_name)
      snps_chr <- snps_chr[-member_ind,]
    }
  }
}

# check if all snps were accounted for
included_variant_id <- c()
for (i in 1:length(snp_clusters)) {
  cluster <- snp_clusters[[i]]
  included_variant_id <- c(included_variant_id, cluster$variant_id)
}
length(included_variant_id)
length(unique(included_variant_id))

cluster_size <- vector('numeric', length = length(snp_clusters))
for (i in 1:length(cluster_size)) {
  cluster_size[i] <- nrow(snp_clusters[[i]])
} 

snp_clusters_clean <- snp_clusters[cluster_size>1]

involved_tissue <- vector('character', length = length(snp_clusters_clean))
all_involved_tissue <- c()
for (i in 1:length(involved_tissue)) {
  tissues <- sort(unique(snp_clusters_clean[[i]]$tissue))
  all_involved_tissue <- union(all_involved_tissue, tissues)
  if (length(tissues)>1) {
    tissues <- paste0(tissues, collapse = '&')
  } 
  involved_tissue[i] <- tissues
}

hit_tissue_tb <- table(involved_tissue)
input <- as.numeric(hit_tissue_tb)
names(input) <- names(hit_tissue_tb)

library(UpSetR)
# Plot
upset_plot <- upset(fromExpression(input), 
                    sets = all_involved_tissue,
                    order.by = "freq", 
                    decreasing = T,
                    mainbar.y.label = "Tissue intersections", 
                    sets.x.label = "Clusters per tissue",
                    mb.ratio = c(.55, 0.45))
                    
png(file='Desktop/circ_gene/plot/final/upset.png',
    width = 5.35, height = 6, units = "in", res = 300) # or other device
upset_plot
dev.off()

# find clusters that hit multiple tissues
hit_tissue_num <- vector('numeric', length = length(snp_clusters_clean))
for (i in 1:length(snp_clusters_clean)) {
  cluster <- snp_clusters_clean[[i]]
  hit_tissue_num[i] <- length(unique(cluster$tissue))
}

clusters_multi_tissue <- snp_clusters_clean[hit_tissue_num>1]
clusters_single_tissue <- snp_clusters_clean[hit_tissue_num==1]

save(list = c('clusters_multi_tissue', 'clusters_single_tissue'),
     file = 'Desktop/circ_gene/data/snp_clusters_weighted_complete_2.RData')




  
  