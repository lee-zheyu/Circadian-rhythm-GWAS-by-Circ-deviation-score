library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)

load('Desktop/circ_gene/data/R_weighted_complete2.RData')

# load annotation database
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
Exons <- exonsBy(txdb, by="tx", use.names=TRUE)
CDS <- cdsBy(txdb, by='tx', use.names=TRUE)
UTR5 <- fiveUTRsByTranscript(txdb, use.names=TRUE)
UTR3 <- threeUTRsByTranscript(txdb, use.names=TRUE)
Introns <- intronsByTranscript(txdb, use.names=TRUE)

# load annotation table
gtf_file <- 'Desktop/circ_gene/data/annotation_table.RData'
load(gtf_file)

snp_location <- vector('list', length = nrow(R_complete))
for (i in 1:nrow(R_complete)) {
  print(i)
  
  # construct range info
  x <- GRanges(R_complete$pos[i])
  
  # map SNP
  cds_coor <- mapToTranscripts(x, CDS)
  utr5_coor <- mapToTranscripts(x,UTR5)
  utr3_coor <- mapToTranscripts(x,UTR3)
  intron_coor <- mapToTranscripts(x, Introns)
  
  # record location
  loc <- c()
  if (length(cds_coor)>0) {
    loc <- c(loc,'CDS')
  }
  if (length(utr5_coor)>0) {
    loc <- c(loc,"5' UTR")
  }
  if (length(utr3_coor)>0) {
    loc <- c(loc,"3' UTR")
  }
  if (length(intron_coor)>0) {
    loc <- c(loc,'Intron')
  }
  if (length(loc)==0) {
    loc <- 'IGR'
  }
  
  snp_location[[i]] <- loc
}

save('snp_location', file = 'Desktop/circ_gene/data/snp_location_weighted_complete2.RData')

load('Desktop/circ_gene/data/snp_location_weighted_complete2.RData')

snp_loc_table <- data.frame(sort(table(unlist(snp_location)), decreasing = T))
colnames(snp_loc_table) <- c('region','frequency')

regions <- as.character(snp_loc_table$region)
regions_ind <- vector('numeric', length = nrow(snp_loc_table))
for (i in 1:length(regions_ind)) {
  regions_ind[i] <- match(snp_loc_table$region[i], regions)
}

snp_loc_table <- snp_loc_table[regions_ind,]


perc <- round(snp_loc_table$frequency/sum(snp_loc_table$frequency), digits = 4)
snp_loc_table <- cbind(snp_loc_table,perc)

snp_loc_table2 <- snp_loc_table %>% 
  mutate(csum = rev(cumsum(rev(perc))), 
         pos = perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), perc/2, pos))

ggplot(snp_loc_table2, aes(x = "" , y = perc, fill = fct_inorder(region))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set2") +
  geom_label_repel(data = snp_loc_table2,
                   aes(y = pos, label = paste0(perc*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Region")) +
  theme_void() + theme(legend.position = 'none') 
  

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'snp_loc.png',
  plot = last_plot(),
  width = 2.2,
  height = 2.2,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

snp_loc_dir <- 'Desktop/circ_gene/data/snp_loc'
snp_loc_files <- list.files(snp_loc_dir)

SNP_loc <- c()
for (i in 1:length(snp_loc_files)) {
  load(file.path(snp_loc_dir, snp_loc_files[i]))
  snp_loc_clean <- snp_loc[snp_loc!='']
  SNP_loc <- c(SNP_loc, snp_loc_clean)
}

snp_loc_table_all <- data.frame(sort(table(SNP_loc), decreasing = T))
colnames(snp_loc_table_all) <- c('region','frequency')

regions_ind <- vector('numeric', length = nrow(snp_loc_table_all))
for (i in 1:length(regions_ind)) {
  regions_ind[i] <- match(snp_loc_table_all$region[i], regions)
}

snp_loc_table_all <- snp_loc_table_all[regions_ind,]

perc <- round(snp_loc_table_all$frequency/sum(snp_loc_table_all$frequency), digits = 4)
snp_loc_table_all <- cbind(snp_loc_table_all,perc)

snp_loc_table_all_2 <- snp_loc_table_all %>% 
  mutate(csum = rev(cumsum(rev(perc))), 
         pos = perc/2 + lead(csum, 1),
         pos = if_else(is.na(pos), perc/2, pos))


ggplot(snp_loc_table_all_2, aes(x = "" , y = perc, fill = fct_inorder(region))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set2") +
  geom_label_repel(data = snp_loc_table_all_2,
                   aes(y = pos, label = paste0(perc*100, "%")),
                   size = 4.5, nudge_x = 1, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Region")) +
  theme_void() + theme(legend.title = element_blank())

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'snp_loc_all.png',
  plot = last_plot(),
  width = 3,
  height = 3,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

p_val_df <- c()
hit_tot <- sum(snp_loc_table$frequency)
all_tot <- sum(snp_loc_table_all$frequency)
for (i in 1:length(regions)) {
  region <- regions[i]
  
  hit_freq <- snp_loc_table$frequency[snp_loc_table$region == region]
  hit_perc <- hit_freq / hit_tot
  all_freq <- snp_loc_table_all$frequency[snp_loc_table_all$region == region]
  all_perc <- all_freq / all_tot
  
  hit_binom_test <- binom.test(hit_freq, hit_tot, p=all_perc, alternative = 'greater')
  
  region_df <- data.frame(region = region, 
                          hit_freq = hit_freq, hit_perc = hit_perc,
                          all_freq = all_freq, all_perc = all_perc,
                          p_val = hit_binom_test$p.value)
  
  p_val_df <- rbind(p_val_df, region_df)
}





