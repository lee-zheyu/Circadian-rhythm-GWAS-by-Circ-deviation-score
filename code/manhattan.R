library(dplyr)
library(ggplot2)
library(ggrepel)

load('Desktop/circ_gene/data/R_weighted_complete2.RData')
load('Desktop/circ_gene/data/snp_clusters_weighted_complete_2.RData')
load('Desktop/circ_gene/data/GTEx_colors.RData')

string2range <- function(pos, delim=' ', region=TRUE) {
  posp <- as.data.frame(do.call(rbind, strsplit(pos, delim)))
  posp[,1] <- posp[,1]
  posp[,2] <- as.numeric(as.character(posp[,2]))
  if(region) {
    posp[,3] <- as.numeric(as.character(posp[,3]))
  } else {
    posp[,3] <- posp[,2]
  }
  return(posp)
}

snp_df <- data.frame(matrix(nrow = nrow(R_complete),ncol = 5))
colnames(snp_df) <- c('SNP','CHR','BP','P','tissue')
snp_df$SNP <- R_complete$RS_id
snp_df$P <- R_complete$p_val
snpsRanges <- string2range(R_complete$pos, delim=":", region=FALSE)
chr <- unlist(lapply(snpsRanges$V1, function(x) {
  gsub('chr','',x)
}))
chr[chr=='X'] <- '24'
snp_df$CHR <- as.numeric(chr)
snp_df$BP <- snpsRanges$V2
snp_df$tissue <- R_complete$tissue

sig = 5e-8 # significant threshold line
sugg = 1e-6 # suggestive threshold line

threshold <- sig
hlight <- c()
for (i in 1:length(clusters_multi_tissue)) {
  snp_cluster <- clusters_multi_tissue[[i]]
  hlight <- c(hlight,snp_cluster$rs_id[1])
}

df.tmp <- snp_df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 

  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(snp_df, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% hlight, "yes", "no")) %>%
  mutate( is_annotate=ifelse(P < threshold, "yes", "no"))

# get chromosome center positions for x-axis
axisdf <- df.tmp %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf$CHR <- as.character(axisdf$CHR)
axisdf$CHR[axisdf$CHR=='24'] <- 'X'

# get chromosome boundaries
start_coor_bg <- unique(df.tmp$tot)
end_coor_bg <- c((start_coor_bg - 1)[-1],max(df.tmp$BPcum))

chr_cum_info <- snp_df %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len)

cluster_rect_multi <- c()
for (i in 1:length(clusters_multi_tissue)) {
  snp_cluster <- clusters_multi_tissue[[i]]
  
  cluster_df <- data.frame(matrix(nrow = nrow(snp_cluster),ncol = 4))
  colnames(cluster_df) <- c('SNP','CHR','BP','P')
  cluster_df$SNP <- snp_cluster$rs_id
  cluster_df$P <- snp_cluster$p_val
  cluster_df$BP <- snp_cluster$pos
  cluster_df$CHR <- unlist(lapply(snp_cluster$chr, function(x){
    gsub('chr','',x)
  }))
  cluster_df$CHR[cluster_df$CHR=='X'] <- '24'
  cluster_df$CHR <- as.numeric(cluster_df$CHR)
  cluster_df <- left_join(cluster_df, chr_cum_info, by=c("CHR"="CHR")) %>%
    mutate(BPcum=BP+tot) 
  
  x_min <- min(cluster_df$BPcum)-5000000
  x_max <- max(cluster_df$BPcum)+5000000
  y_min <- min(-log10(cluster_df$P))-.01
  y_max <- max(-log10(cluster_df$P))+.01
  cluster_rect_multi <- rbind(cluster_rect_multi,
                        data.frame(x_min = x_min, x_max = x_max, 
                                   y_min = y_min, y_max = y_max))
}

cluster_rect_single <- c()
for (i in 1:length(clusters_single_tissue)) {
  snp_cluster <- clusters_single_tissue[[i]]
  
  if (nrow(snp_cluster)>1) {
    cluster_df <- data.frame(matrix(nrow = nrow(snp_cluster),ncol = 4))
    colnames(cluster_df) <- c('SNP','CHR','BP','P')
    cluster_df$SNP <- snp_cluster$rs_id
    cluster_df$P <- snp_cluster$p_val
    cluster_df$BP <- snp_cluster$pos
    cluster_df$CHR <- unlist(lapply(snp_cluster$chr, function(x){
      gsub('chr','',x)
    }))
    cluster_df$CHR[cluster_df$CHR=='X'] <- '24'
    cluster_df$CHR <- as.numeric(cluster_df$CHR)
    cluster_df <- left_join(cluster_df, chr_cum_info, by=c("CHR"="CHR")) %>%
      mutate(BPcum=BP+tot) 
    
    x_min <- min(cluster_df$BPcum)-5000000
    x_max <- max(cluster_df$BPcum)+5000000
    y_min <- min(-log10(cluster_df$P))-.01
    y_max <- max(-log10(cluster_df$P))+.01
    cluster_rect_single <- rbind(cluster_rect_single,
                                 data.frame(x_min = x_min, x_max = x_max, 
                                            y_min = y_min, y_max = y_max))
  }
}

tissues <- unique(df.tmp$tissue)
df.tmp$tissue <- factor(df.tmp$tissue, levels = tissues)

color <- vector('character', length = length(tissues))
for (i in 1:length(color)) {
  tissue <- tissues[i]
  color[i] <- GTEx_colors$color[GTEx_colors$tissue==tissue]
}

top_snp <- R_complete[-log10(R_complete$p_val)>=10,]
top_snp_genes <-  unique(R_complete$gene_symbols[-log10(R_complete$p_val)>=10])
top_snp_genes_clean <- c()
for (i in 1:length(top_snp_genes)) {
  genes <- unlist(strsplit(top_snp_genes[i], ';'))  
  genes <- genes[genes!='']
  top_snp_genes_clean <- union(top_snp_genes_clean, genes)
}

top_snp_genes_df <- data.frame(gene = top_snp_genes_clean)
write.xlsx(top_snp_genes_df, row.names = FALSE,
           'Desktop/circ_gene/data/top_snp_genes.xlsx')

ggplot(df.tmp, aes(x=BPcum, y=-log10(P), color=tissue)) +
  # Show all points
  geom_point(alpha= .5, size = 1) +
  scale_color_manual(values = color) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(7,14)) +
  
  # add shaded backgroud
  annotate("rect", xmin = start_coor_bg, xmax = end_coor_bg, ymin = 7, ymax = 14,
           alpha = .2,fill = c(rep(c('black','grey'),times=11),'black')) +
  
  # add cluster hilight
  annotate("rect", xmin = cluster_rect_single$x_min, xmax = cluster_rect_single$x_max, 
           ymin = cluster_rect_single$y_min, ymax = cluster_rect_single$y_max,
           alpha = .5,fill = 'blue') +
  annotate("rect", xmin = cluster_rect_multi$x_min, xmax = cluster_rect_multi$x_max, 
           ymin = cluster_rect_multi$y_min, ymax = cluster_rect_multi$y_max,
           alpha = .5,fill = 'red') + xlab('Chromosome') +
  theme(legend.title = element_blank(), legend.position = 'bottom', 
        legend.text = element_text(color = 'black', size=7), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = 'black', size=6),
        axis.title = element_text(color = 'black', size=8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color = guide_legend(nrow = 4, byrow = F)) 

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'sig_snp_manhattan.png',
  plot = last_plot(),
  width = 6,
  height = 3.5,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)  

ggplot(df.tmp, aes(x=BPcum, y=-log10(P), color=tissue)) +
  # Show all points
  geom_point(alpha= .5, size = 1) +
  scale_color_manual(values = color) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center, expand = c(0,0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(7,14)) +
  
  # add shaded backgroud
  annotate("rect", xmin = start_coor_bg, xmax = end_coor_bg, ymin = 7, ymax = 14,
           alpha = .2,fill = c(rep(c('black','grey'),times=11),'black')) +
  
  # add cluster hilight
  annotate("rect", xmin = cluster_rect_single$x_min, xmax = cluster_rect_single$x_max, 
           ymin = cluster_rect_single$y_min, ymax = cluster_rect_single$y_max,
           alpha = .5,fill = 'blue') +
  annotate("rect", xmin = cluster_rect_multi$x_min, xmax = cluster_rect_multi$x_max, 
           ymin = cluster_rect_multi$y_min, ymax = cluster_rect_multi$y_max,
           alpha = .5,fill = 'red') + xlab('Chromosome') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color = 'black', size=5),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  guides(color = guide_legend(nrow = 4, byrow = F)) 

ggsave(
  'sig_snp_manhattan2.png',
  plot = last_plot(),
  width = 6,
  height = 2,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
) 
  
clusters <- c(clusters_multi_tissue, clusters_single_tissue)
peak_snps <- c()
for (i in 1:length(clusters)) {
  cluster <- clusters[[i]]
  peak_snps <- c(peak_snps, cluster$rs_id[cluster$peak=='Yes'])
}

sig_snps <- df.tmp[-log10(df.tmp$P)>=10,]
sig_peak_snps <- sig_snps[sig_snps$SNP %in% peak_snps,]

sig_peak_snps_gene_df <- c()
for (i in 1:length(clusters)) {
  cluster <- clusters[[i]]
  peak_rsid <- cluster$rs_id[cluster$peak=='Yes']
  sig_peak_snps_gene <- c()
  if (is.element(peak_rsid, sig_peak_snps$SNP)) {
    rsid <- cluster$rs_id
    snp_records <- R_complete[R_complete$RS_id %in% rsid,]
    for (j in 1:nrow(snp_records)) {
      snp_record <- snp_records[j,]
      if (snp_record$gene_symbols!='') {
        gene_symbols <- unlist(strsplit(snp_record$gene_symbols, ';'))
        gene_types <- unlist(strsplit(snp_record$gene_types, ';'))
        gene_symbols_clean <- gene_symbols[gene_types=="protein_coding"]
        sig_peak_snps_gene <- c(sig_peak_snps_gene, gene_symbols_clean[gene_symbols_clean!='']) 
      }
    }
  }
  genes_unique <- unique(sig_peak_snps_gene)
  genes_clean <- paste(genes_unique[!is.na(genes_unique)], collapse = ';')
  if (length(sig_peak_snps_gene)>0) {
    snp_df <- df.tmp[df.tmp$SNP==peak_rsid,]
    peak_record <- data.frame(rsid = peak_rsid, chr = snp_df$CHR, BP = snp_df$BP, y = -log10(snp_df$P),
                              gene_symbol = genes_clean)
    sig_peak_snps_gene_df <- rbind(sig_peak_snps_gene_df, peak_record)
  }
}

sig_peak_snps_gene_df_ordered <- sig_peak_snps_gene_df[order(sig_peak_snps_gene_df$chr),]







