library(ggplot2)

load('Desktop/circ_gene/data/data/R_weighted_complete2.RData')
load('Desktop/circ_gene/data/data/GTEx_colors.RData')
load('Desktop/circ_gene/data/deviation_by_tissue_weighted_complete.RData')

snp_ct_by_tissue <- data.frame(table(R_complete$tissue))
colnames(snp_ct_by_tissue) <- c('tissue','frequency')
snp_ct_by_tissue <- snp_ct_by_tissue[order(snp_ct_by_tissue$frequency,decreasing = T),]
snp_ct_by_tissue_df <- data.frame(tissue = snp_ct_by_tissue$tissue, quantity = snp_ct_by_tissue$frequency, 
                                  data_type = 'number of discovered SNPs')

tissue_name_complete <- names(deviation_by_tissue)
sample_size_complete <- vector('numeric', length = length(tissue_name_complete))
circ_gene_number <- sample_size_complete
for (i in 1:length(sample_size_complete)) {
  deviation <- deviation_by_tissue[[tissue_name_complete[i]]]
  sample_size_complete[i] <- nrow(deviation)
  circ_gene_number[i] <- ncol(deviation)
}

sample_size_df <- data.frame(tissue = tissue_name_complete, quantity = sample_size_complete, 
                             data_type = 'number of individuals')
circ_gene_number_df <- data.frame(tissue = tissue_name_complete, quantity = circ_gene_number, 
                                  data_type = 'number of circadian genes')

missing_tissue <- sample_size_df$tissue[!sample_size_df$tissue %in% snp_ct_by_tissue_df$tissue]

sample_size_df_clean <- sample_size_df[!sample_size_df$tissue %in% missing_tissue,]
circ_gene_number_df_clean <- circ_gene_number_df[!circ_gene_number_df$tissue %in% missing_tissue,]

snp_ct_df <- rbind(snp_ct_by_tissue_df, sample_size_df_clean, circ_gene_number_df_clean)

snp_ct_df$tissue <- factor(snp_ct_df$tissue, levels = snp_ct_by_tissue$tissue)
snp_ct_df$data_type <- factor(snp_ct_df$data_type, levels = c('number of discovered SNPs', 'number of individuals',
                                                    'number of circadian genes'))

color <- vector('character', length = nrow(snp_ct_df))
for (i in 1:length(color)) {
  tissue <- snp_ct_df$tissue[i]
  color[i] <- GTEx_colors$color[GTEx_colors$tissue==tissue]
}

snp_ct_df_colored <- cbind(snp_ct_df,color)

df_1 <- snp_ct_df_colored[snp_ct_df_colored$data_type=='number of discovered SNPs',]
df_2 <- snp_ct_df_colored[snp_ct_df_colored$data_type=='number of individuals',]
df_3 <- snp_ct_df_colored[snp_ct_df_colored$data_type=='number of circadian genes',]

setwd('Desktop/circ_gene/plot/final/')

p1 <- ggplot(data = df_1, aes(x=tissue, y=quantity, fill=tissue)) +
  geom_bar(stat = "identity") + facet_wrap(~data_type, ncol=1) +
  geom_text(aes(label=quantity), vjust=-.1, size=3) + 
  theme(axis.title = element_blank(), legend.position = 'none',
        axis.text.y = element_text(size = 8, color='black'),
        axis.text.x = element_blank(),
        strip.text = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = color) +
  scale_y_continuous(expand=c(0, 0), limits=c(0, 420))
p1
ggsave(
  'snp_number.png',
  plot = last_plot(),
  width = 5.35,
  height = 1.5,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

p2 <- ggplot(data = df_2, aes(x=tissue, y=quantity, fill=tissue)) +
  geom_bar(stat = "identity") + facet_wrap(~data_type, ncol=1) +
  geom_text(aes(label=quantity), vjust=-.1, size=3) + 
  theme(axis.title = element_blank(), legend.position = 'none',
        axis.text.y = element_text(size = 8, color='black'),
        axis.text.x = element_blank(),
        strip.text = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = color) + scale_y_continuous(expand=c(0, 0), limits=c(0, 740))
p2
ggsave(
  'sample_size.png',
  plot = last_plot(),
  width = 5.35,
  height = 1.5,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

p3 <- ggplot(data = df_3, aes(x=tissue, y=quantity, fill=tissue)) +
  geom_bar(stat = "identity") + facet_wrap(~data_type, ncol=1) +
  geom_text(aes(label=quantity), vjust=-.1, size=3) + 
  theme(axis.title = element_blank(), legend.position = 'none',
        axis.text.y = element_text(size = 8, color='black'),
        axis.text.x = element_text(size = 8, color='black',angle = 90, vjust = 0.5, hjust=1),
        strip.text = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_fill_manual(values = color) + scale_y_continuous(expand=c(0, 0), limits=c(0, 3950))
p3
ggsave(
  'gene_number.png',
  plot = last_plot(),
  width = 5.4,
  height = 3,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

