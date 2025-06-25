library(ggplot2)
library(moments)

load('Desktop/circ_gene/data/deviation_by_tissue_weighted_complete.RData')
load('Desktop/circ_gene/data/GTEx_colors.RData')

deviation_score <- c()
deviation_score_kurtosis <- vector('numeric', length = length(deviation_by_tissue))
names(deviation_score_kurtosis) <- names(deviation_by_tissue)
deviation_score_95th <- deviation_score_kurtosis
for (i in 1:length(deviation_by_tissue)) {
  tissue_name <- names(deviation_by_tissue)[i]
  tissue_data <- deviation_by_tissue[[i]]
  tissue_color <- GTEx_colors$color[GTEx_colors$tissue==tissue_name]
  deviation_sum <- rowSums(tissue_data)
  deviation_score_kurtosis[i] <- kurtosis(deviation_sum)
  deviation_score_tissue <- data.frame(deviation_sum=deviation_sum,
                                       tissue_name=tissue_name,
                                       tissue_color=tissue_color)
  deviation_score <- rbind(deviation_score, deviation_score_tissue)
  
  deviation_score_95th[i] <- quantile(deviation_score_tissue$deviation_sum, probs = .95)
}

ggplot(data = deviation_score, aes(x=deviation_sum, color=tissue_name)) + 
  stat_density(size = 1 ,geom="line",position="identity") +
  scale_color_manual(values = GTEx_colors$color) +
  facet_wrap(~tissue_name, ncol = 4, strip.position = 'bottom', scales = "free") +
  theme(legend.position = 'none',axis.title = element_text(color = 'black', size = 8),
        axis.text = element_text(color = 'black', size = 6),strip.text = element_text(color = 'black',size=6)) +
  xlab('Circadian deviation score') + ylab('Density')

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'dev_score_density.png',
  plot = last_plot(),
  device = "png",
  width = 7,
  height = 7,
  units = 'in',
  dpi = "print",
  bg = NULL,
)


