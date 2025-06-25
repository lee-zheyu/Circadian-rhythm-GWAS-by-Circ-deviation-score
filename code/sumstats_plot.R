library(qqman)
library(ggplot2)

tissue_sumstats_dir <- '~/cmb0/circadian_gene/result/sumstats_final'
tissue_sumstats_files <- list.files(tissue_sumstats_dir)

for (i in 1:length(tissue_sumstats_files)) {
  tissue_file <- tissue_sumstats_files[i]
  print(tissue_file)
  tissue <- gsub('\\.txt','',tissue_file)
  tissue_sumstats <- read.delim(file.path(tissue_sumstats_dir, tissue_file))
  
  CHR <- vector('numeric', length = nrow(tissue_sumstats))
  BP <- CHR
  
  for (j in 1:nrow(tissue_sumstats)) {
    CHR_BP <- unlist(strsplit(tissue_sumstats$Predictor[j], ':'))
    CHR[j] <- CHR_BP[1]
    BP[j] <- CHR_BP[2]
  }
  
  CHR[CHR == 'X'] <- '23'
  CHR <- as.numeric(CHR)
  BP <- as.numeric(BP)
  
  gwasResults <- data.frame(SNP = paste0('SNP', seq(1, nrow(tissue_sumstats))),
                            CHR = CHR, BP = BP, P = tissue_sumstats$p_val)
  
  plots_dir <- '~/cmb0/circadian_gene/result/sumstats_final_plots'
  
  print('Ploting Manhattan')
  man_plot_path <- file.path(plots_dir, paste0(tissue, '_manhattan.png'))
  png(man_plot_path, width = 6, height = 4, units = 'in', res = 300)
  manhattan(gwasResults, main = "Manhattan Plot", cex = 0.6, cex.axis = 0.9,
            chrlabs = c(1:22,'X'))
  dev.off()
  print('Manhattan saved.')
  
  print('Ploting Q-Q')
  qq_plot_path <- file.path(plots_dir, paste0(tissue, '_qq.png'))
  png(qq_plot_path, width = 5, height = 4, units = 'in', res = 300)
  qq(gwasResults$P)
  dev.off()
  print('Q-Q saved.')
}

