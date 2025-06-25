library(GenomicFeatures)
library(ggplot2)

load('Desktop/circ_gene/data/R_weighted_complete2.RData')

## import GTF
gtf_file_path <- 'Desktop/circ_gene/data/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz'
gtf <- makeTxDbFromGFF(gtf_file_path)

##drop non-standard chromosomes
gtf <- keepStandardChromosomes(gtf,pruning='coarse')
all_transcripts <- transcripts(gtf, use.names = TRUE)

snp_location <- vector('list', length = nrow(R_complete))
names(snp_location) <- R_complete$variant_id
non_tx_ind <- c()
for (i in 1:nrow(R_complete)) {
  print(i)
  
  # construct range info
  x <- GRanges(R_complete$pos[i])
  
  overlap <- findOverlaps(x, all_transcripts, type = 'within')
  
  if (length(overlap) == 0) {
    non_tx_ind <- c(non_tx_ind, i)
    next
  }
  
  snp_loc <- c()
  # incorrect_mapping <- c()
  for (j in 1:length(overlap)) {
    tx_id <- names(all_transcripts[overlap@to[j]])
    
    tx <- all_transcripts[tx_id]
    
    if (tx@strand@values == '+') {
      boundary_5 <- tx@ranges@start
      boundary_3 <- tx@ranges@start + tx@ranges@width - 1
    } else {
      boundary_3 <- tx@ranges@start
      boundary_5 <- tx@ranges@start + tx@ranges@width - 1
    }
    
    prime5_dist <- abs(x@ranges@start - boundary_5)
    prime3_dist <- abs(boundary_3 - x@ranges@start)
    
    if (prime5_dist < prime3_dist) {
      utr_dist <- data.frame(distance = prime5_dist, end = "5' end")
    } else if (prime5_dist > prime3_dist) {
      utr_dist <- data.frame(distance = prime3_dist, end = "3' end")
    } else if (prime5_dist == prime3_dist) {
      next
    }
    
    snp_loc <- rbind(snp_loc, utr_dist)
  }
  snp_location[[i]] <- snp_loc 
  
}
save(list = c('snp_location', 'non_tx_ind'), 
     file = 'Desktop/circ_gene/data/sig_snp_loc_2.RData')

load('Desktop/circ_gene/data/sig_snp_loc_2.RData')
sig_snp_loc_df <- c()
snp_location_clean <- snp_location[-non_tx_ind]
for (i in 1:length(snp_location_clean)) {
  sig_snp_loc_df <- rbind(sig_snp_loc_df, snp_location_clean[[i]])
}
sig_snp_loc_df2 <- cbind(sig_snp_loc_df, group = 'Circ-SNPs')

random_snp_loc_dir <- 'Desktop/circ_gene/data/random_snp_loc'
random_snp_loc_files <- list.files(random_snp_loc_dir)

random_snp_loc <- c()
for (i in 1:length(random_snp_loc_files)) {
  load(file.path(random_snp_loc_dir, random_snp_loc_files[i]))
  snp_location_clean <- snp_location[-non_tx_ind]
  random_snp_loc <- c(random_snp_loc, snp_location_clean)
}

random_snp_loc_df <- c()
for (i in 1:length(random_snp_loc)) {
  random_snp_loc_df <- rbind(random_snp_loc_df, random_snp_loc[[i]])
}
save(random_snp_loc_df, file = 'Desktop/circ_gene/data/random_snp_loc_df.RData')

load('Desktop/circ_gene/data/random_snp_loc_df.RData')
random_snp_loc_df2 <- cbind(random_snp_loc_df, group = 'random SNPs')

set.seed(25)
remaining_random_ind <- seq(1, nrow(random_snp_loc_df2))
random_smooth_5_DF <- c()
random_smooth_3_DF <- c()
threshold_max <- 0
iterations <- 100
spacer_vec <- vector('numeric', length = iterations)
spacer_folds <- 10
for (i in 1:iterations) {
  random_ind <- sample(x = remaining_random_ind, size = nrow(sig_snp_loc_df2), replace = FALSE)
  random_snp_loc_df3 <- random_snp_loc_df2[remaining_random_ind[random_ind],]

  snp_loc_df <- rbind(sig_snp_loc_df2, random_snp_loc_df3)
  threshold <- max(snp_loc_df$distance)
  if (threshold > threshold_max) {
    threshold_max <- threshold
  }
  # padding <- threshold / 10
  snp_loc_df_filtered <- snp_loc_df[snp_loc_df$distance<=threshold,]
  
  # modify distances
  spacer <- threshold * spacer_folds
  spacer_vec[i] <- spacer
  mod5prime <- snp_loc_df_filtered[snp_loc_df_filtered$end=="5' end",]
  mod3prime <- snp_loc_df_filtered[snp_loc_df_filtered$end=="3' end",]
  mod3prime$distance <- spacer-mod3prime$distance
  
  modDF <- rbind(mod5prime,mod3prime)
  
  modDF_random <- modDF[modDF$group == 'random SNPs',]
  random_density <- density(modDF_random$distance, bw = 1/(range(modDF_random$distance)[2]))
  random_density_df <- data.frame(x = random_density$x, y = random_density$y, group = 'random SNPs')
  mod5prime_random <- mod5prime[mod5prime$group == 'random SNPs',]
  mod5prime_random_max <- max(mod5prime_random$distance)
  random_density_5 <- data.frame(x = random_density$x[random_density$x <= mod5prime_random_max],
                                 y = random_density$y[random_density$x <= mod5prime_random_max])
  random_smooth_5 <- lowess(x=random_density_5$x, y=random_density_5$y)
  random_smooth_5_order <- order(random_smooth_5$x, decreasing = FALSE)
  random_smooth_5_df <- data.frame(x = random_smooth_5$x[random_smooth_5_order],
                                   y = random_smooth_5$y[random_smooth_5_order])
  random_smooth_5_DF <- rbind(random_smooth_5_DF, cbind(random_smooth_5_df, batch = i))
  
  mod3prime_random <- mod3prime[mod3prime$group == 'random SNPs',]
  mod3prime_random_min <- min(mod3prime_random$distance)
  random_density_3 <- data.frame(x = random_density$x[random_density$x >= mod3prime_random_min],
                                 y = random_density$y[random_density$x >= mod3prime_random_min])
  random_smooth_3 <- lowess(x=random_density_3$x, y=random_density_3$y)
  random_smooth_3_order <- order(random_smooth_3$x, decreasing = FALSE)
  random_smooth_3_df <- data.frame(x = random_smooth_3$x[random_smooth_3_order],
                                   y = random_smooth_3$y[random_smooth_3_order])
  random_smooth_3_DF <- rbind(random_smooth_3_DF, cbind(random_smooth_3_df, batch = i))
  
}

save(list = c('random_smooth_5_DF', 'random_smooth_3_DF'),
     file = 'Desktop/circ_gene/data/random_smooth_1000.RData')

load('Desktop/circ_gene/data/random_smooth_1000.RData')

random_low_5 <- Inf
random_high_5 <- -Inf
random_low_5_ind <- 0
random_high_5_ind <- 0
for (i in unique(random_smooth_5_DF$batch)) {
  batch_data <- random_smooth_5_DF[random_smooth_5_DF$batch == i,]
  batch_representative <- batch_data$y[1]
  if (batch_representative > random_high_5) {
    random_high_5_ind <- i
    random_high_5 <- batch_representative
  }
  
  if (batch_representative < random_low_5) {
    random_low_5_ind <- i
    random_low_5 <- batch_representative
  }
}

random_5_high <- random_smooth_5_DF[random_smooth_5_DF$batch==random_high_5_ind,]
random_5_low <- random_smooth_5_DF[random_smooth_5_DF$batch==random_low_5_ind,]

random_low_3 <- Inf
random_high_3 <- -Inf
random_low_3_ind <- 0
random_high_3_ind <- 0
for (i in unique(random_smooth_3_DF$batch)) {
  batch_data <- random_smooth_3_DF[random_smooth_3_DF$batch == i,]
  batch_representative <- batch_data$y[nrow(batch_data)]
  if (batch_representative > random_high_3) {
    random_high_3_ind <- i
    random_high_3 <- batch_representative
  }
  
  if (batch_representative < random_low_3) {
    random_low_3_ind <- i
    random_low_3 <- batch_representative
  }
}

random_3_high <- random_smooth_3_DF[random_smooth_3_DF$batch==random_high_3_ind,]
random_3_low <- random_smooth_3_DF[random_smooth_3_DF$batch==random_low_3_ind,]

# modify distances
spacer_max <- threshold_max * spacer_folds
mod5prime_sig <- sig_snp_loc_df2[sig_snp_loc_df2$end=="5' end",]
mod3prime_sig <- sig_snp_loc_df2[sig_snp_loc_df2$end=="3' end",]
min_3prime_dist_sig <- min(mod3prime_sig$distance)
mod3prime_sig$distance <- spacer_max-mod3prime_sig$distance

modDf_circ <- rbind(mod5prime_sig,mod3prime_sig)

sig_density <- density(modDf_circ$distance, bw = 1/(range(modDf_circ$distance)[2]))
sig_density_df <- data.frame(x = sig_density$x, y = sig_density$y, group = 'Circ-SNPs')
mod5prime_sig <- mod5prime_sig[mod5prime_sig$group == 'Circ-SNPs',]
mod5prime_sig_max <- max(mod5prime_sig$distance)
sig_density_5 <- data.frame(x = sig_density$x[sig_density$x <= mod5prime_sig_max],
                            y = sig_density$y[sig_density$x <= mod5prime_sig_max])
sig_smooth_5 <- lowess(x=sig_density_5$x, y=sig_density_5$y)
sig_smooth_5_order <- order(sig_smooth_5$x, decreasing = FALSE)
sig_smooth_5_df <- data.frame(x = sig_smooth_5$x[sig_smooth_5_order],
                              y = sig_smooth_5$y[sig_smooth_5_order],
                              batch = 1)

mod3prime_sig <- mod3prime_sig[mod3prime_sig$group == 'Circ-SNPs',]
mod3prime_sig_min <- min(mod3prime_sig$distance)
sig_density_3 <- data.frame(x = sig_density$x[sig_density$x >= mod3prime_sig_min],
                            y = sig_density$y[sig_density$x >= mod3prime_sig_min])
sig_smooth_3 <- lowess(x=sig_density_3$x, y=sig_density_3$y)
sig_smooth_3_order <- order(sig_smooth_3$x, decreasing = FALSE)
sig_smooth_3_df <- data.frame(x = sig_smooth_3$x[sig_smooth_3_order],
                              y = sig_smooth_3$y[sig_smooth_3_order],
                              batch = 1)

sig_true_rows <- sig_density_df$x <= mod5prime_sig_max | sig_density_df$x >= mod3prime_sig_min
sig_true_rows_ind <- seq(1, nrow(sig_density_df))[sig_true_rows]
sig_density_df_clean <- sig_density_df[sig_true_rows_ind,]

farthest_loc <- c(max(sig_smooth_3_df$x), max(random_3_high$x), max(random_3_low$x))
names(farthest_loc) <- c('sig_smooth_3_df', 'random_3_high', 'random_3_low')
loc_order <- order(farthest_loc, decreasing = TRUE)
spacers <- c(spacer_max, spacer_vec[random_3_high$batch[1]], spacer_vec[random_3_low$batch[1]])
spacer_final <- spacers[loc_order[1]]

for (i in 2:length(loc_order)) {
  reference <- get(names(farthest_loc)[loc_order[1]])
  raw_curve <- get(names(farthest_loc)[loc_order[i]])
  
  x_offest <- max(reference$x) - max(raw_curve$x)
  new_curve <- data.frame(x = raw_curve$x + x_offest, y = raw_curve$y, batch = raw_curve$batch)
  assign(names(farthest_loc)[loc_order[i]], new_curve)
}

tick_breaks <- c(0, max(c(sig_smooth_5_df$x, random_5_high$x, random_5_low$x)),
                 min(c(sig_smooth_3_df$x, random_3_high$x, random_3_low$x)) , 
                 spacer_final ) # - spacer_offset

tick_labels <- c("5' end", 
                 formatC(tick_breaks[2], format = "e", digits = 2),
                 formatC(tick_breaks[4] - tick_breaks[3], format = "e", digits = 2),
                 "3' end")

random_df <- rbind(rbind(cbind(random_5_high, group = 'high'), cbind(random_5_low, group = 'low')),
                   rbind(cbind(random_3_high, group = 'high'), cbind(random_3_low, group = 'low')))

ggplot() +
  geom_point(data = sig_density_df_clean, aes(x = x, y = y), colour = '#F8766D', alpha = .5) +
  geom_line(data = random_5_high, aes(x = x, y = y), colour = "blue", size = .2, alpha = .6) +
  geom_line(data = random_5_low, aes(x = x, y = y), colour = "blue", size = .2, alpha = .6) +
  geom_line(data = random_3_high, aes(x = x, y = y), colour = "blue", size = .2, alpha = .6) +
  geom_line(data = random_3_low, aes(x = x, y = y), colour = "blue", size = .2, alpha = .6) +
  geom_area(data = random_df, aes(x, y, fill = group), alpha = .2) +
  scale_fill_manual(values=c("blue", "grey")) + 
  geom_line(data = sig_smooth_5_df, aes(x = x, y = y), colour = "red", size = .5, alpha = .6) +
  geom_line(data = sig_smooth_3_df, aes(x = x, y = y), colour = "red", size = .5, alpha = .6) +
  scale_x_continuous(breaks=tick_breaks, labels= tick_labels) +
  theme(axis.text.x = element_text(color=c('#00A08A','#00A08A',"#FF0000","#FF0000"),
                                   size=10),
        axis.title = element_blank(), axis.text = element_text(size = 10),
        legend.position = 'none')

setwd('Desktop/circ_gene/plot/final/')
ggsave(
  'tx_boundary_dist.png',
  plot = last_plot(),
  width = 5,
  height = 3,
  units = 'in',
  device = "png",
  dpi = "print",
  bg = NULL,
)

sig_snp_loc_df_5 <- sig_snp_loc_df[sig_snp_loc_df$end == "5' end",]
sig_snp_loc_df_3 <- sig_snp_loc_df[sig_snp_loc_df$end == "3' end",]

random_snp_loc_df_5 <- random_snp_loc_df[random_snp_loc_df$end == "5' end",]
random_snp_loc_df_3 <- random_snp_loc_df[random_snp_loc_df$end == "3' end",]

summary(sig_snp_loc_df_5$distance)
summary(random_snp_loc_df_5$distance)

summary(sig_snp_loc_df_3$distance)
summary(random_snp_loc_df_3$distance)











