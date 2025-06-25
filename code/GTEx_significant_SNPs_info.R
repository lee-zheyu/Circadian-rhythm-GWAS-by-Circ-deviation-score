library(rsnps)
library(dplyr)

# add info for significant SNPs
load('Desktop/circ_gene/data/significant_SNP.RData')

# add rs id
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
snp_pos <- unlist(lapply(R$pos, function(x){
  gsub('chr','',x)
}))
rs_id <- vector('list', length = nrow(R))
for (i in 1:length(rs_id)) {
  # print(i)
  snp_grange <- GRanges(snp_pos[i])
  overlap <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, snp_grange)
  if (length(overlap)==0) {
    next
  } else {
    rs_id[[i]] <- overlap@elementMetadata@listData$RefSNP_id
  }
}

rs_id_len <- vector('numeric', length = length(rs_id))
for (i in 1:length(rs_id_len)) {
  rs_id_len[i] <- length(rs_id[[i]])
}

no_rs_pos <- snp_pos[rs_id_len==0]

missing_rs_id <- c('rs193097739', 'rs184192643', 'rs138892563', 'rs146743744')

no_rs_ind <- seq(1,length(rs_id))[rs_id_len==0]
for (i in 1:length(no_rs_ind)) {
  rs_id[[no_rs_ind[i]]] <- missing_rs_id[i]
}

rs_id <- rs_id[rs_id_len!=0]
R <- R[rs_id_len!=0,]

rs_id_len_recheck <- vector('numeric', length = length(rs_id))
for (i in 1:length(rs_id_len_recheck)) {
  rs_id_len_recheck[i] <- length(rs_id[[i]])
}
unique(rs_id_len_recheck)

RS_id <- unlist(rs_id)

R_id <- cbind(RS_id,R)

# manually map SNPs to genes
library(GenomicRanges)
library(writexl)

source('Desktop/circ_gene/data/diseaseSNPsFunction.R')

gtf <- read.table('C:/Users/zheyuli/research/pathway_SNP/data/gencode.v45.chr_patch_hapl_scaff.annotation.gtf.gz', 
                  header = F, stringsAsFactors = F, sep = '\t')
gtfGene <- gtf[gtf[,3]=="gene", c(1,4,5)]
geneNames <- unlist(lapply(gtf[gtf[,3]=="gene", 9], function(x) {
  y <- strsplit(x, ';')[[1]][1]
  gsub('gene_id ','',y)
}))

# convert genes to GRanges
gtfGranges <- range2GRanges(gtfGene)
names(gtfGranges) <- geneNames

Gene <- vector('character',length = nrow(R_id))
gene_num <- vector('numeric', length = nrow(R_id))
missing_ind <- c()
for (i in 1:length(Gene)) {
  print(i)
  
  # convert snp coor to GRange
  snpCoor <- R_id$pos[i]
  snpsRanges <- string2range(snpCoor, delim=":", region=FALSE)
  snpsGranges <- range2GRanges(snpsRanges)
  names(snpsGranges) <- snpCoor
  
  # find overlap of snps and genes
  overlap <- GenomicRanges::findOverlaps(snpsGranges, gtfGranges, maxgap = -1, type = 'any')
  
  # make vector of SNPs to genes
  hit_gene <- names(gtfGranges[overlap@to])
  
  if (length(hit_gene)==0) {
    missing_ind <- c(missing_ind,i)
    next
  }
  
  hit_gene_unversioned <- unlist(lapply(hit_gene, function(x){
    gsub('\\.\\d+','',x)
  }))
  
  if (length(hit_gene)==1) {
    Gene[i] <- hit_gene_unversioned
    gene_num[i] <- 1
  } else {
    geneIDs <- unique(hit_gene_unversioned)
    Gene[i] <- paste(geneIDs,collapse = ';')
    gene_num[i] <- length(geneIDs)
  }
}


R_genes <- cbind(R_id, Gene)

# check if gene is circadian
load('Desktop/circ_gene/data/deviation_by_tissue_weighted_complete.RData')
R_genes <- R_genes[order(R_genes$tissue),]

# global
circadian_genes <- c()
for (i in 1:length(deviation_by_tissue)) {
  tissue_circadian_genes <- colnames(deviation_by_tissue[[i]])
  # tissue_circadian_genes <- unlist(lapply(tissue_circadian_genes_versioned, function(x){
  #   gsub('\\.\\d+','',x)
  # }))
  circadian_genes <- c(circadian_genes, tissue_circadian_genes)
}
circadian_genes <- unique(circadian_genes)

is_circadian_global <- vector('numeric', length = nrow(R_genes))
for (j in 1:length(is_circadian_global)) {
  hit_gene <- R_genes$Gene[j]
  if (hit_gene=='') {
    is_circadian_global[j] <- 0
    next
  } else {
    hit_gene_split_version <- unlist(strsplit(hit_gene,';'))
    hit_gene_split <- unlist(lapply(hit_gene_split_version, function(x){
      gsub('\\.\\d+','',x)
    }))
    if (sum(hit_gene_split %in% circadian_genes)==0) {
      is_circadian_global[j] <- 1
    } else {
      is_circadian_global[j] <- 2
    }
  }
}
R_genes$is_circadian <- is_circadian_global

# check if hit genes are tissue specific circadian genes
significant_tissues <- unique(R_genes$tissue)
is_circadian_tissue <- c()
for (i in 1:length(significant_tissues)) {
  tissue <- significant_tissues[i]
  circadian_genes <- colnames(deviation_by_tissue[[tissue]])
  R_sub <- R_genes[R_genes$tissue==tissue,]
  is_circadian <- vector('numeric',length = nrow(R_sub))
  for (j in 1:length(is_circadian)) {
    hit_gene <- R_sub$Gene[j]
    if (hit_gene=='') {
      is_circadian[j] <- 0
      next
    } else {
      hit_gene_split_version <- unlist(strsplit(hit_gene,';'))
      hit_gene_split <- unlist(lapply(hit_gene_split_version, function(x){
        gsub('\\.\\d+','',x)
      }))
      if (sum(hit_gene_split %in% circadian_genes)==0) {
        is_circadian[j] <- 1
      } else {
        is_circadian[j] <- 2
      }
    }
  }
  is_circadian_tissue <- c(is_circadian_tissue,is_circadian)
}
R_genes <- cbind(R_genes,is_circadian_tissue)

# check whether the SNP is eQTL
eQTL_dir <- "Desktop/circ_gene/data/GTEx_Analysis_v8_eQTL"
files_e <- list.files(eQTL_dir)
files_sig_var_e <- files_e[grepl('signif_variant_gene_pairs',files_e)]

sig_var_id_e <- c()
for (i in 1:length(files_sig_var_e)) {
  tissue_data <- read.delim(file.path(eQTL_dir,files_sig_var_e[i]), header = T)
  sig_var_id_e <- union(sig_var_id_e, tissue_data$variant_id)
}

is_eQTL <- vector('numeric',length = nrow(R))
for (i in 1:length(is_eQTL)) {
  variant_id <- R$variant_id[i]
  if (is.element(variant_id,sig_var_id_e)) {
    is_eQTL[i] <- 1
  }
}

sQTL_dir <- "Desktop/circ_gene/data/GTEx_Analysis_v8_sQTL"
files_s <- list.files(sQTL_dir)
files_sig_var_s <- files_s[grepl('.v8.sqtl_signifpairs.txt.gz',files_s)]

sig_var_id_s <- c()
for (i in 1:length(files_sig_var_s)) {
  tissue_data <- read.delim(file.path(sQTL_dir,files_sig_var_s[i]), header = T)
  sig_var_id_s <- union(sig_var_id_s, tissue_data$variant_id)
}

is_sQTL <- vector('numeric',length = nrow(R))
for (i in 1:length(is_sQTL)) {
  variant_id <- R$variant_id[i]
  if (is.element(variant_id,sig_var_id_s)) {
    is_sQTL[i] <- 1
  }
}

edQTL_dir <- "Desktop/circ_gene/data/GTEx_Analysis_v8_edQTL"
files_ed <- list.files(edQTL_dir)
files_sig_var_ed <- files_ed[grepl('signif_variant_',files_ed)]

sig_var_id_ed <- c()
for (i in 1:length(files_sig_var_ed)) {
  tissue_data <- read.delim(file.path(edQTL_dir,files_sig_var_ed[i]), header = T)
  sig_var_id_ed <- union(sig_var_id_ed, tissue_data$variant_id)
}

is_edQTL <- vector('numeric',length = nrow(R))
for (i in 1:length(is_edQTL)) {
  variant_id <- R$variant_id[i]
  if (is.element(variant_id,sig_var_id_ed)) {
    is_edQTL[i] <- 1
  }
}

R_QTL <- cbind(R_genes, is_eQTL, is_sQTL, is_edQTL)

# add gene symbols and types
library(biomaRt)

ensembl_human <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# collect gene ids
gene_ids <- c()
for (i in 1:nrow(R_QTL)) {
  g <- R_QTL$Gene[i]
  if (g == '') {
    next
  }
  gene_ids <- c(gene_ids,unlist(strsplit(g,';')))
}
gene_ids <- unique(gene_ids)

# convert gene ids to gene symbols
gene_info <- getBM(attributes = c("ensembl_gene_id","external_gene_name",
                                  "gene_biotype"),
                   filters = "ensembl_gene_id", values = gene_ids,
                   mart = ensembl_human)
colnames(gene_info) <- c('gene_id','gene_symbol','gene_type')
gene_symbols <- vector('character', length = nrow(R_QTL))
gene_types <- gene_symbols
for (i in 1:nrow(R_QTL)) {
  hit_gene <- R_QTL$Gene[i]
  if (hit_gene=='') {
    next
  } else {
    hit_gene_split_version <- unlist(strsplit(hit_gene,';'))
    hit_gene_split <- unlist(lapply(hit_gene_split_version, function(x){
      gsub('\\.\\d+','',x)
    }))
    symbol <- vector('character', length = length(hit_gene_split))
    types <- symbol
    for (j in 1:length(hit_gene_split)) {
      ind <- match(hit_gene_split[j],gene_info$gene_id)
      if (is.na(ind)) {
        symbol[j] <- 'missing'
        types[j] <- 'missing'
      } else {
        symbol[j] <- gene_info$gene_symbol[ind]
        types[j] <- gene_info$gene_type[ind]
      }
    }
    gene_symbols[i] <- paste0(symbol, collapse = ';')
    gene_types[i] <- paste0(types, collapse = ';')
  }
}

R_complete <- cbind(R_QTL,gene_symbols,gene_types)

save('R_complete', file = 'Desktop/circ_gene/data/R_weighted_complete2.RData')

load('Desktop/circ_gene/data/R_weighted_complete2.RData')
write_xlsx(R_complete[,-c(7, 14:18)],
           "Desktop/circ_gene/data/circ_SNPs_info.xlsx")
