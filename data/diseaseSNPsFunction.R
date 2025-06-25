# return symbol and color for p-val
pValMark <- function(pVal) {
  mark <- c()
  markColor <- c()
  if (pVal>0.01 & pVal<0.05) {
    mark <- '*'
    markColor <- 'gray80'
  } else if(pVal<=0.01 & pVal>0.001) {
    mark <- '**'
    markColor <- 'gray60'
  } else if(pVal<=0.001 & pVal>0.0001) {
    mark <- '***'
    markColor <- 'gray40'
  } else if(pVal<=0.0001 & pVal>0.00001) {
    mark <- '****'
    markColor <- 'gray20'
  } else if(pVal<= 0.00001){
    mark <- '*****'
    markColor <- 'black'
  } else {
    mark <- ''
    markColor <- 'white'
  }
  
  return(list(M = mark,C = markColor))
}


# return q-val for a rsid
qVal <- function(rsid,GTEx) {
  GTExInd <- seq(1,nrow(GTEx))
  rsidInd <- GTExInd[rsid==GTEx$rs_id_dbSNP151_GRCh38p7]
  qVal <- c()
  if (isEmpty(rsidInd)==FALSE) {
    qVal <- GTEx[rsidInd,]$qval
  }
  return(qVal)
}

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

range2GRanges <- function(df) {
  require(GenomicRanges)
  require(IRanges)
  gr <- GenomicRanges::GRanges(
    seqnames = df[,1],
    ranges=IRanges(start = df[,2], end = df[,3])
  )
  return(gr)
}



