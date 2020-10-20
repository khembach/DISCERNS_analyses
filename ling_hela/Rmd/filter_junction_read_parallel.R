filter_jr_parallel <- function(bam, yield_size = 200000, 
                                  lib_type = "PE", stranded = "reverse") {
  if (!lib_type %in% c("SE", "PE")) {
    stop('Parameter lib_type has to be either "SE" or "PE".')
  }
  if (!stranded %in% c("unstranded", "forward", "reverse")) {
    stop('Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  }
  bf <- BamFile(bam, yieldSize = yield_size)
  
  YIELD <- function(x, ...) {
    flag0 <- scanBamFlag(isProperPair = TRUE)
    param <- ScanBamParam(what = c("qname", "flag"), flag = flag0)
    readGAlignments(x, param = param)
  }
  
  YIELD_SE <- function(x, ...) {
    param <- ScanBamParam(what = c("qname", "flag"))
    readGAlignments(x, param = param)
  }
  
  MAP <- function(reads, stranded, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag,
                               bitnames = c("isMinusStrand", "isFirstMateRead"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    mcols(reads)$isFirstMateRead <- bfbm[ , "isFirstMateRead"]
    
    if (stranded == "forward") {
      ## isMinusStrand==0 & isFirstMateRead==0 --> "-"
      ## isMinusStrand==0 & isFirstMateRead==1 --> "+"
      ## isMinusStrand==1 & isFirstMateRead==0 --> "+"
      ## isMinusStrand==1 & isFirstMateRead==1 --> "-" 
      strand(reads) <- ifelse(mcols(reads)$isMinusStrand +
                                mcols(reads)$isFirstMateRead == 0, "-",
                              ifelse(mcols(reads)$isMinusStrand +
                                       mcols(reads)$isFirstMateRead == 1, "+",
                                     "-"))
      reads
    } else {  # reverse
      ## isMinusStrand==0 & isFirstMateRead==0 --> "+"
      ## isMinusStrand==0 & isFirstMateRead==1 --> "-"
      ## isMinusStrand==1 & isFirstMateRead==0 --> "-"
      ## isMinusStrand==1 & isFirstMateRead==1 --> "+"
      strand(reads) <- ifelse(mcols(reads)$isMinusStrand +
                                mcols(reads)$isFirstMateRead == 0, "+",
                              ifelse(mcols(reads)$isMinusStrand +
                                       mcols(reads)$isFirstMateRead == 1, "-",
                                     "+"))
      reads
    }
  } 
  
  MAP_SE_REVERSE <- function(reads, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag, bitnames = c("isMinusStrand"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    strand(reads) <- ifelse(mcols(reads)$isMinusStrand == 0, "-", "+")
    reads
  }
  
  MAP_UNSTRANDED <- function(reads, ...) {
    reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
  }
  
  DONE <- function(value) {
    length(value) == 0L
  }
  
  if (lib_type == "PE") {
    if( stranded == "unstranded") {
      reduceByYield(bf, YIELD, MAP_UNSTRANDED, REDUCE = c, DONE, 
                    parallel = TRUE)
    } else {
      reduceByYield(bf, YIELD, MAP, REDUCE = c, DONE, parallel = TRUE, 
                    stranded = stranded)
    }
  } else if (stranded != "reverse") { # SE
    reduceByYield(bf, YIELD_SE, MAP_UNSTRANDED, REDUCE = c, DONE, 
                  parallel = TRUE)
  } else { 
    reduceByYield(bf, YIELD_SE, MAP_SE_REVERSE, REDUCE = c, DONE, 
                  parallel = TRUE)
  }
}




filter_jr <- function(bam, yield_size = 200000, 
                                  lib_type = "PE", stranded = "reverse") {
  if (!lib_type %in% c("SE", "PE")) {
    stop('Parameter lib_type has to be either "SE" or "PE".')
  }
  if (!stranded %in% c("unstranded", "forward", "reverse")) {
    stop('Parameter stranded has to be one of "unstranded", "forward" or "reverse".')
  }
  bf <- BamFile(bam, yieldSize = yield_size)
  
  YIELD <- function(x, ...) {
    flag0 <- scanBamFlag(isProperPair = TRUE)
    param <- ScanBamParam(what = c("qname", "flag"), flag = flag0)
    readGAlignments(x, param = param)
  }
  
  YIELD_SE <- function(x, ...) {
    param <- ScanBamParam(what = c("qname", "flag"))
    readGAlignments(x, param = param)
  }
  
  MAP <- function(reads, stranded, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag,
                               bitnames = c("isMinusStrand", "isFirstMateRead"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    mcols(reads)$isFirstMateRead <- bfbm[ , "isFirstMateRead"]
    
    if (stranded == "forward") {
      ## isMinusStrand==0 & isFirstMateRead==0 --> "-"
      ## isMinusStrand==0 & isFirstMateRead==1 --> "+"
      ## isMinusStrand==1 & isFirstMateRead==0 --> "+"
      ## isMinusStrand==1 & isFirstMateRead==1 --> "-" 
      strand(reads) <- ifelse(mcols(reads)$isMinusStrand +
                                mcols(reads)$isFirstMateRead == 0, "-",
                              ifelse(mcols(reads)$isMinusStrand +
                                       mcols(reads)$isFirstMateRead == 1, "+",
                                     "-"))
      reads
    } else {  # reverse
      ## isMinusStrand==0 & isFirstMateRead==0 --> "+"
      ## isMinusStrand==0 & isFirstMateRead==1 --> "-"
      ## isMinusStrand==1 & isFirstMateRead==0 --> "-"
      ## isMinusStrand==1 & isFirstMateRead==1 --> "+"
      strand(reads) <- ifelse(mcols(reads)$isMinusStrand +
                                mcols(reads)$isFirstMateRead == 0, "+",
                              ifelse(mcols(reads)$isMinusStrand +
                                       mcols(reads)$isFirstMateRead == 1, "-",
                                     "+"))
      reads
    }
  } 
  
  MAP_SE_REVERSE <- function(reads, ...) {
    reads <- reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
    
    bfbm <- bamFlagAsBitMatrix(mcols(reads)$flag, bitnames = c("isMinusStrand"))
    mcols(reads)$isMinusStrand <- bfbm[ , "isMinusStrand"]
    strand(reads) <- ifelse(mcols(reads)$isMinusStrand == 0, "-", "+")
    reads
  }
  
  MAP_UNSTRANDED <- function(reads, ...) {
    reads[ cigarOpTable(cigar(reads))[ ,"N"] > 0, ]
  }
  
  DONE <- function(value) {
    length(value) == 0L
  }
  
  if (lib_type == "PE") {
    if( stranded == "unstranded") {
      reduceByYield(bf, YIELD, MAP_UNSTRANDED, REDUCE = c, DONE, 
                    parallel = FALSE)
    } else {
      reduceByYield(bf, YIELD, MAP, REDUCE = c, DONE, parallel = FALSE, 
                    stranded = stranded)
    }
  } else if (stranded != "reverse") { # SE
    reduceByYield(bf, YIELD_SE, MAP_UNSTRANDED, REDUCE = c, DONE, 
                  parallel = FALSE)
  } else { 
    reduceByYield(bf, YIELD_SE, MAP_SE_REVERSE, REDUCE = c, DONE, 
                  parallel = FALSE)
  }
}

