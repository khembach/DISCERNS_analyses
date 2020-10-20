# 07.07.2016 Katharina Hembach
## This script reads the human_ce_location_hg19_GRCh38.txt file and export a gff file with the GRCh38 locations of the cryptic exons from Ling et al.


dat <- read.table("human_ce_location_hg19_GRCh38.txt", header=TRUE, stringsAsFactors=FALSE)
chromosomes <- unlist(lapply(strsplit(dat$location_GRCh38, split=":"), function(i) i[1]))
positions <- unlist(lapply(strsplit(dat$location_GRCh38, split=":"), function(i) i[2]))
starts <- as.integer(unlist( lapply(strsplit(positions, "-"), function(i) i[1]) ) )
ends <- as.integer(unlist( lapply(strsplit(positions, "-"), function(i) i[2]) ) )
cryptic_exons <- data.frame(name=dat$geneSymbol, chrom=chromosomes, start=starts, end=ends)

library(rtracklayer)
library(GenomicRanges)

targetRanges <- IRanges(start=cryptic_exons$start, end=cryptic_exons$end)
targetTrack <- GRanges(seqnames=cryptic_exons$chrom, ranges=targetRanges )

export(targetTrack, "human_ce_GRCh38.gff")
write.table(cryptic_exons, "ce_GRCH38_start_end.txt", quote=FALSE, row.names=FALSE, sep="\t")