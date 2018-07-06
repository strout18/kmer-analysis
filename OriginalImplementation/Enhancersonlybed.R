fasta = read.table("/xxxxxx/enhanceristhisright.bed")
colnames(fasta) = c("chr", "start", "stop")
fasta = fasta[,1:3]
fasta$all = paste(fasta$chr, fasta$start, fasta$stop, sep = " ")
max10 = read.table("/xxxxx/R_FinalFiltered.txt")
colnames(max10) = c("geneid", "CHROM", "START", "STOP", "foldchange")
max10 = max10[-1,]
max10$all = paste(max10$CHROM, max10$START, max10$STOP, sep = " ")
ultimate = merge(fasta, max10, by.x = "all", by.y = "all")
ultimate = ultimate[, c("geneid", "CHROM", "START", "STOP", "foldchange")]
write.table(ultimate, file = "/xxxxxx/R_EnhancersFiltered.txt", quote = FALSE, sep = "\t", row.names = FALSE)
