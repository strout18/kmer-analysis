fpkm = read.table("/Users/stella/Desktop/Genes.fpkm_table.txt")
rownames(fpkm) = fpkm[,1]
fpkm = fpkm[,-1]
colnames(fpkm) = c("0hr_rep0", "0hr_rep1", "4hr_rep0", "4hr_rep1")
fpkm = fpkm[-1,]
class(fpkm$`0hr_rep0`)
fpkm$`0hr_rep0` = as.numeric(as.character(fpkm$`0hr_rep0`))
fpkm$`0hr_rep1` = as.numeric(as.character(fpkm$`0hr_rep1`))
fpkm$`4hr_rep0` = as.numeric(as.character(fpkm$`4hr_rep0`))
fpkm$`4hr_rep1` = as.numeric(as.character(fpkm$`4hr_rep1`))
fpkm$ctrlavg = rowMeans(fpkm[,c(1,2)])
fpkm$testavg = rowMeans(fpkm[,c(3,4)])
fpkm$ctrlavg[fpkm$`0hr_rep0`==0 | fpkm$`0hr_rep1`==0] = fpkm$`0hr_rep1`[fpkm$`0hr_rep0`==0 | fpkm$`0hr_rep1`==0] + fpkm$`0hr_rep0`[fpkm$`0hr_rep0`==0 | fpkm$`0hr_rep1`==0]
fpkm$testavg[fpkm$`4hr_rep0`==0 | fpkm$`4hr_rep1`==0] = fpkm$`4hr_rep1`[fpkm$`4hr_rep0`==0 | fpkm$`4hr_rep1`==0] + fpkm$`4hr_rep0`[fpkm$`4hr_rep0`==0 | fpkm$`4hr_rep1`==0]
fpkm$foldchange = 0
fpkm$foldchange[fpkm$ctrlavg != 0 & fpkm$testavg !=0] = log2(fpkm$testavg[fpkm$ctrlavg!=0 & fpkm$testavg != 0]/fpkm$ctrlavg[fpkm$ctrlavg!=0 & fpkm$testavg!=0])
library(readxl)
e2g = read_xlsx("/Users/stella/Desktop/Lin Lab/HAEC_P65_TNF_4hr_1_peaks_AllEnhancers_ENHANCER_TO_GENE.xlsx")
fpkm$geneid = rownames(fpkm)
final = merge(fpkm, e2g, by.x = "geneid", by.y = "CLOSEST_GENE")
final = subset(final, select = c("geneid", "foldchange", "CHROM", "START", "STOP"))
final = final[,c(1,3,4,5,2)]
#write.table(final, file = "/Users/stella/Desktop/Lin Lab/R_FinalOut.txt", quote = FALSE, sep = "\t", row.names = FALSE)
max10 = fpkm[,c(1:4)]
max10$max = apply(max10, 1, max)
max10 = max10[max10$max > 10, ]
max10$geneid = row.names(max10)
max10final = merge(final, max10, by.x = "geneid", by.y= "geneid")
max10final = max10final[,c(1:5)]
downreg = max10final[,c(1, 6, 7)]
downreg = downreg[downreg$foldchange < 0, ]
max10final = max10final[abs(max10final$foldchange) > 1, ]
write.table(max10final, file = "/Users/stella/Desktop/Lin Lab/R_FinalFiltered.txt", quote = FALSE, sep = "\t", row.names = FALSE)
