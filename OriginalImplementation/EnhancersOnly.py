peakfasta = "/xxxxxx/enhancersonly.fa"
rtable = "/xxxxxx/R_EnhancersFiltered.txt"
import csv
gene2enhancer = {}
enhancer2gene = {}
gene_fc = {}
with open(rtable) as f:
    data = csv.reader(f, delimiter = "\t")
    for row in data:
        if row[0] != "geneid":
            i = 0
            nonunique = True
            while nonunique:
                geneid = str(i) + row[0]
                if geneid not in gene2enhancer.keys():
                    nonunique = False
                i += 1
            gene2enhancer[geneid] = row[1] + ":" + str(row[2]) + "-" + str(row[3])
            enhancer2gene[gene2enhancer[geneid]] = geneid
            gene_fc[geneid] = str(row[4])
    f.close()
j = 0
for key in gene2enhancer.keys():
    print (key, gene2enhancer[key])
    if j>10:
        break
    j += 1
j = 0
for key in enhancer2gene.keys():
    print (key, enhancer2gene[key])
    if j>10:
        break
    j += 1
gene_seq = {}
with open(peakfasta) as f:
    for line in f:
        if line[0] == ">":
            chromid = line[1:len(line)-1] #removes newline char
            if chromid in enhancer2gene:
                gene_seq[enhancer2gene[chromid]] = next(f).upper()
    f.close()
i = 0
for key in gene_seq:
    print (key, gene_seq[key])
    i += 1
    if i>10:
        break
import Reusable_Kmer_Generator as freqfunc
import EnhancerDiversity as ed
lstall5mers = ed.enhancer_diversity(5, gene_seq, gene_fc, False, 'A')
lstall6mers = ed.enhancer_diversity(6, gene_seq, gene_fc, False, 'A')
