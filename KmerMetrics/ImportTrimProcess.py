peakfasta = "" #PUT YOUR FASTA FILE HERE
rtable = "" #PUT YOUR DATA TABLE HERE 
import csv 
gene2enhancer = {} #dict in form geneID: chr+start+stop
enhancer2gene = {} 
gene_fc = {} #dict in form geneID: fold change
with open(rtable) as f:
    data = csv.reader(f, delimiter = "\t")
    for row in data:
        if row[0] != "geneid":
            i = 0
            #the following code is put in place in case you have multiple sequences assigned to the same gene. It adds 1 to a number 
            #in front of the gene ID until a unique gene ID is generated so you can keep your sequences separate
            #for example: the first time the program sees the id "VCAM1" it assigns that particular sequence the gene ID "0VCAM1", 
            #the 2nd time the current sequence is assigned "1VCAM1", etc.
            nonunique = True
            while nonunique: #loops through until it generates a unique gene ID
                geneid = str(i) + row[0] 
                if geneid not in gene2enhancer.keys(): 
                    nonunique = False
                i += 1
            gene2enhancer[geneid] = row[1] + ":" + str(row[2]) + "-" + str(row[3])
            enhancer2gene[gene2enhancer[geneid]] = geneid
            gene_fc[geneid] = str(row[4])
    f.close()
j = 0
for key in gene2enhancer.keys(): #optional to get a sense of the output
    print (key, gene2enhancer[key])
    if j>10: 
        break
    j += 1
j = 0
for key in enhancer2gene.keys(): #optional to get a sense of the output
    print (key, enhancer2gene[key])
    if j>10:
        break
    j += 1
gene_seq = {} #dict in form geneID, sequence
with open(peakfasta) as f: #trims FASTA
    for line in f:
        if line[0] == ">":
            chromid = line[1:len(line)-1] #removes newline char
            if chromid in enhancer2gene:
                gene_seq[enhancer2gene[chromid]] = next(f).upper()
    f.close()
i = 0
for key in gene_seq: #optional to get a sense of the output
    print (key, gene_seq[key])
    i += 1
    if i>10:
        break
