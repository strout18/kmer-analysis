def nuc2(updownreg, nuc, kmerlength):
        #updownreg = list of either genes more upregulated or genes more downregulated
        #gene_seq should be dictionary in form gene: sequence
        nuc = nuc.upper()
        avg_in_gene = []
        for x in updownreg.keys():
            giq = {} #gene in question
            kmeravg = []
            freqfunc.freq_compile_from_string(giq, gene_seq[x], kmerlength)
            for y in giq.keys():
                count = 0
                for z in y:
                    if z == nuc:
                        count += 1
                kmeravg.append(count)
            geneavg = sum(kmeravg) / len(kmeravg)
            avg_in_gene.append(geneavg)
        print ('avg: ', sum(avg_in_gene)/len(avg_in_gene))
        print (len(avg_in_gene))
