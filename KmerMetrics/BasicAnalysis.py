#### SORT K-MERS BY GENE EXPRESSION ####
#params: kmerlength (int), gene_seq (dict in form geneID, sequence- COMPILE THIS FROM IMPORTTRIMPROCESS.PY!!!!), show_updownkmers (bool
#showing all k-mers and whether they appear more in up or downregulated genes)
def basic_analysis(kmerlength, gene_seq, show_updownkmers = True):
    print ('Kmers of length', kmerlength)
    freq_dict = {}
    for key in gene_seq:
        freqfunc.freq_compile_from_string(freq_dict, gene_seq[key], kmerlength) #REUSABLEKMERGENERATOR.PY SHOULD BE IMPORTED AS FREQFUNC
        #IN MAIN PROGRAM (you can import as something else but then must change 'freqfunc' in this module as well)
    sorted_all = freqfunc.fsortbycount(freq_dict)
    upreg = {} #dict of genes and their corresponding fold changes
    downreg = {}
    for key in gene_fc.keys():
        if float(gene_fc[key]) > 1.5:
            upreg[key] = float(gene_fc[key])
        elif float(gene_fc[key]) < -1.5:
            downreg[key] = float(gene_fc[key])
    print (len(upreg.keys()), len(downreg.keys()))
    upregfreq = {} #dict of kmers and corresponding kmer frequencies for upregulated genes
    for key in upreg.keys():
        freqfunc.freq_compile_from_string(upregfreq, gene_seq[key], kmerlength)
    upregfreq_sorted = freqfunc.fsortbycount(upregfreq) #list of tuples (kmer:freq) sorted by freq in upreg genes
    downregfreq = {}
    for key in downreg.keys():
        freqfunc.freq_compile_from_string(downregfreq, gene_seq[key], kmerlength)
    downregfreq_sorted = freqfunc.fsortbycount(downregfreq)
    upregdiff = {}
    downregdiff = {}
    #### COMPARE FREQUENCIES FOR EACH K-MER - SUBFUNCTION ####
    def freq_compare(lst1, lst2, dct1, dct2): #lst1, dct1 = sortedlst and freqdict for upreg genes, lst2, dict2 = same for downreg genes
        #up and down can be switched, but lst1 and dct1 must always have same gene regulation (up or down)
        dct1sum = sum(dct1.values())
        dct2sum = sum(dct2.values())
        for i in range(len(lst1)): #loops through upreg
            kmer = lst1[i][0]
            pct1 = int(lst1[i][1]) / dct1sum #percent k-mer frequency comprises of sum of all k-mer frequencies
            if kmer in dct2.keys():
                pct2 = int(dct2[kmer]) / dct2sum
            else:
                pct2 = 0
            if pct1 > pct2:
                upregdiff[kmer] = pct1
        for i in range(0, len(lst2)): #loops through downreg
            kmer = lst2[i][0]
            pct1 = int(lst2[i][1]) / dct2sum
            if kmer in dct1.keys():
                pct2 = int(dct1[kmer]) / dct1sum
            else:
                pct2 = 0
            if pct1 > pct2:
                downregdiff[kmer] = pct1
    freq_compare(upregfreq_sorted, downregfreq_sorted, upregfreq, downregfreq)
    upregdiffsorted = freqfunc.fsortbycount(upregdiff)
    downregdiffsorted = freqfunc.fsortbycount(downregdiff)
    print (len(upregdiffsorted))
    if show_updownkmers:
        for x in upregdiffsorted: #KMERS THAT COMPRISE A GREATER PERCENT OF ALL KMERS IN UPREGULATED GENES THAN 
            #THEY DO IN DOWNREGULATED GENES
            print ('Up:', x[0], ':', x[1])

        for x in downregdiffsorted:
            print ('Down:', x[0], ':', x[1])
    outputlst = [freq_dict, upregfreq, downregfreq, upregdiff, downregdiff] #useful to store in a var in main program 
    return outputlst
