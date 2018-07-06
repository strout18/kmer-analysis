
# coding: utf-8

# In[ ]:

import Reusable_Kmer_Generator as freqfunc
def enhancer_diversity(kmerlength, gene_seq, gene_fc,  show_updownkmers = True, nuc = 'T'):
    print ('Kmers of length', kmerlength)
    freq_dict = {}
    for key in gene_seq:
        freqfunc.freq_compile_from_string(freq_dict, gene_seq[key], kmerlength)
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
    #NOTE: THIS IS IMPERFECT BC IT GOES BY TOP CUTOFF OF GIVEN LIST, AND FORGOES THOSE OF OPP LIST
    #POSSIBLE FIX: DO IT FOR ENTIRE LENGTH OF LIST
    #LST MUST BE SORTED
    upregdiff = {}
    downregdiff = {}
    def freq_compare(lst1, lst2, dct1, dct2): #IF LST = UPREG, OPPDICT MUST EQUAL DOWNREG
        dct1sum = sum(dct1.values())
        dct2sum = sum(dct2.values())
        for i in range(len(lst1)): #loops through upreg
            kmer = lst1[i][0]
            pct1 = int(lst1[i][1]) / dct1sum
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
    count = 0
    for x in upregdiffsorted:
        if len(set(x[0])) > 2:
            count += 1
    print (count, 'diverse k-mers out of', len(upregdiffsorted), 'that comprise a greater percent of all kmers in upregulated genes.')
    print ("%s / %s = %s" % (count, len(upregdiffsorted), count/ len(upregdiffsorted)))
    count = 0
    for x in downregdiffsorted:
        if len(set(x[0])) > 2:
            count += 1
    print (count, 'diverse k-mers out of', len(downregdiffsorted), 'that comprise a greater percent of all kmers in downregulated genes.')
    try:
        print ("%s / %s = %s" % (count, len(downregdiffsorted), count/ len(downregdiffsorted)))  
    except ZeroDivisionError:
        pass
    def nuc_content(sortedlst, cutoff, nuc): #dct = upregdiff or downregdiff, nuc = 'A', 'C' ,'T', 'G' 
        #UPREG SHOULD BE SORTED GREATEST TO LEAST, DOWNREG LEAST TO GREATEST
        nuc = nuc.upper()
        nuccount = []
        for i in range(cutoff):
            count = 0
            for x in sortedlst[i][0]:
                if x.upper() == nuc:
                    count +=1
                nuccount.append(count)
        avg = sum(nuccount) / len(nuccount)
        print ('Average number of %s in %s-mer: %s' %(nuc, kmerlength, avg))
    numcutoff = 200
    if len(upregdiffsorted) < 200 or len(downregdiffsorted) < 200:
        if len(upregdiffsorted) > len(downregdiffsorted):
            numcutoff = len(downregdiffsorted)
        else:
            numcutoff = len(upregdiffsorted)
    print ('Cutoff:', numcutoff)
    print ('Upregulated')
    nuc_content(upregdiffsorted, numcutoff, nuc)
    print ('Downregulated')
    nuc_content(downregdiffsorted, numcutoff, nuc)
    outputlst = [freq_dict, upregfreq, downregfreq, upregdiff, downregdiff] 
    return outputlst

