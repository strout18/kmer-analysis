#params: udiffsorted (list of sorted tuples from BasicAnalysis.py for upreg genes), ddiffsorted (same for downreg genes) numcutoff (# of 
#k-mers from most frequent to cutoff to analyze nuc_content for), nuc ('A', 'C', 'T', or 'G')
def nuc_content(udiffsorted, ddiffsorted, nuc, numcutoff = 200): 
	nuc = nuc.upper()
	nuccount = []
	if len(udiffsorted) < 200 or len(ddiffsorted) < 200:
		#if either list has fewer k-mers than cutoff, it goes with lowest # as new cutoff for both lists
		if len(udiffsorted) > len(ddiffsorted):
			cutoff = len(ddiffsorted)
		else:
			cutoff = len(udiffsorted)
	else:
		cutoff = numcutoff
	for i in range(cutoff):
	    count = 0
	    for x in sortedlst[i][0]:
		if x.upper() == nuc:
		    count +=1
	    nuccount.append(count)
	avg = sum(nuccount) / len(nuccount)
	print ('Cutoff:', cutoff) 
	print ('Average number of %s in %s-mer: %s' %(nuc, kmerlength, avg))
