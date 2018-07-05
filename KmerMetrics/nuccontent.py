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
