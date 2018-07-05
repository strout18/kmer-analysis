#params: udiffsorted (sorted list of tuples from BasicAnalysis.py of k-mers appearing more in upregulated genes), ddiffsorted (same for 
#downreg genes)
def diversity(udiffsorted, ddiffsorted):
  count = 0
  for x in udiffsorted:
      if len(set(x[0])) > 2: #if set of unique nucs in k-mer > 2
          count += 1
  print (count, 'diverse k-mers out of', len(udiffsorted), 'that comprise a greater percent of all kmers in upregulated genes.')
  try:
    print ("%s / %s = %s" % (count, len(udiffsorted), count/ len(udiffsorted))) #return # diverse kmers/ all kmers
  except ZeroDivisionError:
    pass
  count = 0
  for x in ddiffsorted:
      if len(set(x[0])) > 2:
          count += 1
  print (count, 'diverse k-mers out of', len(ddiffsorted), 'that comprise a greater percent of all kmers in downregulated genes.')
  try:
      print ("%s / %s = %s" % (count, len(ddiffsorted), count/ len(ddiffsorted)))  
  except ZeroDivisionError:
      pass
