#params: udiffsorted (sorted list of tuples from BasicAnalysis.py of k-mers appearing more in upregulated genes), ddiffsorted (same for 
#downreg genes)
def diversity(udiffsorted, ddiffsorted):
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
