
# coding: utf-8

# In[1]:

#TO USE THESE FUNCTIONS YOU MUST HAVE THE FOLLOWING MODULES IMPORTED IN YOUR MAIN PROGRAM: RANDOM, RE, PLOTLY,
#PLOTLY.PLOTLY, PLOTLY.GRAPH_OBJS


# In[ ]:

#NOTES: UPDATE LISTMOSTFREQ SO WRITING TO A FILE IS OPTIONAL
#UPDATE PLOTLY SO IT EITHER PRINTS INLINE IN THE RUNNING NOTEBOOK OR OPENS URL TO NEW GRAPH
#UPDATE COMPILE FROM FILE SO IT ALSO IGNORES 'N'

# In[5]:

import re
import plotly.plotly as py
import plotly.graph_objs as go
import plotly
import random


# In[6]:

#READING FILES
def trim(untrimmedfasta, trimmedout): #trims FASTA file saved as a .txt by removing chr data and making all nucleotides capitalized
    with open(untrimmedfasta) as f:
        with open(trimmedout, 'w') as f1:
            for line in f:
                if line[0] != '>':
                    f1.write(line.upper())
            f1.close()
            f.close()


# In[3]:

#kmerdict = EMPTY dict 
def freq_compile_from_file(kmerdict, trimmedfasta, length):
#COMPILING DICTIONARY
    with open(trimmedfasta) as f:
        for line in f: #one line corresponds to one sequence at that chromosome location
            for i in range(0, len(line) - length):
                kmer = line[i:i+length] #makes k-mer out of next k letters
                if kmer not in kmerdict:
                    kmerdict[kmer] = 1
                else:
                    kmerdict[kmer] += 1
        f.close()
        print ("Finished!")


# In[ ]:

def freq_compile_from_string(kmerdict, string, length):
    for i in range(0, len(string) - length):
        kmer = string[i:i+length]
        if 'N' not in kmer.upper():
            if kmer not in kmerdict:
                kmerdict[kmer] = 1
            else:
                kmerdict[kmer] += 1
   # if i == len(string) - length - 1:
    #    print ("Finished!")


# In[11]:

#SORTING 
def fsortbycount(kdict): #compiles a sorted list (tuple) of all k-mers with highest count
    sortbycount = sorted(kdict.items(), key = lambda x: x[1], reverse = True)
    return sortbycount
def listall(sortedlst): #list all k-mers found
    for x in sortedlst:
        print (x[0], ':', str(x[1]))
def listmostfreq(sortedlst, cutoff, fmostfreq): #cutoff is number of kmers to display - e.g. top 25 means cutoff = 25
    i = 0
    with open(fmostfreq, 'w') as f:
        while i< cutoff:
            print (sortedlst[i][0], ':', str(sortedlst[i][1]))
            f.write(sortedlst[i][0])
            if i!= cutoff - 1:
                f.write('\n')
            i += 1
        f.close()


# In[14]:

#PARAMS: SORTED LIST OF KMERS, FREQ DICT, MOTIF, KMER LENGTH
def motif_analysis(sortedlst, kdict, motif, kmerlength): 
    sums = 0
    length = 0
    motif_kmers = []
    for i in range(0, len(motif)-kmerlength):
        kmer = ''
        for j in range(i, i + kmerlength):
            kmer = kmer + motif[j] #creates regular expression for k-mer from motif
        poss_kmers = [x[0] for x in sortedlst if re.search(kmer, x[0])] #loops thru list, finds k-mers that satisfy regex
        print (kmer)
        for x in poss_kmers:
            print (x, ':', kdict[x])
            sums += kdict[x]
            length += 1
            motif_kmers.append(x)
        print()
    print ('Mean:', sums/length) #prints mean frequency of nfkb kmers
    sorted_motif_kmers = sorted(motif_kmers, key = lambda x:kdict[x], reverse = True)
    return sorted_motif_kmers


# In[4]:

#PLOTLY
def plotwithmotiflines(sortedall, sorted_motif, kdict, graphtitle, url, gmode = 'lines'):
    names = [x[0] for x in sortedall] #plots kmer vs frequency
    values = [x[1] for x in sortedall]
    annlst = []
    shapelst = []
    count = 1
    for x in sorted_motif: #creates dotted line and label for nfkb kmers
        anndict = {}
        shapedict = {}
        anndict['x'] = x
        anndict['text'] = x
        anndict['showarrow'] = False 
        shapedict['type'] = 'line'
        shapedict['x0'] = x
        shapedict['x1'] = x
        shapedict['y0'] = int(kdict[x])
        shapedict['line'] = {
            'color': "#ff9400",
            'dash': "dash"
        }
        if count == 1: #positions line and label
            shapedict['y1'], anndict['y'] = 34500, 35000
        elif count == 2:
            shapedict['y1'], anndict['y'] = 1000, 700
        elif count ==3:
            shapedict['y1'], anndict['y'] = 29500, 30000
        elif count ==4:
            shapedict['y1'], anndict['y'] = 3000, 2500
        elif count == 5:
            shapedict['y1'], anndict['y'] = 24500, 25000
        else:
            shapedict['y1'], anndict['y'] = 5500, 5000
            count = 0
        count +=1
        annlst.append(anndict)
        shapelst.append(shapedict)


    plotly.tools.set_credentials_file(username='StellaTrout', api_key='sygEgpmtf2kyR8TYxVFV')

    trace0 = go.Scatter(
        x= names,
        y= values,
        mode = gmode
    )
    trace1 = go.Scatter( #plots motif kmers in different color than other points
        x = [x for x in sorted_motif],
        y = [kdict[x] for x in sorted_motif],
        mode = "markers",
        marker = {'color': "ff3d74"}
    )
    data = [trace0, trace1]
    layout = go.Layout( #formats axes
        xaxis = dict(
            showticklabels = False,
            title = "K-Mers"
        ),
        yaxis = dict(
            title = "Frequencies"
        ),
        title = graphtitle,
        annotations = annlst,
        shapes = shapelst
    )
    fig = go.Figure(data=data, layout=layout)
    py.iplot(fig, filename = url)


# In[16]:

#PLOTLY
def plotwithmotif(sortedall, sorted_motif, kdict, graphtitle, url, gmode = 'lines'):
    names = [x[0] for x in sortedall] #plots kmer vs frequency
    values = [x[1] for x in sortedall]

    trace0 = go.Scatter(
        x= names,
        y= values,
        mode = gmode
    )
    trace1 = go.Scatter( #plots nfkb kmers in different color than other points
        x = [x for x in sorted_motif],
        y = [kdict[x] for x in sorted_motif],
        mode = "markers",
        marker = {'color': "ff3d74"}
    )
    data = [trace0, trace1]
    layout = go.Layout( #formats axes
        xaxis = dict(
            showticklabels = False,
            title = "K-Mers"
        ),
        yaxis = dict(
            title = "Frequencies"
        ),
        title = graphtitle
    )
    fig = go.Figure(data=data, layout=layout)
    py.iplot(fig, filename = url)


# In[ ]:

#PLOTLY
def plotbasic(sortedall, graphtitle, url, gmode = 'lines'):
    names = [x[0] for x in sortedall] #plots kmer vs frequency
    values = [x[1] for x in sortedall]
    trace0 = go.Scatter(
        x= names,
        y= values,
        mode = gmode
    )
    data = [trace0]
    layout = go.Layout( #formats axes
        xaxis = dict(
            showticklabels = False,
            title = "K-Mers"
        ),
        yaxis = dict(
            title = "Frequencies"
        ),
        title = graphtitle
    )
    fig = go.Figure(data=data, layout=layout)
    py.iplot(fig, filename = url)


# In[8]:

def filerandseq(trimmedfasta, shuffledfile, dinucshuff, timestoshuffle):
    def dinucl(lst):
        i = 0
        dlst = []
        while i< len(lst):
            dlst.append(lst[i] + lst[i + 1]) #combines element with element after it
            i += 2
        return dlst

    def shuffleACopy(x):
            b = x[:] # make a copy of the keys
            random.shuffle(b) # shuffle the copy
            return b
    with open(trimmedfasta) as f:
        with open(shuffledfile, 'w') as f1:
            for line in f:
                temp = list(line[0:len(line)- 1])
                if dinucshuff:
                    dinucl_temp = dinucl(temp) #makes a list of nucleotide pairs
                    shuffled_kmers = shuffleACopy(dinucl_temp)    
                else:
                    shuffled_kmers = shuffleACopy(temp)
                for i in range(0, timestoshuffle):
                    shuffled_kmers = shuffleACopy(shuffled_kmers)
                shuffled_kmers.append('\n')
                f1.write(''.join(shuffled_kmers))
            f1.close()
            f.close()


# In[ ]:

def stringrandseq(string, dinucshuff, timestoshuffle):
    def dinucl(lst):
        i = 0
        dlst = []
        while i< len(lst):
            dlst.append(lst[i] + lst[i + 1]) #combines element with element after it
            i += 2
        return dlst

    def shuffleACopy(x):
            b = x[:] # make a copy of the keys
            random.shuffle(b) # shuffle the copy
            return b
    temp = list(line[0:len(line)- 1])
    if dinucshuff:
        dinucl_temp = dinucl(temp) #makes a list of nucleotide pairs
        shuffled_kmers = shuffleACopy(dinucl_temp)    
    else:
        shuffled_kmers = shuffleACopy(temp)
    for i in range(0, timestoshuffle):
        shuffled_kmers = shuffleACopy(shuffled_kmers)
    shuffled_kmers_str = ''.join(shuffled_kmers)
    return shuffled_kmers_str


# In[17]:

#COMPILE DATA ON KMERS FROM RANDOM FILE
def filerandseq_basicanalysis(shuffledfile, randdict):
    with open(shuffledfile) as f:
        for line in f: #one line corresponds to one sequence at that chromosome location
            for i in range(0, len(line) - 5):
                randseq_kmer = line[i:i+5] #makes 5-mer out of next five letters
                if randseq_kmer not in randdict:
                    randdict[randseq_kmer] = 1
                else:
                    randdict[randseq_kmer] += 1
        f.close()
    print ("Finished!")
    randseq_sortedlst = fsortbycount(randdict) #sorts rand kmers by frequency
    return randseq_sortedlst
def randseq_motifanalysis(randsortedlst, kdict, motif, kmerlength):
    randseq_motif_sorted = motif_analysis(randsortedlst, kdict, motif, kmerlength)
    return randseq_motif_sorted


# In[ ]:

#COMPILE DATA ON KMERS FROM RANDOM FILE
def stringrandseq_basicanalysis(randdict, string, length):
    for i in range(0, len(string) - length):
        randseq_kmer = line[i:i+length] #makes 5-mer out of next five letters
        if randseq_kmer not in randdict:
            randdict[randseq_kmer] = 1
        else:
            randdict[randseq_kmer] += 1
    print ("Finished!")
    randseq_sortedlst = fsortbycount(randseq_kmerdict) #sorts rand kmers by frequency
    return randseq_sortedlst

