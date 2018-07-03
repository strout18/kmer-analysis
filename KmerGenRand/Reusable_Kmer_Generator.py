#### IMPORT STATEMENTS ####
import re
import plotly.plotly as py
import plotly.graph_objs as go
import plotly
import random


#### TRIM FUNCTION ####
#params: untrimmedfasta (filename of original fasta), trimmedout (empty file where you want the trimmed sequences to be saved)

def trim(untrimmedfasta, trimmedout): #trims FASTA file saved as a .txt by removing chr data and making all nucleotides capitalized
    with open(untrimmedfasta) as f:
        with open(trimmedout, 'w') as f1:
            for line in f:
                if line[0] != '>':
                    f1.write(line.upper())
            f1.close()
            f.close()

#### COMPILE FREQUENCY DATA FROM FILE ####
#inputs: params (dictionary - likely empty unless you're compiling one dict from mult files), trimmedfasta (file containing sequences 
#w/o chr data), length (length of k-mers desired)

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


#### COMPILE FREQUENCY DATA FROM STRING ####
#params: kmerdict (dictionary - empty unless you're compiling one dict from multiple strings), string (input sequence), length (desired 
#k-mer length)

def freq_compile_from_string(kmerdict, string, length):
    for i in range(0, len(string) - length):
        kmer = string[i:i+length]
        if 'N' not in kmer.upper():
            if kmer not in kmerdict:
                kmerdict[kmer] = 1
            else:
                kmerdict[kmer] += 1

#### SORT FREQUENCY DICTIONARY ####
#params: kdict (dictionary already filled with k-mers and corresponding frequencies)

def fsortbycount(kdict): #compiles a sorted list (tuple) of all k-mers with highest count in form (k-mer, count)
    sortbycount = sorted(kdict.items(), key = lambda x: x[1], reverse = True)
    return sortbycount

#### LIST ALL K-MERS IN FREQUENCY DICTIONARY ####
#params: sortedlst (list of tuples obtained from fsortbycount)

def listall(sortedlst): #list all k-mers found
    for x in sortedlst:
        print (x[0], ':', str(x[1]))

#### LIST MOST FREQUENTLY OCCURRING K-MERS ####
#params: sortedlst (list of tuples obtained from fsortbycount), cutoff (number of kmers to display- e.g. top 25 means cutoff = 25), 
#write_file (bool of whether or not to write this list to a file), fmostfreq (optional string corresponding to a file to write to)

def listmostfreq(sortedlst, cutoff, write_file, fmostfreq = ''): 
    i = 0
    if write_file:
        with open(fmostfreq, 'w') as f:
            while i< cutoff:
                print (sortedlst[i][0], ':', str(sortedlst[i][1]))
                f.write(sortedlst[i][0])
                if i!= cutoff - 1:
                    f.write('\n')
                i += 1
            f.close()
    else:
        while i< cutoff:
                print (sortedlst[i][0], ':', str(sortedlst[i][1]))
                i += 1

#### COMPARE FREQUENCIES OF K-MERS FROM A DESIRED MOTIF TO ALL K-MER FREQUENCIES ####
#params: sortedlst (list of tuples obtained from fsortbycount), kdict (completed dict of frequencies for all k-mers), motif in regex list
#form (SEE MOTIFINPUTHELP.MD), k-merlength (length of k-mers in frequency dictionary)

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


#### USE PLOTLY TO CREATE GRAPH OF ALL K-MERS AND COMPARE WITH MOTIF K-MERS WITH LINES TRACED TO MOTIF LABELS ####
#YOU MUST HAVE A PLOTLY ACCOUNT TO USE THIS FUNCTION (https://plot.ly/)
#params: sortedall (list of tuples of all k-mers obtained from fsortbycount), sorted_motif (list of k-mers occurring in given motif, this 
#is the result of the function motif_analysis), kdict (dict of all k-mers and their counts), graphtitle (desired graph title), 
#url (desired filename of the graph on the plot.ly website), gmode (graph mode, 'lines', 'markers', or 'lines+markers')

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

	
    plotly.tools.set_credentials_file(username='_____', api_key='______') #INSERT YOUR USERNAME AND API KEY HERE

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


#### USE PLOTLY TO CREATE GRAPH OF ALL K-MERS AND COMPARE WITH MOTIF K-MERS #### 
#YOU MUST HAVE A PLOTLY ACCOUNT TO USE THIS FUNCTION (https://plot.ly/)
#params: see plotwithmotiflines

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


#### USE PLOTLY TO GRAPH K-MERS AND THEIR CORRESPONDING FREQUENCIES ####
#YOU MUST HAVE A PLOTLY ACCOUNT TO USE THIS FUNCTION (https://plot.ly/)
#params: see plotwithmotiflines

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


#### READ IN A FILE OF SEQUENCES, RANDOMLY SHUFFLE BY LINE AND WRITE THESE RANDOMIZED SEQUENCES TO A NEW FILE ####
#params: trimmedfasta (filename of sequences w/o chr data), shuffledfile (location of file to write shuffled sequences to), 
#dinucshuff (bool, if True maintain dinucleotide frequencies when shuffling), timestoshuffle (integer # of times to shuffle)

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


#### READ IN A SEQUENCE STRING, RANDOMLY SHUFFLE AND RETURN IT ####
#params: string (input sequence w/o chr data), see filerandseq

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

#### COMPILE DATA ON K-MERS FROM RANDOM FILE ####
#params: shuffledfile (filename of file containing shuffled sequences), randdict (empty dictionary unless compiling from mult files),
#kmerlength: desired number of nucleotides in k-mer

def filerandseq_basicanalysis(shuffledfile, randdict, kmerlength):
    with open(shuffledfile) as f:
        for line in f: #one line corresponds to one sequence at that chromosome location
            for i in range(0, len(line) - kmerlength):
                randseq_kmer = line[i:i+kmerlength] #makes k-mer out of next five letters
                if randseq_kmer not in randdict:
                    randdict[randseq_kmer] = 1
                else:
                    randdict[randseq_kmer] += 1
        f.close()
    print ("Finished!")
    randseq_sortedlst = fsortbycount(randdict) #sorts rand kmers by frequency
    return randseq_sortedlst

#### COMPARE FREQUENCIES OF K-MERS FROM A GIVEN MOTIF TO ALL K-MER FREQUENCIES ####
#params: see motif_analysis

def randseq_motifanalysis(randsortedlst, kdict, motif, kmerlength):
    randseq_motif_sorted = motif_analysis(randsortedlst, kdict, motif, kmerlength)
    return randseq_motif_sorted

#### COMPILE DATA ON K-MERS FROM A SHUFFLED STRING ####
#params: string ( string containing shuffled sequence), randdict (empty dict unless compiling one dict from mult strings), length (k-mer
#length)
def stringrandseq_basicanalysis(string, randdict, length):
    for i in range(0, len(string) - length):
        randseq_kmer = line[i:i+length] #makes 5-mer out of next five letters
        if randseq_kmer not in randdict:
            randdict[randseq_kmer] = 1
        else:
            randdict[randseq_kmer] += 1
    randseq_sortedlst = fsortbycount(randseq_kmerdict) #sorts rand kmers by frequency
    return randseq_sortedlst

