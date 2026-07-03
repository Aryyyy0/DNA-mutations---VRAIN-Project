##PROBABILISTIC AUTOMATA FOR GENOMIC ANALYSIS, this is just a test
import probabilistic_automata as pa 
import itertools
from pyranges.readers import read_gtf
from processing_fasta import *
from collections import defaultdict
import pandas

class pdfa():
    def __init__(self,states,transitions,final,initial):
        self.states = states
        self.transitions = transitions #in the form {q1 : {char : (q2,prob)}}
        self.final = final
        self.initial = initial 

    def accepts(self,seq):
        current = self.transitions[self.initial][seq[0]][0]
        prob = self.transitions[self.initial][seq[0]][1]
        for i in range(1,len(seq)):
            if seq[i] not in self.transitions[current].keys():
                return (False,0)

            prob = prob * self.transitions[current][seq[i]][1]
            current = self.transitions[current][seq[i]][0]
            
        
        if current in self.final:
            if prob >= 0.8:
                return (True,prob)
            return (False,prob)

        return (False,prob)

def all_kmers(l):
    return set(''.join(t) for t in l)

def get_kmers(seqs,k):
    kmers = set()
    for i in range(len(seqs)):
        kmers.update([seqs[i][j:j+k] for j in range(len(seqs[i])-k+1)])
    return set(kmers)

cwd = Path(os.getcwd())
gtf = read_gtf(str(cwd)+'/data/gencode.v49lift37.basic.annotation_protein_coding.gtf',as_df = True)
chromosome = 'chr22'
gtf = gtf[gtf['Chromosome']==chromosome]
mut = list(cwd.glob('mutated/*.fasta'))

output = []

for gene in range(20):

    R = [Fasta(str(mut[j])).dna[gtf['Start'].iloc[gene]:gtf['End'].iloc[gene]] for j in range(1,len(mut))]
    k = 10
    all_k = all_kmers(itertools.product('ACGT',repeat=k))
    #sigma = ['A','C','T','G']
    I = set([R[i][:k-1] for i in range(len(R))]) #initial kmers
    F = set([R[i][len(R[i])-k+1:len(R[i])] for i in range(len(R))]) #terminal kmers
    T = set(all_k.difference(get_kmers(R,k)))
    start =  "" #for our purpose, the initial state is the empty string
    delta = dict()
    prob = defaultdict(lambda : defaultdict(int))
    Q = {""}
    
    for a in I:
        a = str(a)
        for j in range(len(a)):
            if j == 0:
                Q.add(a[0])
                if "" in delta.keys():
                    delta[""][a[0]] = (a[0],)
                    prob[""][a[0]] +=1
                else:
                    delta[""] = {a[0]:(a[0],)}
                    prob[""] = {a[0]:1}
            else:
                Q.add(a[0:j+1])
                if a[0:j] in delta.keys():
                    delta[a[0:j]][a[j]] = (a[0:j+1],)
                    prob[a[0:j]][a[j]] +=1
                else:
                    delta[a[0:j]] = {a[j]:(a[0:j+1],)}
                    prob[a[0:j]] = {a[j]:1}


    for x in all_k.difference(T):
        x = str(x)
        Q.add(x[1:k])
        if x[0:k-1] in delta.keys():
            delta[x[0:k-1]][x[k-1]] = (x[1:k],)
            if x[k-1] not in prob[x[0:k-1]].keys():
                prob[x[0:k-1]][x[k-1]] =1
            else: 
                prob[x[0:k-1]][x[k-1]] +=1
        else:
            delta[x[0:k-1]] = {x[k-1]:(x[1:k],)}
            prob[x[0:k-1]] = {x[k-1]:1}

    Q_f = F

    for i in prob.keys():
        count = 0
        for j in prob[i].keys():
            count += prob[i][j]
        for r in delta[i].keys():
            delta[i][r] = delta[i][r] + (prob[i][r]/count,)

    #print(delta[list(delta.keys())[1]])
    automata = pdfa(Q,delta,Q_f,start)

    test = Fasta(str(mut[0])).dna[gtf['Start'].iloc[gene]:gtf['End'].iloc[gene]]
    output.append(automata.accepts(test))

print(pandas.DataFrame(output,index=gtf['gene_id'][:20]))