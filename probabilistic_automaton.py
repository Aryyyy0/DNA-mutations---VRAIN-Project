import itertools
from pyranges.readers import read_gtf
from processing_fasta import *
from collections import defaultdict
import pandas as pd
import math

"""
This file contains the implementation of a probabilistic deterministic finite automaton (PDFA) and 
its application to parse sequences of nucleotides. The PDFA is defined by its states, transitions, final states, and initial state. 
The transitions are represented as a dictionary where each state maps to another dictionary of characters and their corresponding next state and transition probability.
"""


class pdfa():
    def __init__(self,states,transitions,final,initial):
        self.states = states
        self.transitions = transitions #in the form {q1 : {char : (q2,prob)}}
        self.final = final
        self.initial = initial 
        
  
    def accepts(self,seq):
        """ 
        function parsing a sequence through the automata, 
                 if the sequence is accepted, return True and the score of acceptance, 
                 else return False and the score of rejection (sum of log(probabilities of the transitions))
                 
        """
        current = self.transitions[self.initial][seq[0]][0] # the state we transition to from the initial state with the first character of the sequence
        prob = self.transitions[self.initial][seq[0]][1] # the probability of transitioning from the initial state to the current state with the first character of the sequence
        for i in range(1,len(seq)):
            if seq[i] not in self.transitions[current].keys():  #rejection with -inf probability
                return (False,float('-inf'))

            # update the probability of acceptance by summing on logarithms of the probabilities of the transitions
            prob = prob + math.log(self.transitions[current][seq[i]][1] )
            current = self.transitions[current][seq[i]][0]
            
        # if we finish parsing the sequence and we are in a final state, return True and the probability of acceptance
        if current in self.final:
            return (True,prob)

        return (False,prob)

def all_kmers(l):
    return set(''.join(t) for t in l)

def get_kmers(seqs,k):
    kmers = set()
    for i in range(len(seqs)):
        kmers.update([seqs[i][j:j+k] for j in range(len(seqs[i])-k+1)])
    return set(kmers)

"""
       - Building 20 automata for 20 genes on the 22nd chromosome, 
         and testing them on the corresponding sequences of the first fasta file in the mutated folder.
       - Output is a dataframe with the gene ids as index and the acceptance status and score for each gene.
       - the closest the score is to 0, the more likely the sequence is to be accepted by the automata.
 
"""

# '/data/gencode.v49lift37.basic.annotation_protein_coding.gtf'

def automata_builder(gtf_file_path: str, chromosome: str, gene_number: int):
    """
        parameters:
            gtf_file_path: path to the gtf file containing the gene annotations
            chromosome: chromosome number to filter the gtf file
            gene_number: index of the gene to build the automata for    
        returns:
            automata: a probabilistic deterministic finite automaton (PDFA) built from the sequences of the specified gene
    """
    cwd = Path(os.getcwd())
    gtf = read_gtf(str(cwd)+gtf_file_path,as_df = True)
    gtf = gtf[gtf['Chromosome']==chromosome]
    mut = list(cwd.glob('mutated/*.fasta'))
    
    R = [Fasta(str(mut[j])).dna[gtf['Start'].iloc[gene_number]:gtf['End'].iloc[gene_number]] for j in range(1,len(mut))] # Training Data Preparation 
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
    
    return automata

def automata_tester(automata,gtf_file_path: str , chromosome: str, gene_number: int):

    cwd = Path(os.getcwd())
    gtf = read_gtf(str(cwd)+gtf_file_path,as_df = True)
    gtf = gtf[gtf['Chromosome']==chromosome]
    mut = list(cwd.glob('mutated/*.fasta'))
    
    test = Fasta(str(mut[0])).dna[gtf['Start'].iloc[gene_number]:gtf['End'].iloc[gene_number]]
    
    return automata.accepts(test)

#-------------------------------------------- TESTING THE AUTOMATA BUILDER AND TESTER FUNCTIONS --------------------------------------------
if __name__ == '__main__':
    gtf_file_path = '/data/gencode.v49lift37.basic.annotation_protein_coding.gtf'
    output = []
    for gene_i in range(20):
        automata_i = automata_builder(gtf_file_path,'chr22',gene_i)
        test_i = automata_tester(automata_i,gtf_file_path,'chr22', gene_i)
        output.append(test_i)

    # Pandas DataFrame to display the results :  gene_id |acceptance status | acceptance score
    cwd = Path(os.getcwd())
    gtf = read_gtf(str(cwd)+'/data/gencode.v49lift37.basic.annotation_protein_coding.gtf',as_df = True)
    gtf = gtf[gtf['Chromosome']=='chr22']
    print(pd.DataFrame(output,index=gtf['gene_id'][:20]))