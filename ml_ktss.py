import math 
import numpy as np 
from automata.fa.dfa import DFA 
from collections import defaultdict
import itertools

def all_kmers(l):
    return set(''.join(t) for t in l)

def get_kmers(seqs,k):
    kmers = set()
    for i in range(len(seqs)):
        kmers.update([seqs[i][j:j+k] for j in range(len(seqs[i])-k+1)])
    return set(kmers)


class ktss_language():
    def __init__(self,R,k,all_k): #R set of sequences, therefore set of DNA strings, one for each individual, spanning the gene
        self.sigma = ['A','C','G','T',""]
        self.k = k
        self.I = set([R[i][:k-1] for i in range(len(R))]) #initial kmers
        self.F = set([R[i][len(R[i])-k+1:len(R[i])] for i in range(len(R))]) #terminal kmers
        
        
        self.all = all_k
        self.T = set(self.all.difference(get_kmers(R,k)))  #'prohibited' kmers
        
    def DFA_infer(self):
        k = self.k
        sigma = self.sigma
        I = self.I 
        F = self.F 
        T = self.T 
        Q = {""}
        q_0 =  "" #for our purpose, the initial state is the empty string
        delta = dict()
        for a in I:
            a = str(a)
            for j in range(len(a)):
                if j == 0:
                    Q.add(a[0])
                    if "" in delta.keys():
                        delta[""][a[0]] = a[0]
                    else:
                        delta[""] = {a[0]:a[0]}
                else:
                    Q.add(a[0:j+1])
                    if a[0:j] in delta.keys():
                        delta[a[0:j]][a[j]] = a[0:j+1]
                    else:
                        delta[a[0:j]] = {a[j]:a[0:j+1]}
        #print('ok')
        for x in self.all.difference(T):
            x = str(x)
            Q.add(x[1:k])
            if x[0:k-1] in delta.keys():
                delta[x[0:k-1]][x[k-1]] = x[1:k]
            else:
                delta[x[0:k-1]] = {x[k-1]:x[1:k]}

        for q in Q:
            if q not in delta.keys():
                delta[q] = {"":q}

        #print('ok')
        Q_f = F
        return DFA(states=set(Q),input_symbols=sigma,transitions=delta,initial_state=q_0,final_states=Q_f,allow_partial=True)