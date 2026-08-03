import os 
import subprocess
from pathlib import Path
import math
from collections import defaultdict

#from decimal import Decimal

def log2(n):
    
    if n == 0: #the math.log2() method would raise a domain error
        return 0
    else:
        return math.log2(n)
def calculate_entropy(kmers:list): #calculate entropy for each gene of each patient
    entropy = 0
    count = len(kmers)
    dict_freq = defaultdict(int) 
    for k in kmers:
        dict_freq[k] += 1
    for freq in dict_freq.keys():
        entropy -= (dict_freq[freq]/count)*log2(dict_freq[freq]/count)
    return entropy

def calculate_entropy(kmers:list): 
    """  Float version of the entropy calculation method. Used in the entropy_vector class to compute the entropy of a list of k-mers."""
    entropy = 0.0
    count = len(kmers)
    dict_freq = defaultdict(int) 
    
    for k in kmers:
        dict_freq[k] += 1
        
    for freq in dict_freq.keys():
        # FORCE FLOAT DIVISION HERE:
        probability = float(dict_freq[freq]) / count
        entropy -= probability * log2(probability)
        
    return entropy 

def calculate_entropy_gene(kmers_sample:dict): #calculate entropy for each gene of each patient
    #print('Calculating gene entropy')
    dict_entropy = dict()
    for elem in kmers_sample.keys():
        dict_entropy[elem] = calculate_entropy(kmers_sample[elem]) #normalization to compare results 
    return  dict_entropy

def calculate_entropy_windows(kmers_sample:dict,size_window): #same as gene entropy, but on smaller portions of the gene at a time
        # the idea is to plot the average increase/decrease in entropy of mut. samples compared to the reference
    #print('Calculating windows entropy')
    dict_entropy = dict()
    for elem in kmers_sample.keys():
        dict_entropy[elem] = list()
        l = int(len(kmers_sample[elem])/size_window)
        if l == 0:
            dict_entropy[elem].append(calculate_entropy(kmers_sample[elem]))
        else: 
            for i in range(l):
                dict_entropy[elem].append(calculate_entropy(kmers_sample[elem][i*size_window:(i+1)*size_window]))   
    return dict_entropy


def calculate_entropy_kmers(kmers_dict:dict): #calculate the entropy for each position (by collecting kmers in that position from the different samples)
    #be aware: it's a bit slow
    #print('Calculating kmers entropy')
    dict_entropy = dict()
    for elem in kmers_dict['ref'].keys():
        k_list = list()
        dict_entropy[elem] = list()
        for j in range(len(kmers_dict['ref'][elem])):
            k_list = [kmers_dict[i][elem][j] for i in kmers_dict.keys()]
            dict_entropy[elem].append(calculate_entropy(k_list))
            
    return dict_entropy  


def distances(kmers,sample_name, entropy_sample,entropy_ref):
    #print('Calculating distances')
    dict_sc_distance = dict() #segment compositional distance using Shannon divergence
    dict_euclidean = dict() #euclidean distance from prob. distributions
    dict_manhattan = dict() #manhattan distance from prob. distributions
    for elem in kmers['ref'].keys():
        dict_freq_union = dict()
        k_list = list(set(kmers['ref'][elem] + kmers[sample_name][elem]))
        dict_freq_union['ref'] = {i:0 for i in k_list}
        dict_freq_union[sample_name] = {i:0 for i in k_list}
        for k in kmers['ref'][elem]:
            dict_freq_union['ref'][k] += 1
        for q in kmers[sample_name][elem]:
            dict_freq_union[sample_name][q] +=1
        
        entropy = 0
        eucl = 0
        man = 0
        probs = dict()
        probs['ref'] = {k: dict_freq_union['ref'][k]/len(kmers['ref'][elem]) for k in k_list}
        probs[sample_name] = {q: dict_freq_union[sample_name][q]/len(kmers[sample_name][elem]) for q in k_list}


        for kmer in k_list:
            entropy -= ((probs['ref'][kmer]+probs[sample_name][kmer])/2)*log2((probs['ref'][kmer]+probs[sample_name][kmer])/2)
            eucl += (probs['ref'][kmer]-probs[sample_name][kmer])**2
            man += abs(probs['ref'][kmer]-probs[sample_name][kmer])

        sc_distance = entropy-entropy_sample[elem]/2-entropy_ref[elem]/2
        if abs(sc_distance) <= 0.0000000001: #dealing with float point operations error
            sc_distance = 0
        
        dict_sc_distance[elem] = math.sqrt(sc_distance)
        dict_euclidean[elem] = math.sqrt(eucl)
        dict_manhattan[elem] = man
    
    return (dict_sc_distance, dict_euclidean, dict_manhattan)

def windows_distances(kmers,sample_name,size_window,entropy_sample,entropy_ref):
    #print('Calculating windows distances')
    dict_sc_distance = dict() #segment compositional distance using Shannon divergence
    dict_euclidean = dict() #euclidean distance from prob. distributions
    dict_manhattan = dict() #manhattan distance from prob. distributions
    for elem in entropy_ref.keys():
        l = int(len(kmers['ref'][elem])/size_window)
        dict_sc_distance[elem] = list()
        dict_euclidean[elem] = list()
        dict_manhattan[elem] = list()
        
        for i in range(l):
            dict_freq_union = dict()
            k_list = list(set(kmers['ref'][elem][i*size_window:(i+1)*size_window] + kmers[sample_name][elem][i*size_window:(i+1)*size_window]))
            dict_freq_union['ref'] = {i:0 for i in k_list}
            dict_freq_union[sample_name] = {i:0 for i in k_list}
            for k in kmers['ref'][elem][i*size_window:(i+1)*size_window]:
                dict_freq_union['ref'][k] += 1
            for q in kmers[sample_name][elem][i*size_window:(i+1)*size_window]:
                dict_freq_union[sample_name][q] +=1
            
            entropy = 0
            eucl = 0
            man = 0
            probs = dict()
            probs['ref'] = {k: dict_freq_union['ref'][k]/len(kmers['ref'][elem][i*size_window:(i+1)*size_window]) for k in k_list}
            probs[sample_name] = {q: dict_freq_union[sample_name][q]/len(kmers[sample_name][elem][i*size_window:(i+1)*size_window]) for q in k_list}
            for kmer in k_list:
                entropy -= ((probs['ref'][kmer]+probs[sample_name][kmer])/2)*log2((probs['ref'][kmer]+probs[sample_name][kmer])/2)
                eucl += (probs['ref'][kmer]-probs[sample_name][kmer])**2
                man += abs(probs['ref'][kmer]-probs[sample_name][kmer])

            sc_distance = entropy-entropy_sample[elem][i]/2-entropy_ref[elem][i]/2
            if abs(sc_distance) <= 0.0000000001:
                sc_distance = 0
            
            dict_sc_distance[elem].append(math.sqrt(sc_distance))
            dict_euclidean[elem].append(math.sqrt(eucl))
            dict_manhattan[elem].append(man)

    return (dict_sc_distance,dict_euclidean,dict_manhattan)