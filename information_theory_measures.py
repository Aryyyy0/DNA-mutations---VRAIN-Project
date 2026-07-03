import os 
import subprocess
from pathlib import Path
import math
from collections import defaultdict
#from decimal import Decimal

def log2(n):
    if n == 0:
        return 0
    else:
        return math.log2(n)

def calculate_entropy(kmers:list): #calculate entropy for each gene of each patient; this is already pretty fast
    entropy = 0
    count = len(kmers)
    dict_freq = defaultdict(int) 
    for k in kmers:
        dict_freq[k] += 1
    for freq in dict_freq.keys():
        entropy -= (dict_freq[freq]/count)*log2(dict_freq[freq]/count)
    return entropy

def calculate_entropy_gene(kmers_sample:dict): #calculate entropy for each gene of each patient; this is already pretty fast
    print('Calculating gene entropy')
    dict_entropy = dict()
    for gene in kmers_sample.keys():
        dict_entropy[gene] = calculate_entropy(kmers_sample[gene]) #normalization to compare results 
    return  dict_entropy

def calculate_entropy_windows(kmers_sample:dict,size_window): #same as gene entropy, but on smaller portions of the gene at a time
        # don't worry that the values of the last window are smaller, you will compare the same portion across different samples therefore they will be on the same scale
        # the idea is to plot the average increase/decrease in entropy of mut. samples compared to the reference
    print('Calculating windows entropy')
    dict_entropy = dict()
    for gene in kmers_sample.keys():
        dict_entropy[gene] = list()
        l = int(len(kmers_sample[gene])/size_window)
        if l == 0:
            dict_entropy[gene].append(calculate_entropy(kmers_sample[gene]))
        else: 
            for i in range(l):
                dict_entropy[gene].append(calculate_entropy(kmers_sample[gene][i*size_window:(i+1)*size_window]))   
    return dict_entropy


def calculate_entropy_kmers(kmers_dict:dict): #calculate the entropy for each position (by collecting kmers in that position from the different samples)
    #be aware: it's a bit slow
    print('Calculating kmers entropy')
    dict_entropy = dict()
    for gene in kmers_dict['ref'].keys():
        k_list = list()
        dict_entropy[gene] = list()
        for j in range(len(kmers_dict['ref'][gene])):
            k_list = [kmers_dict[i][gene][j] for i in kmers_dict.keys()]
            dict_entropy[gene].append(calculate_entropy(k_list))
            
    return dict_entropy  


def distances(kmers,sample_name, entropy_sample,entropy_ref):
    print('Calculating compositional distance')
    dict_sc_distance = dict() #segment compositional distance using Shannon divergence
    dict_euclidean = dict() #euclidean distance from prob. distributions
    dict_manhattan = dict() #manhattan distance from prob. distributions
    for gene in kmers['ref'].keys():
        dict_freq_union = dict()
        k_list = list(set(kmers['ref'][gene] + kmers[sample_name][gene]))
        dict_freq_union['ref'] = {i:0 for i in k_list}
        dict_freq_union[sample_name] = {i:0 for i in k_list}
        for k in kmers['ref'][gene]:
            dict_freq_union['ref'][k] += 1
        for q in kmers[sample_name][gene]:
            dict_freq_union[sample_name][q] +=1
        
        entropy = 0
        eucl = 0
        man = 0
        probs = dict()
        probs['ref'] = {k: dict_freq_union['ref'][k]/len(kmers['ref'][gene]) for k in k_list}
        probs[sample_name] = {q: dict_freq_union[sample_name][q]/len(kmers[sample_name][gene]) for q in k_list}


        for kmer in k_list:
            entropy -= ((probs['ref'][kmer]+probs[sample_name][kmer])/2)*log2((probs['ref'][kmer]+probs[sample_name][kmer])/2)
            eucl += (probs['ref'][kmer]-probs[sample_name][kmer])**2
            man += abs(probs['ref'][kmer]-probs[sample_name][kmer])

        sc_distance = entropy-entropy_sample[gene]/2-entropy_ref[gene]/2
        if abs(sc_distance) <= 0.0000000001:
            sc_distance = 0
        
        dict_sc_distance[gene] = math.sqrt(sc_distance)
        dict_euclidean[gene] = math.sqrt(eucl)
        dict_manhattan[gene] = man
    
    return (dict_sc_distance, dict_euclidean, dict_manhattan)

def windows_distances(kmers,sample_name,size_window,entropy_sample,entropy_ref):
    print('Calculating windows distances')
    dict_sc_distance = dict() #segment compositional distance using Shannon divergence
    dict_euclidean = dict() #euclidean distance from prob. distributions
    dict_manhattan = dict() #manhattan distance from prob. distributions
    for gene in entropy_ref.keys():
        l = int(len(kmers['ref'][gene])/size_window)
        dict_sc_distance[gene] = list()
        dict_euclidean[gene] = list()
        dict_manhattan[gene] = list()
        
        for i in range(l):
            dict_freq_union = dict()
            k_list = list(set(kmers['ref'][gene][i*size_window:(i+1)*size_window] + kmers[sample_name][gene][i*size_window:(i+1)*size_window]))
            dict_freq_union['ref'] = {i:0 for i in k_list}
            dict_freq_union[sample_name] = {i:0 for i in k_list}
            for k in kmers['ref'][gene][i*size_window:(i+1)*size_window]:
                dict_freq_union['ref'][k] += 1
            for q in kmers[sample_name][gene][i*size_window:(i+1)*size_window]:
                dict_freq_union[sample_name][q] +=1
            
            entropy = 0
            eucl = 0
            man = 0
            probs = dict()
            probs['ref'] = {k: dict_freq_union['ref'][k]/len(kmers['ref'][gene][i*size_window:(i+1)*size_window]) for k in k_list}
            probs[sample_name] = {q: dict_freq_union[sample_name][q]/len(kmers[sample_name][gene][i*size_window:(i+1)*size_window]) for q in k_list}
            for kmer in k_list:
                entropy -= ((probs['ref'][kmer]+probs[sample_name][kmer])/2)*log2((probs['ref'][kmer]+probs[sample_name][kmer])/2)
                eucl += (probs['ref'][kmer]-probs[sample_name][kmer])**2
                man += abs(probs['ref'][kmer]-probs[sample_name][kmer])

            sc_distance = entropy-entropy_sample[gene][i]/2-entropy_ref[gene][i]/2
            if abs(sc_distance) <= 0.0000000001:
                sc_distance = 0
            
            dict_sc_distance[gene].append(math.sqrt(sc_distance))
            dict_euclidean[gene].append(math.sqrt(eucl))
            dict_manhattan[gene].append(man)

    return (dict_sc_distance,dict_euclidean,dict_manhattan)