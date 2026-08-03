"""
    This file is for the Entropy-Based Vectorization of DNA Sequences  
     
    An entropy vector is a list of entropy measures representing the sequence. 
    1000 long samples.
    100 long windows ( here represented as window_length )
    1 entropy measure per window  -> 10 entropy measures per sample
    
    tested in the sampling_seq.ipnyb file
      
"""
from curses import window


import sampling
import math
import entropy_measures


def get_kmers(sequence, k):
    """
        Extracts the k-mers of length k from the sequence.
        Returns a list of k-mers.
    
    """
    kmers = [] 
    for i in range(len(sequence)-k+1):
        kmers.append(sequence[i:i+k])
    return kmers # Return the list


class entropy_vector:
    def __init__(self, sample: sampling.mutated_sample, entropies: list ):
        self.sample = sample
        self.entropies = []
    
    def vectorize (self, window_length, k):
        """
        parameters:
            window_length: length of the window to chunk the sample into
            k: length of the k-mers to extract from each window
        Returns:  
            vector of entropy measures for the sample. -> one entropy measure per window
            vector will be found in the entropies attribute of the class instance.
            
        Note the windows are non-overlapping 
        
        """
        #Chunks the sample into wind of length window_length and computes the entropy for each window.
        
        for i in range(0,1000, window_length):
            wind_i = self.sample.sequence[i:i+window_length] 
            if len(wind_i) < window_length:
                break
            # Exttracts the kmers of length k for the window --> [list of k-mers] 
            list_kmers = get_kmers(wind_i, k)
            # Calls the entropy_measures  --- calculate_entropy(kmers:list) ----  to compute the entropy 
            # calculate_entropy returns a float value representing the entropy of the list of k-mers
            self.entropies.append(entropy_measures.calculate_entropy(list_kmers)) 
        
        
    
        
        
