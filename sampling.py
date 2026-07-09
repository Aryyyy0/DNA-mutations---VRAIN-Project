import itertools
from pyranges.readers import read_gtf
from processing_fasta import *
from collections import defaultdict
import pandas as pd
import math
import probabilistic_automaton as pa
import random


class mutated_sample:
   def __init__(self, id: int, chromosome: str, location: int, sequence: list, automata: pa.pdfa):
        """
        the mutated gene sample has:
         - id 
         - chromosome
         - location of the gene on the chromosome
         - automata : the trained DPFA used to generate it. form :
                            self.states = states
                            self.transitions = transitions #in the form {q1 : {char : (q2,prob)}}
                            self.final = final
                            self.initial = initial 
         """
        self.id = id
        self.chromosome = chromosome
        self.location = location
        self.sequence = sequence
        self.automata = automata
        
   def generate_mutated_sample(self):
        """
        generates a mutated sample using:
        - the trained DPFA
        
        """
        current = self.automata.initial
        sequence = ''
        while current not in self.automata.final:
            next_symbol, next_state = choose_random_transition(self.automata.transitions[current])
            sequence += next_symbol
            #print(sequence)
            current = next_state
        self.sequence = sequence
        print("The generated sequence is : ", self.sequence)

def choose_random_transition(possible_transitions: dict):
    """
    automata.transitions[q] is the dictionary of possible transitions : {char : (q2,prob)}
    
    """
    r = random.random() # Tosses a coin between 0.0 and 1.0
    cumulative_prob = 0.0
    
    for symbol, (next_state, transition_prob) in possible_transitions.items():
        cumulative_prob += transition_prob
        if r < cumulative_prob:
            return symbol, next_state
    
        
    

        
    
    