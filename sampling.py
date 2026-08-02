
from pyranges.readers import read_gtf
from processing_fasta import *
from collections import defaultdict
import pandas as pd
import probabilistic_automaton as pa
import random


class mutated_sample:
   def __init__(self, id: int, chromosome: str, location: int, sequence: list, automata: pa.pdfa, status: bool):
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
        self.status = True # default value
        
   def generate_mutated_sample(self):
        """
        Method filling the sequence field of the mutated sample instance.
        generates a mutated sample using the trained DPFA. No control on the length of the sample.
        
        """
        current = self.automata.initial
        sequence = ''
        while current not in self.automata.final:
            next_symbol, next_state = choose_random_transition(self.automata.transitions[current])
            sequence += next_symbol
            #print(sequence)
            current = next_state
        self.sequence = sequence
        #print("The generated sequence is : ", self.sequence)
        
   def generate_mutated_sample_length_control(self, target_length):
        """
        Generates a sequence of exactly 'target_length' using the trained DPFA.
        """
        counter = 0
        status = True
        current = self.automata.initial
        sequence = ''
        
        while counter < target_length and status:
            possible_transitions = self.automata.transitions.get(current, {}) 
            
            # SAFETY CHECK: If the automaton hits a dead-end (leaf) before 1000
            if not possible_transitions:
                status = False
                break
              
            # Use standard RandomQ extraction naturally
            next_symbol, next_state = choose_random_transition(possible_transitions)
            
            sequence += next_symbol
            current = next_state
            counter += 1
            
        self.sequence = sequence
        self.status = status

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
    
