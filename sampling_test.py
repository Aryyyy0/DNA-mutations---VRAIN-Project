
import pandas
from pathlib import Path
import random
from processing_fasta import *
from entropy_measures import *
from ml_ktss import *
import argparse
import itertools
import numpy as np
from pyranges.readers import read_gtf
import probabilistic_automaton as pa
import sampling 


""" This file is to test the generating model of mutated sequences.In the sampling.py file
    Implemented using the trained DPFA probabilistic_automaton.py """
    
#parameters used for DPFA or/and sequence generation
gtf_file_path = '/data/gencode.v49lift37.basic.annotation_protein_coding.gtf'
chromosome = 'chr22'
location = 10 #its a number because on the gtf dataframe we are going to have genes ordered


my_automata = pa.automata_builder(gtf_file_path, chromosome, gene_number = location) # creating the automata:

print("The automata has been created successfully, ")
print(my_automata.final)
print(my_automata.initial)
#print(my_automata.states)
print("now we are going to generate sequences using the automata")

#---------------------------------------------------------------Generating test 1 -------------------------------------------------------
test1 =sampling.mutated_sample(id = 1, chromosome = chromosome, location = location, sequence = [], automata = my_automata)
test1.generate_mutated_sample()
print("The generated sequence is : ", test1.sequence)
print("The acceptance status and score of the generated sequence is : ", my_automata.accepts(test1.sequence))
print("The length of the generated sequence is : ", len(test1.sequence))




