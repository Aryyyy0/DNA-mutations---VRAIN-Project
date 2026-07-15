
import pandas as pd
import os
import probabilistic_automaton as pa
import sampling 


""" This file is to test the generating model of mutated sequences.In the sampling.py file
    Implemented using the trained DPFA probabilistic_automaton.py """
  
"""
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
"""
#---------------------------------------------------------------Generating test 1 -------------------------------------------------------
"""test1 =sampling.mutated_sample(id = 1, chromosome = chromosome, location = location, sequence = [], automata = my_automata)
test1.generate_mutated_sample()
print("The generated sequence is : ", test1.sequence)
print("The acceptance status and score of the generated sequence is : ", my_automata.accepts(test1.sequence))
print("The length of the generated sequence is : ", len(test1.sequence))
"""
#---------------------------------------------------------------Generating test 2 -------------------------------------------------------
""" creating N new samples and saving them in a dataframe with their acceptance status and score, and their length, and saving it as a csv file in the newly_generated_sequences folder 
N = 10
chromosome = 'chr21'
location = 20
pd.DataFrame.generated_sequences = []
for i in range(N):
    test =sampling.mutated_sample(id = i, chromosome = chromosome, location = location, sequence = [], automata = my_automata)
    test.generate_mutated_sample()
    acceptance_status, score = my_automata.accepts(test.sequence)
    pd.DataFrame.generated_sequences.append([test.id, len(test.sequence), acceptance_status, score])
    
print("The generated sequences are : ", pd.DataFrame.generated_sequences) """

# maybe orientate the generation in a folder called newly_generated_sequences
# see the scores of acceptance of the generated sequences, and see if they are accepted or not by the automata 
# create a data frame with the generated sequences their length and their acceptance status and score, and save it as a csv file in the newly_generated_sequences folder
#
#-----------------------Generating 100 sequences for each location on each chromosome and saving them in a database --------


database_rows = [] #CSV file 
chr_list = ['chr22'] # for the moment running only on chr22, but we can run it on all chromosomes later
gtf_file_path = '/data/gencode.v49lift37.basic.annotation_protein_coding.gtf' # to build on the automata
for chr in chr_list:
    os.makedirs(f"generated_sequences/{chr}", exist_ok=True) #folder for each chromosome
    
    for location in range(1, 10): # admitting that there are 10 diff locations 
        my_automata = pa.automata_builder(gtf_file_path, chr, gene_number = location) # creating the automata
        N = 100 # number of samples to generate
        for i in range(N):
            sample = sampling.mutated_sample(id = i, chromosome = chr, location = location, sequence = [], automata = my_automata)
            sample.generate_mutated_sample()
            acceptance_status, score = my_automata.accepts(sample.sequence)

            # create a text file for the generated sequence
            file_name = f"{location}_sample_{i}.fasta"
            file_path = f"generated_sequences/{chr}/{file_name}"
            # save the generated sequence in a fasta file
            with open(file_path, "w") as f:
                f.write(f">{location}_sample_{i}\n") # Standard FASTA header
                f.write(sample.sequence)
            
            # save the metadata AND the file path to our database list
            database_rows.append({
                "sample_id": i,
                "chromosome": chr,
                "gene_location": location,
                #"start_loc": start_loc,
                #"end_loc": end_loc,  
                "sequence_length": len(sample.sequence),
                "acceptance_status": acceptance_status,
                "acceptance_score": score,
                "file_path": file_path
            })

            # export the clean, readable database to a CSV
            df = pd.DataFrame(database_rows)
            df.to_csv("mutated_samples_database.csv", index=False)
            
my_database = pd.read_csv("mutated_samples_database.csv")
print(my_database.head()) # visualize the first lignes of the database 

