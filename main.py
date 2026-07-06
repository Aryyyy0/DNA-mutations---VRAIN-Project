import os 
import subprocess
import pandas
from pathlib import Path
import random
from processing_fasta import *
from entropy_measures import *
#from plots import *
from ml_ktss import *
import argparse
import itertools
import numpy as np
from pyranges.readers import read_gtf

# The Command Line Interface:
""" Example of usage:
    To run the Entropy analysis:

    python main.py --vcf my_snps.vcf --annot my_genes.txt entropy --kmer 5 --samples patient1,patient2

    To run the KTSS analysis:

    python main.py --vcf my_snps.vcf --annot my_genes.txt ktss --chromosome chr21
"""
# It is Mandatory to specify the VCF file and the annotation file. The user can choose to calculate entropy measures or to perform machine learning using automata.
parser = argparse.ArgumentParser(prog='main.py')
parser.add_argument('--vcf',help='Name of the VCF file',type=str,required=True)
parser.add_argument('--annot','-a',help='Name of the annotation file',type=str,required=True)
# Entropy measures
subparser = parser.add_subparsers(dest='command')
entropy = subparser.add_parser('entropy')
entropy.add_argument('--kmer','-k',type=int,help='k-mer size',default=7)
entropy.add_argument('--window','-w',type=int,help='Window size',default=100)
entropy.add_argument('--stats','-s',type=bool,help='Element statistics plot', default=False)
entropy.add_argument('--chromosome','-c',type=str,help='Chromosome name, in the form chrXX')
entropy.add_argument('--filtered','-f',type=bool,help='Whether the VCF file is filtered (synonymous variants out)',default=False)
entropy.add_argument('--samples','-sm',help='List of sample ids',required=True)
# Automata
ktss = subparser.add_parser('ktss')
ktss.add_argument('--kmer','-k',type=int,help='k-mer size',default=10)
ktss.add_argument('--chromosome','-c', type=str,help='Chromosome name, in the form chrXX')
# Reading all the arguments from the command line 
args = parser.parse_args()


#the code supposes the presence of a folder 'data' where FASTAs, VCFs and the annotation are present
cwd = Path(os.getcwd())
#considering that the annotation file can either be a GTF or a BED file
if "gtf" in args.annot:
    annot = read_gtf(str(cwd)+'/data/'+args.annot,as_df = True)
else:
    annot = pandas.read_table(str(cwd)+'/data/'+args.annot,columns=['Chromosome','Start','End','ID','Element','Type'])


with open(str(cwd)+'/data/'+args.samples,'r') as file: #need to write it 
    s = file.read()
sample_list = s.split()

#box_plot(annot,chrom_list,'stats/genes_pc_boxplot.pdf')

#ENTROPY MEASURES
if args.command == 'entropy':
    if args.chromosome:
        chrom_list = list(cwd.glob('data/'+args.chromosome+'.fasta')) 
    else:
        chrom_list = list(cwd.glob('data/*.fasta'))

    #ITERATING INSIDE THE CHROMOSOMES
    for i in range(len(chrom_list)): ## to save up space, everything will be calculated one chromosome at a time
                                    ## after the plots have been produced, data will be deleted
        chromosome = str(chrom_list[i].stem)
        if os.path.isdir('plots/'+chromosome) == False:
            subprocess.run('mkdir plots/'+chromosome,shell=True)
        if os.path.isdir('results/'+chromosome) == False:
            subprocess.run('mkdir results/'+chromosome,shell=True)
        if os.path.isdir('mutated')==False:
            subprocess.run('mkdir mutated',shell=True)

        plots = 'plots/'+chromosome
        res = 'results/'+chromosome

        if args.filtered == True:
            subprocess.run('mkdir plots/'+chromosome+'/filtered',shell=True)
            subprocess.run('mkdir results/'+chromosome+'/filtered',shell=True)
            plots = 'plots/'+chromosome+'/filtered'
            res = 'results/'+chromosome+'/filtered'
        
        #genes_stats(annot[annot['Chromosome']==chromosome],chromosome,'stats/stats_'+chromosome)
        [createFASTA(str(cwd)+'/data/'+args.vcf,str(chrom_list[i]),'mutated/'+chromosome,sample=sample_id) for sample_id in sample_list]
        print('Working on '+chromosome)
        
        k = args.kmer
        window = args.window
        ref = Fasta(str(chrom_list[i]))
        
        #CALCULATING THE ENTROPY IN EACH SAMPLE SEQUENCE
        mut = list(cwd.glob('mutated/*.fasta')) # product of the createFASTA call 
        print('k = '+str(k))

        entropy_gene = dict()
        seg_comp_distance = dict()
        eucl_distance = dict()
        manh_distance = dict()
        entropy_gene_windows = dict()
        w_seg_comp_distance = dict()
        w_eucl_distance = dict()
        w_manh_distance = dict()
        

        dict_ref = ref.getKmers(annot,chromosome,k,window)
        lengths = {gene : len(dict_ref[gene]) for gene in dict_ref.keys()}
        entropy_windows_ref = calculate_entropy_windows(dict_ref,window)
        entropy_ref = calculate_entropy_gene(dict_ref)
        
        for j in range(len(mut)): #iterate over the samples to calculate kmers and entropy
            seq = Fasta(str(mut[j]))
            sample = mut[j].stem  
            dict_kmers = dict()
            
            dict_kmers = seq.getKmers(annot,chromosome,k,window)
            entropy_gene[sample] = calculate_entropy_gene(dict_kmers)
            entropy_gene_windows[sample] = calculate_entropy_windows(dict_kmers,window)
            (seg_comp_distance[sample],eucl_distance[sample],manh_distance[sample]) = distances({'ref':dict_ref,sample:dict_kmers},sample,entropy_gene[sample],entropy_ref)
            (w_seg_comp_distance[sample],w_eucl_distance[sample],w_manh_distance[sample]) = windows_distances({'ref':dict_ref,sample:dict_kmers},sample,window,entropy_gene_windows[sample],entropy_windows_ref)
                
        #ENTROPY DIFFERENCE BETWEEN THE FASTA WITH ALL MUTATIONS AND THE REFERENCE
        '''createFASTA(str(cwd)+'/data/'+args.vcf,str(chrom_list[i]),'mutated/'+chromosome) #creates fasta where all mutations from all patients are applied
        dict_kmers['all'] = Fasta(str(cwd.glob('mutated/*all.fasta'))).getKmers(annot,chromosome,k,window)
        entropy_all = {'all': calculate_entropy_gene(dict_kmers['all'])}
        entropy_windows_all = {'all': calculate_entropy_windows(dict_kmers['all'])}
        diff = pandas.DataFrame.from_dict(entropy_gene_plot(annot[annot['Chromosome']==chromosome],entropy_all,entropy_ref,chromosome,''plots/'+chromosome+'/gene_entropy_all_'+chromosome+'_k'+str(k)+'.pdf'),orient='index')
        diff.to_csv(''results/'+chromosome+'/avg_gene_all_'+chromosome+'_k'+k+'.txt', sep='\t', index=True, header=False)
        '''
        dict_kmers = dict() #freeing some RAM for the plots

        #AVERAGE ENTROPY DIFFERENCE BETWEEN EACH SAMPLE AND THE REFERENCE
        pandas.DataFrame.from_dict(entropy_gene_plot(annot[annot['Chromosome']==chromosome],lengths,entropy_gene,entropy_ref,chromosome,plots+'/norm_pc_gene_entropy_diff_'+chromosome+'_k'+str(k)+'.pdf'),orient='index').to_csv(res+'/avg_pc_gene_'+chromosome+'_k'+str(k)+'.txt', sep='\t', index=True, header=False)
        pandas.DataFrame.from_dict(entropy_windows_plot(annot[annot['Chromosome']==chromosome],entropy_gene_windows,entropy_windows_ref,chromosome,plots+'/window_entropy_diff_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.pdf'),orient='index').to_csv(res+'/avg_window_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.txt', sep='\t', index=True, header=False)
        
        #ENTROPY DISTANCES
        pandas.DataFrame.from_dict(distance_plot(annot[annot['Chromosome']==chromosome],seg_comp_distance,lengths,chromosome,plots+'/norm_pc_seg_comp_distance_'+chromosome+'_k'+str(k)+'.pdf','Normalized segmental compositional distance'),orient='index').to_csv(res+'/avg_pc_scd_'+chromosome+'_k'+str(k)+'.txt',sep='\t',index=True,header=False)
        pandas.DataFrame.from_dict(distance_plot(annot[annot['Chromosome']==chromosome],eucl_distance,lengths,chromosome,plots+'/norm_pc_eucl_distance_'+chromosome+'_k'+str(k)+'.pdf','Normalized Euclidean distance'),orient='index').to_csv(res+'/avg_pc_eucl_'+chromosome+'_k'+str(k)+'.txt',sep='\t',index=True,header=False)
        pandas.DataFrame.from_dict(distance_plot(annot[annot['Chromosome']==chromosome],manh_distance,lengths,chromosome,plots+'/norm_pc_manh_distance_'+chromosome+'_k'+str(k)+'.pdf','Normalized Manhattan distance'),orient='index').to_csv(res+'/avg_pc_man_'+chromosome+'_k'+str(k)+'.txt',sep='\t',index=True,header=False)
        
        #WINDOWS ENTROPY DISTANCES 
        pandas.DataFrame.from_dict(windows_distance_plot(annot[annot['Chromosome']==chromosome],w_seg_comp_distance,chromosome,plots+'/seg_comp_distance_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.pdf','Windows segmental compositional distance'),orient='index').to_csv(res+'/avg_scd_w_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.txt',sep='\t',index=True,header=False)
        pandas.DataFrame.from_dict(windows_distance_plot(annot[annot['Chromosome']==chromosome],w_eucl_distance,chromosome,plots+'/eucl_distance_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.pdf','Windows Euclidean distance'),orient='index').to_csv(res+'/avg_eucl_w_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.txt',sep='\t',index=True,header=False)
        pandas.DataFrame.from_dict(windows_distance_plot(annot[annot['Chromosome']==chromosome],w_manh_distance,chromosome,plots+'/manh_distance_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.pdf','Windows Manhattan distance'),orient='index').to_csv(res+'/avg_man_w_'+chromosome+'_k'+str(k)+'_w'+str(window)+'.txt',sep='\t',index=True,header=False)
        
        #ENTROPY VALUES AT EACH POSITION OF EACH GENE
        '''del dict_kmers['ref']
        entropy_kmers = calculate_entropy_kmers(dict_kmers)  
        entropy_kmers_plot(annot[annot['Chromosome']==chromosome],entropy_kmers,chromosome,''plots/'+chromosome+'/kmers_entropy_'+chromosome+'_k'+str(k)+'.pdf') 
        pandas.DataFrame.from_dict(entropy_kmers).to_csv('results/'+chromosome+'/kmers_entropy_'+chromosome+'_k'+str(k)+'.txt')
        '''

        if args.stats: #statistics when calculating gene entropy, useful due to the gene lenght bias issue
            distr_variants_genes(annot,str(cwd)+'/data/'+args.vcf,pandas.read_table('results/'+chromosome+'/avg_pc_gene_'+chromosome+'_k'+str(k)+'.txt',sep='\t',header=None,names=['gene','value']),'stats/variant_distribution_'+chromosome+'_k'+str(k)+'.pdf',50,chromosome)
        
        #REMOVING MUTATED FASTA FILES (they take too much space)
        for file in cwd.glob('mutated/*'):
            os.remove(file)
        

#MACHINE LEARNING USING AUTOMATA 
if args.command == 'ktss':
    k = args.kmer
    if args.chromosome:
        chrom_list = list(cwd.glob('data/'+args.chromosome+'.fasta')) 
    else:
        chrom_list = list(cwd.glob('data/*.fasta'))

    all_k = all_kmers(itertools.product('ACGT',repeat=k))

    for i in range(len(chrom_list)): 
        chromosome = str(chrom_list[i].stem)
        #[createFASTA(str(cwd)+'/data/'+args.vcf,str(chrom_list[i]),'mutated/'+chromosome,sample=sample_id) for sample_id in sample_list]
        mut = list(cwd.glob('mutated/*.fasta'))
        annot_c = annot[annot['Chromosome']==chromosome]
        error = []
       
        for gene in range(20):
            p_set = None
            p_set = [Fasta(str(mut[j])).dna[annot_c['Start'].iloc[gene]:annot_c['End'].iloc[gene]] for j in range(len(mut))] #get DNA sequences from gene coordinates
            #print(p_set)
            n_set = None
            print(annot_c['gene_id'].iloc[gene])
            error_p = []
            #error_n = []

            for n in range(len(mut)): #leave-one-out validation
                p_automata = None
                p_automata = ktss_language(p_set[0:n]+p_set[n+1:len(p_set)],k,all_k).DFA_infer() #automata for patients
                #n_automata = ktss_language([n_set[:n],n_set[n+1:]],k).DFA_infer() # for non-patients individuals
                error_p.append(1-np.average(p_automata.accepts_input(p_set[n])))#, not(n_automata.accepts_input(p_set[i])))
                #error_n.append(1-np.average(not(p_automata.accepts_input(n_set[n])),n_automata.accepts_input(n_set[i])))
            
            error.append(np.average(error_p))#,error_n])
            
        pandas.DataFrame(error,index=annot_c['gene_id'][:20]).to_csv('results/'+chromosome+'/error_ktss'+chromosome+'_k'+str(k)+'.txt',sep='\t')

print('Done')