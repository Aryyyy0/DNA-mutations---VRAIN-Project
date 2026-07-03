import os
import subprocess
from pathlib import Path
import pyranges
from Bio import SeqIO

def createFASTA(vcf,chromosome,out:str,sample=None):
    if sample == None:
        subprocess.run('bcftools consensus -f ' + chromosome +' '+ vcf + ' -o '+ out+'_all.fasta',shell=True)
        return None
    subprocess.run('bcftools consensus -f ' + chromosome +' '+ vcf + ' -o '+ out+'_'+sample+'.fasta' + ' -s '+sample,shell=True)
    return None

class Fasta():
    def __init__(self, fasta):
        self.record = SeqIO.read(fasta, "fasta")
        self.dna = self.record.seq
        self.sample = self.record.description
    def getKmers(self,annot,chr,k,size_window):
        print('Producing the kmers')
        kmers = dict()
        annot = annot[annot['Chromosome'] == chr]
        
        for elem in range(len(annot.index)):
            add = size_window - (annot['End'].iloc[elem]-annot['Start'].iloc[elem]-1)%size_window #so that also the last window has the correct size
            if add==size_window:
                add = 0
            if 'gene_id' in annot:
                kmers[annot['gene_id'].iloc[elem]] = [str(self.dna[i:i+k]) for i in range(annot['Start'].iloc[elem],annot['End'].iloc[elem]+add-1) ]
            else:
                kmers[annot['Element'].iloc[elem]] = [str(self.dna[i:i+k]) for i in range(annot['Start'].iloc[elem],annot['End'].iloc[elem]+add-1) ]
                
        
        return kmers
