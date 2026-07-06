import matplotlib.pyplot as plt
import numpy as np 
import math
import vcfpy
import pandas as pd

def genes_stats(annot,chr,file):
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 7)
    lengths = [annot['End'].iloc[i]-annot['Start'].iloc[i] for i in range(len(annot.index))]
    counts, bins = np.histogram(lengths,bins=1000)
    plt.stairs(counts, bins)
    ax.set_ylabel('Count')
    ax.set_xlabel('Length')
    ax.set_title('Gene lenght distribution on '+chr)
    plt.savefig(file+'.pdf')
    fig.clf(True)
    counts, bins = np.histogram(lengths,bins=1000,range=(0,25000))
    plt.stairs(counts, bins)
    ax.set_ylabel('Count')
    ax.set_xlabel('Length')
    ax.set_title('Gene lenght distribution on '+chr+' range (0,200000)')
    plt.savefig(file+'_magn.pdf')
    fig.clf(True)
    plt.close('all')
    return None

def box_plot(annot,chrs,file):
    lengths = list()
    labels = list()
    for i in range(len(chrs)):
        chr = str(chrs[i].stem)
        labels.append(chr)
        lengths.append([math.log10(annot[annot['Chromosome']==chr]['End'].iloc[i]-annot[annot['Chromosome']==chr]['Start'].iloc[i]) for i in range(len(annot[annot['Chromosome']==chr].index))])
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 7)
    plt.boxplot(lengths)
    ax.set_ylabel('Log10 gene size')
    ax.set_xlabel('Chromosome')
    plt.savefig(file)
    fig.clf(True)
    plt.close('all')
    return None


def entropy_plot(annot,entropy:dict,label,file):
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 7)
    ax.bar(annot['Start'],list(entropy.values()),width=1)
    ax.set_ylabel('Normalized Entropy difference')
    ax.set_xlabel('Position')
    ax.set_title(label)
    plt.savefig(file)
    fig.clf(True)
    plt.close('all')
    return None

def entropy_gene_plot(annot,l,entropy_samples,entropy_ref,chr,file):
    #print('Producing gene plot')
    avg_diff = dict()
    for gene in entropy_ref.keys():
        avg_diff[gene] = np.average([entropy_samples[i][gene]*math.log2(l[gene])-entropy_ref[gene]*math.log2(l[gene]) for i in entropy_samples.keys()])
    entropy_plot(annot, avg_diff,'Average gene normalized entropy difference on '+chr, file)  
    return avg_diff

def entropy_windows_plot(annot,entropy_samples,entropy_ref,chr,file):
    #print('Producing windows plot')
    fig,ax = plt.subplots()
    avg_diff = dict()
    plt.figure(figsize=(12,7))
    for gene in entropy_ref.keys():
        start = int(annot.loc[annot['gene_id']==gene,'Start'].iloc[0])
        #print(start)
        avg_diff[gene] = np.average([np.subtract(entropy_samples[i][gene],entropy_ref[gene]) for i in entropy_samples.keys()],axis=0)
        plt.bar(range(start,start+len(avg_diff[gene])),avg_diff[gene],color='C0')
    plt.ylabel('Entropy difference')
    plt.xlabel('Position')
    plt.title('Average window entropy difference on '+chr)
    #plt.savefig(file) 
    fig.clf(True)
    plt.close('all')  
    return avg_diff

def entropy_kmers_plot(annot,entropy,chr,file):
    fig,ax = plt.subplots()
    displ = 0
    for gene in entropy.keys():
        plt.bar(range(annot['Start'].iloc(gene),annot['Start'].iloc(gene)+len(entropy[gene])),entropy[gene])
        displ += len(entropy[gene])
    plt.figure(figsize=(12,7))
    plt.ylabel('Entropy')
    plt.xlabel('Position')
    plt.title('Positional entropy for genes at '+chr)
    plt.savefig(file)
    fig.clf(True)
    plt.close('all')
    return None

def distance_plot(annot,distance,l,chr,file,d_type):
    #print('Producing distance plot')
    avg_distance = dict()
    for gene in list(annot['gene_id']):
        avg_distance[gene] = np.average([distance[i][gene]*math.log2(l[gene]) for i in distance.keys()])
    fig, ax = plt.subplots()
    ax.bar(annot['Start'],list(avg_distance.values()),width=1)
    fig.set_size_inches(12, 7)
    ax.set_ylabel(d_type)
    ax.set_xlabel('Position')
    ax.set_title('Average '+d_type+ ' on '+chr)
    plt.savefig(file)
    fig.clf(True)
    plt.close()
    return avg_distance

def windows_distance_plot(annot,distance,chr,file,d_type):
    #print('Producing window distance plot')
    fig,ax = plt.subplots()
    avg_distance = dict()
    plt.figure(figsize=(12,7))
    for gene in list(annot['gene_id']):
        start = int(annot.loc[annot['gene_id']==gene,'Start'].iloc[0])
        #print(start)
        avg_distance[gene] = np.average([distance[i][gene] for i in distance.keys()],axis=0)
        plt.bar(range(start,start+len(avg_distance[gene])),avg_distance[gene],color='C0')
    plt.ylabel(d_type)
    plt.xlabel('Position')
    plt.title('Average '+d_type+' on '+chr)
    #plt.savefig(file) 
    fig.clf(True)
    plt.close('all')  
    return avg_distance

def distr_variants_genes(annot, vcf, res, file,cut,chr):
    res = res.sort_values('value',ascending=False).head(cut)
    fig,ax = plt.subplots()
    count = dict()
    length = dict()
    with vcfpy.Reader.from_path(vcf) as reader:
        for gene in res.loc[:,'gene']:
            annot_gene = annot[annot['gene_id']==gene]
            index = list(annot_gene.index)[0]
            length[gene] = annot_gene.at[index,'End']-annot_gene.at[index,'Start']-1
            chromosome = str(annot_gene.at[index,'Chromosome'])[3:len(str(annot_gene.at[index,'Chromosome']))]
            #print(chromosome)
            count[gene] = len(pd.DataFrame(reader.fetch(chromosome,annot_gene.at[index,'Start'],annot_gene.at[index,'End'])).index)

    #print(count.values())
    fig.set_size_inches(12, 7)
    fig.subplots_adjust(right=0.8)
    p1 = ax.bar(range(1,len(count)+1),count.values(),label='Variants count')
    ax2 = ax.twinx()
    ax3 = ax.twinx()
    ax3.spines.right.set_position(("axes", 1.1))
    ax2.set_ylim(0,max(res['value'])+0.01)
    p2, = ax2.plot(range(1,len(count)+1), res['value'],marker='.',lw = 1,color='green',label='Normalized Entropy')
    ax.set_ylabel('Number of variants')
    ax.set_title('Variants distribution top entropy genes '+chr)
    ax2.set_ylabel('Normalized Entropy')
    ax3.set_ylim(0,max(length.values())+10000)
    p3, = ax3.plot(range(1,len(count)+1),length.values(),label='Gene length',color='red',lw=1,marker='.')
    ax3.set_ylabel('Gene length')
    ax.set_xlabel('Ranked genes')
    ax.legend(handles=[p1,p2,p3])
    plt.savefig(file) 
    fig.clf(True)
    plt.close('all')
    return None
  
