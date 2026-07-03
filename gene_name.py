##obtain gene name from gene id
from pyranges.readers import read_gtf
import argparse
import pandas

parser = argparse.ArgumentParser(prog='gene_name.py')
parser.add_argument('-l',help='File of gene ids, one for each row',required=True)
annot = read_gtf('data/gencode.v49lift37.basic.annotation_protein_coding.gtf',as_df=True)
args = parser.parse_args()

ids = pandas.read_table(args.l)

names = [annot[annot['gene_id']==ids['gene_id'].iloc[i]]['gene_name'] for i in range(len(ids['gene_id']))]

ids['gene_names'] = names 
ids.to_csv('names_'+args.l,header=False,index=False,sep='\t')