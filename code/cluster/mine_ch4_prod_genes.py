import pandas as pd
import numpy as np
from Bio import SeqIO
from pysam import FastaFile
from dppd import dppd
import glob
import os
import sys

'''
description
'''

print('- taking in arguments from standard input')
m_smithii_1_filepath = sys.argv[1]
m_smithii_2_filepath = sys.argv[2]
m_smithii_3_filepath = sys.argv[3]
m_smithii_4_filepath = sys.argv[4]
output_filepath = sys.argv[5]

print('- iterating over various M. smithii genomes to obtain genes involved in methane production')

def genome_parse(filepath):
	for sequence in  SeqIO.parse(filepath, 'fasta'):
		gene_name = sequence.description.split(' ')[1]
		if any(s in gene_name for s in ('fwdA', 'fwdB', 'fwdD', 'fwdE', 'fwdF', 'fwdG', 'ftr', 
										'hmd', 'mch', 'mcrA', 'mcrB', 'mcrG', 'mer', 'mtrA', 
										'mtrB', 'mtrC', 'mtrD', 'mtrE', 'mtrF', 'mtrG')):
			# ch4_prod_genes[gene_name] = sequence.seq
			current_entry = '>' + sequence.description.split(' ')[0] + ';' + sequence.description.split(' ')[1]
			current_seq = sequence.seq

			with open(output_filepath, 'a') as f:
				f.write(current_entry)
				f.write('\n')
				f.write(str(current_seq))
				f.write('\n')

genome_parse(m_smithii_1_filepath)
genome_parse(m_smithii_2_filepath)
genome_parse(m_smithii_3_filepath)
genome_parse(m_smithii_4_filepath)
