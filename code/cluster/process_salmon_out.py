import pandas as pd
import numpy as np
from dppd import dppd
import glob
import os
import sys

'''
description
'''

print('- initialize tpm and raw count matricies')
quant_file = pd.read_csv('results/salmon_out/SRR5164008/quant.sf', sep = '\t', header = 0)
tpm = pd.DataFrame(quant_file['Name'])
counts = pd.DataFrame(quant_file['Name'])

print('- iterate over quant files and fill tpm and count matricies')
# filenames = [i for i in glob.glob('results/salmon_out/*/quant.sf')]
# quants = [pd.read_csv(file, sep = "\t") for file in filenames]

for file in glob.glob('results/salmon_out/*/quant.sf'):
	
	# get current run and quant info
	current_run = file.split('/')[2]
	current_quant = pd.read_csv(file, sep = "\t")
	
	# add current info to master matricies
	tpm[current_run] = current_quant['TPM']
	counts[current_run] = current_quant['NumReads']

print('- saving TPM and count matricies')
tpm.to_csv('results/salmon_processed/tpm.tsv', sep = '\t', index = False)
counts.to_csv('results/salmon_processed/counts.tsv', sep = '\t', index = False)
