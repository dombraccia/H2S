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
dataset = sys.argv[1].split('/')[1]
print('-- dataset name: ' + dataset)
function = sys.argv[2]
print('-- function name: ' + function)
a_quant_file = glob.glob('results/salmon_out/' + dataset + '/' + function + '/*/quant.sf')[0]
quant_file = pd.read_csv(a_quant_file, sep = '\t', header = 0)
# quant_file = pd.read_csv('results/salmon_out/SRR5164008/quant.sf', sep = '\t', header = 0)
tpm = pd.DataFrame(quant_file[['Name', 'Length', 'EffectiveLength']])
counts = pd.DataFrame(quant_file[['Name', 'Length', 'EffectiveLength']])

print('- iterate over quant files and fill tpm and count matricies')
# filenames = [i for i in glob.glob('results/salmon_out/*/quant.sf')]
# quants = [pd.read_csv(file, sep = "\t") for file in filenames]

for file in glob.glob('results/salmon_out/' + dataset + '/' + function + '/*/quant.sf'):
	
	# get current run and quant info
	current_run = file.split('/')[4]
	current_quant = pd.read_csv(file, sep = "\t")
	
	# add current info to master matricies
	tpm[current_run] = current_quant['TPM']
	counts[current_run] = current_quant['NumReads']

print('- saving TPM and count matricies')

print('-- TPM file path:  ' + 'results/salmon_processed/' + dataset + '/' + function + '_' + 'tpm.tsv')
tpm.to_csv('results/salmon_processed/' + dataset + '/' + function + '_' + 'tpm.tsv', sep = '\t', index = False)
print('-- counts file path:  ' + 'results/salmon_processed/' + dataset + '/' + function + '_' + 'counts.tsv', sep = '\t')
counts.to_csv('results/salmon_processed/' + dataset + '/' + function + '_' + 'counts.tsv', sep = '\t', index = False)
