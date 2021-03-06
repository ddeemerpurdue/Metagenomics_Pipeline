'''
Author: Dane
Snakemake specific script that makes a call to filter_fast.py to filter
assemblies to only contain contigs >= {params.length}
'''
# Step 1
from filterSeqLength import *

# Initialize input parameters from snakemake into a list
# Input is list of all the split files
in_list = str(snakemake.input.samples)
in_list = in_list.split()
# Params
param = int(snakemake.params.length)
# Outputs
out_list = str(snakemake.output.outputs)
out_list = out_list.split()

for i, o in zip(in_list, out_list):
    filter_fasta(i, param, o)

print('Done!')
