# Step 1
from filter_fast import *

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
    write_fasta(i, param, o)

print('Done!')
