# Step 1:

import os
import sys

# Initialize input parameters from snakemake into a list
in_list = str(snakemake.input.lists)
out_list = str(snakemake.output)
# Put the input inside of a list
in_list = in_list.split()
out_list = out_list.split()


# Step 1: Split each file into 10 smaller files
for file in in_list:
    prefix = file.rsplit('.',1)[0] + str('_')
    cmd = f"split -l$((`wc -l < {file}`/9)) {file} {prefix}"
#    cmd = f"split -n 10 {file} {prefix}"
    os.system(cmd)

