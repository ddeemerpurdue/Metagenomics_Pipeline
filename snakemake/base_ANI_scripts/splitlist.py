'''
Author: Dane
Snakemake-specific python file that takes in a file, in_list, that contains
each .fa file path for the filtered assembly (each line is a .fa entry).
The output is a file with the same name as in_list except instead of ending
in .txt it ends in _{splits}.txt. The goal is to split the in_list file into
n separate files in order to parallelize downstream processes.
'''

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

