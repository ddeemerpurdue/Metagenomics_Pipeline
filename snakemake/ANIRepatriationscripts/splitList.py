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
out_directory = str(snakemake.params.directory)
split_size = int(snakemake.params.split_size)
# Put the input inside of a list
in_list = in_list.split()
out_list = out_list.split()


# Step 1: Split each file into N smaller files
for file in in_list:
    prefix = out_directory + '/' + os.path.basename(file).split('.AllContigsList.txt')[0] + str('_')
    cmd = f"split -l$((`wc -l < {file}`/ {split_size})) {file} {prefix}"
#    cmd = f"split -n 10 {file} {prefix}"
    print(f"Splitting {file} into {split_size} equal chunks.")
    os.system(cmd)

