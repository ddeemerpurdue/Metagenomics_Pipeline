'''
Author: Dane
Module that is snakemake-specific and takes in 
'''
import os
import sys

# Initialize input parameters from snakemake into a list
in_contigs = str(snakemake.input)
out_list = str(snakemake.output)
# Put the input inside of a list
in_contigs = in_contigs.split()
out_list = out_list.split()

print(f"{in_contigs}\n{out_list}")

# in_contigs contains a list of all the *.complete.tkn files created from
# split_filtered_contigs rule, while out_list contains output files where
# each line is the path to a contig
for in_contig, out in zip(in_contigs, out_list):
    # Remove basename to get location of .fa files
    directory = in_contig.rsplit('/', 1)[0]
    line = []
    # Go through every .fa file in the directory
    for file in os.listdir(directory):
        if file.endswith('.fa'):
            line.append(file)
        else:
            pass
    line = sorted(line)
    # Write every .fa file from the directory to the output file (lists/{sample}.{length}.txt
    with open(out, 'w') as o:
        for l in line:
            o.write(str(directory) + '/' + str(l))
            o.write('\n')
