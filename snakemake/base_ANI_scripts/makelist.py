# Step 1:

import os
import sys

# Initialize input parameters from snakemake into a list
in_contigs = str(snakemake.input)
out_list = str(snakemake.output)
# Put the input inside of a list
in_contigs = in_contigs.split()
out_list = out_list.split()

print(f"{in_contigs}\n{out_list}")


for directory, out in zip(in_contigs, out_list):
    d = directory.rsplit('/', 1)[0]
    line = []
    print(d)
    print(out)
    for file in os.listdir(d):
        if file.endswith('.fa'):
            line.append(file)
        else:
            pass
    line = sorted(line)
    with open(out, 'w') as o:
        for l in line:
            o.write(str(d) + '/' + str(l))
            o.write('\n')

