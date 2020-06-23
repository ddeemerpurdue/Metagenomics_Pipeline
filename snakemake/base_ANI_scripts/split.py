import os
import subprocess

# Initialize input parameters from snakemake into a list
in_contigs = str(snakemake.input)
paramers = snakemake.params.folder
out_contigs = str(snakemake.output)
# Put the input inside of a list
in_contigs = in_contigs.split()
out_contigs = out_contigs.split()

for i, p in zip(in_contigs, paramers):
    print(f"{i}\n{p}")
    print(type(i))
    print(type(p))
    cmd = f"./split_mfa.sh {i} {p}"
    os.system(cmd)

for i in out_contigs:
    print(i)
    with open(i, 'w') as o:
        o.write(f"{i} successfully split!")
