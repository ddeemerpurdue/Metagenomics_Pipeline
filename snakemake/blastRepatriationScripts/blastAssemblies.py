'''
Given 2 folders, one for a set of bins and another for a set of assemblies
(labeled - {bin_number}.fasta), blast binners against one another.

Example usage:
$ python blastAssemblies.py -a ReferenceFiles/particle/
-b Fastas/ParticleTRA90R99ANIT95M20/
'''
import os
import sys


def blast_all_files(binsLoc, assemblyLoc):
    '''

    '''
    assembly_files = []
    bin_files = []
    for file in os.listdir(assemblyLoc):
        if file.endswith('.fasta'):
            full_filename = os.path.join(assemblyLoc, file)
            assembly_files.append(full_filename)
    for file in os.listdir(binsLoc):
        if file.endswith('.fasta'):
            if file.
            full_filename = os.path.join(binsLoc, file)
            bin_files.append(full_filename)
    # Sort them
    assembly_files = sorted(assembly_files)
    bin_files = sorted(bin_files)
    # Create strings to use as commands
    for bin_file, assem_file in zip(bin_files, assembly_files):
        bin_number = os.path.basename(assem_file).split(".")[0]
        command = f"blastn -query {bin_file} -subject {assem_file} -outfmt '6' -out blast_bin_p/blastn.{bin_number}.txt"
        os.system(command)


if __name__ == "__main__":
    blast_all_files(sys.argv[1], sys.argv[2])
