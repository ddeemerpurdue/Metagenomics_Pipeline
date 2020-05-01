import os
import sys
import subprocess

# file = sys.argv[1]
""" So if we import this into another module, you still
want to get the information from the file, but not execute
the code """


def read_sample_file(file):
    with open(file) as f:
        line = f.readline().split('\t')
        if line[0].lower() == 'assembly':
            assembly = str(line[1].strip())
        else:
            print('Error, assembly line not formatted correctly.')
        line = f.readline().split('\t')
        if line[0] == 'abundance':
            abundance = str(line[1].strip())
        elif line[0] == 'abund_list':
            abundance = line[1]
        else:
            print('Error: abundance file not formatted correctly.')
        line = f.readline().split('\t')
        if line[0] == 'output':
            output = str(line[1].strip())
        else:
            print('Error: output file not formatted correctly.')
    print('Read in file information.')
    return assembly, abundance, output


if __name__ == "__main__":
    os.system('module load bioinfo')
    os.system('module load MaxBin/2.2.3')
    assembly, abundance, output = read_sample_file()
    subprocess.check_call(['run_MaxBin.pl', '-contig', assembly, '-abund',
                           abundance, '-out', output])
    print('Finished!')
