"""
Program that takes a directory as input and then finds
all files that end with *.abundance.txt and creates a
new abund_list.txt for MaxBin input

Example useage:
$ python make_abund_list.py -i <path/to/input/dir/> -o <outputfile.txt>
"""
import os

# import argparse

''' Arguments '''
# parser = argparse.ArgumentParser(description="Parser")
# parser.add_argument("-i", "--Input",
#                     help="Directory containing .abundance.txt",
#                     required=True)
# parser.add_argument("-o", "--Output",
#                     help="Output file to write to",
#                     required=True)
# argument = parser.parse_args()
''' Arguments '''
infile = str(snakemake.input)
outfile = str(snakemake.output)

wdir = os.getcwd()

mylist = []
infile = infile.split()
print(infile)
for f in infile:
    print(f)
    if f.endswith('abundance.txt'):
        mylist.append(f)

print(mylist)

with open(outfile, 'w') as o:
    for f in mylist:
        o.write(str(f) + '\n')
