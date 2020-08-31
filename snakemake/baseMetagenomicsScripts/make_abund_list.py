"""
Program that takes a directory as input and then finds
all files that end with *.abundance.txt and creates a
new abund_list.txt for MaxBin input

Example useage:
$ python make_abund_list.py -i <path/to/input/dir/> -o <outputfile.txt>
"""
import os

infile = str(snakemake.input)
outfile = str(snakemake.output)


def write_abundance_list(snakemake_in, snakemake_out):
    mylist = []
    infile = snakemake_in.split()
    for f in infile:
        if str(os.path.isfile(f)):
            mylist.append(f)
        else:
            return f"File: {f} does not exist!!!"
    with open(outfile, 'w') as o:
        for f in mylist:
            o.write(str(f) + '\n')
    return 0


if __name__ == "__main__":
    write_abundance_list(infile, outfile)
