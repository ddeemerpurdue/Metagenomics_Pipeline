'''
Program designed to take in as input a raw assembly fasta file,
a standardized ANI file, and the sample name (column 0 sample in
ANI file). Output are {sample}.{bin}.fasta files.
'''

import argparse


def get_contig_names_from_recycled(anifile, sample):
    contig_names = {}
    with open(anifile) as i:
        line = i.readline()  # Skip header
        line = i.readline()
        while line:
            if sample == line.split('\t')[0]:
                contig = line.split('\t')[2]
                contig_bin = line.split('\t')[4].strip()
                contig_names[contig] = contig_bin
            else:
                pass
            line = i.readline()
    return contig_names


def write_to_fasta(anifile, sample, assembly):
    seqlist = get_contig_names_from_recycled(anifile, sample)
    with open(assembly) as assem:
        line = assem.readline().strip()
        while line:
            match = line.strip('>')
            if (line.startswith('>') and match in seqlist):
                writefile = f"{sample}.{seqlist[match]}.fasta"
                with open(writefile, 'a') as out:
                    out.write(line + '\n')
                    line = assem.readline()
                    while not line.startswith('>'):
                        out.write(line)
                        line = assem.readline()
            else:
                line = assem.readline().strip()


if __name__ == "__main__":
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-a", "--ANI", help="ANI master file", required=True)
    parser.add_argument("-s", "--Assembly", help="Raw Assembly", required=True)
    parser.add_argument("-q", "--Query", help="Query Sample", required=True)
    argument = parser.parse_args()
    write_to_fasta(argument.ANI, argument.Query, argument.Assembly)
