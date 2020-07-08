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
                contig_ref = line.split('\t')[1]
                contig = line.split('\t')[2]
                contig_bin = line.split('\t')[4].strip()
                contig_names[contig] = [contig_bin, contig_ref]
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
                writefile = f"{sample}.{seqlist[match][0]}.fasta"
                with open(writefile, 'a') as out:
                    out.write(line + '\n')
                    line = assem.readline()
                    while not line.startswith('>'):
                        out.write(line)
                        line = assem.readline()
            else:
                line = assem.readline().strip()


def write_new_ani(anifile, sample, oldani):
    output = str(oldani.split('.')[0]) + str(".new.txt")
    new_contigs = get_contig_names_from_recycled(anifile, sample)
    with open(output, 'w') as o:
        with open(oldani) as oldani:
            line = oldani.readline()
            while line:
                old_ref = line.split('\t')[1]
                old_node = line.split('\t')[2]
                if old_node in new_contigs:
                    newref = new_contigs[old_node][1]
                    if newref == old_ref:
                        writeline = "\t".join(line.split('\t')[0:7])
                        writeline = writeline + "\t" + new_contigs[old_node][0] + "\t"
                        writeline = writeline + line.split('\t')[8] + "\n"
                        o.write(writeline)
                    else:
                        pass
                else:
                    o.write(line)
                line = oldani.readline()
    return "Done"


if __name__ == "__main__":
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-a", "--ANI", help="ANI master file", required=True)
    parser.add_argument("-s", "--Assembly", help="Raw Assembly", required=True)
    parser.add_argument("-q", "--Query", help="Query Sample", required=True)
    parser.add_argument("-o", "--OldANI", help="Original ANI file",
                        required=False)
    argument = parser.parse_args()
    write_to_fasta(argument.ANI, argument.Query, argument.Assembly)
    if argument.OldANI:
        write_new_ani(argument.ANI, argument.Query, argument.OldANI)
    else:
        pass
