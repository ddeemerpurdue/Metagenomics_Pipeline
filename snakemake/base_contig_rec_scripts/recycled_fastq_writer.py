'''
Program designed to take in as input a raw assembly fasta file,
a standardized ANI file, and the sample name (column 0 sample in
ANI file). Output are {sample}.{bin}.fasta files.
Note: This file should also be found in the Contig_Purifier repository

Example usage:
$ python recycledBinFastaWriter.py -r repatriatedFile.txt -s assemblyFile.fasta -q sampleName
'''

import argparse


def get_contig_names_from_recycled(repat, sample):
    '''
    Function to read in a repatriated output and write
    a dictionary in the form:
    {node_identifier:[bin_id, reference_sample]}
    '''
    contig_names = {}
    with open(repat) as i:
        line = i.readline()  # Skip header
        line = i.readline()
        while line:
            line = line.split('\t')
            if sample == line[0]:
                contig_ref = line[1]
                contig = line[2]
                if contig.startswith('>'):
                    pass
                else:
                    contig = '>' + contig
                contig_bin = line[4].strip()
                contig_names[contig] = [contig_bin, contig_ref]
            else:
                pass
            line = i.readline()
    return contig_names


def write_to_fasta(repat, sample, assembly):
    '''
    Function that searches an assembly based on a dictionary
    with node identifiers as keys and bin/reference as values
    and writes a new fasta file
    '''
    seqdic = get_contig_names_from_recycled(repat, sample)
    with open(assembly) as assem:
        line = assem.readline().strip()
        while line:
            if (line.startswith('>') and line in seqdic):
                writefile = f"{sample}.{seqdic[line][0]}.fasta"
                with open(writefile, 'a') as out:
                    out.write(line + '\n')
                    line = assem.readline()
                    while not line.startswith('>'):
                        out.write(line)
                        line = assem.readline()
            else:
                line = assem.readline().strip()


def write_new_ani(repat, sample, oldani):
    output = str(oldani.rsplit('.', 1)[0]) + str(".new.txt")
    new_contigs = get_contig_names_from_recycled(repat, sample)
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
    parser.add_argument("-r", "--Repat", help="Repatriated file",
                        required=True)
    parser.add_argument("-s", "--Assembly", help="Raw Assembly", required=True)
    parser.add_argument("-q", "--Query", help="Query Sample", required=True)
    parser.add_argument("-o", "--OldANI", help="Original ANI file",
                        required=False)
    argument = parser.parse_args()
    write_to_fasta(argument.Repat, argument.Query, argument.Assembly)
    if argument.OldANI:
        write_new_ani(argument.Repat, argument.Query, argument.OldANI)
    else:
        pass
