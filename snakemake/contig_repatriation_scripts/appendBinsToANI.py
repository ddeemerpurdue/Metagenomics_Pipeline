'''
Author: Dane
Program designed to take in the raw input from the fastANI snakemake
pipeline and standardize it with bin information from RefineM output
(or whatever the most up-to-date output is from).
Input requires a tab-delimited list of the binfiles as well.
Example useage:
$ python add_bins.py -a ani.txt -o out.txt -b bin1.txt bin2.txt

Note: bin files need to be in the specified format: sample.Bins.txt
'''

import argparse

parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-b", "--Bins", required=True, nargs="*",
                    help="One or more files with contig name and bin info")

parser.add_argument("-a", "--ANI", help="ANI file created from FASTANI",
                    required=True)
parser.add_argument("-o", "--Output", help="Standardized ANI file to write to",
                    required=True)
argument = parser.parse_args()


binlist = []
for b in argument.Bins:
    binlist.append(b)
print(binlist)


def get_bindict(binfile):
    # Start with query bins
    bindict = {}
    with open(binfile) as qf:
        line = qf.readline()
        while line:
            bin_name = line.split('\t')[0].split('.')[1].split('_')[0]
            node = line.split('\t')[1].strip()
            node = node.strip('>')
            bindict[node] = bin_name
            line = qf.readline()
    return bindict


def get_all_bindict(binlist):
    all_dict = {}
    for file in binlist:
        dic = get_bindict(file)
        key = str(file).split('.')[0]
        all_dict[key] = dic
    return all_dict


def append_bins_to_ani(binfile, anifile, output):
    writelist = []
    bins = get_all_bindict(binfile)
    # Parse through ANI file
    with open(output, 'w') as o:
        with open(anifile) as af:
            constant_line = af.readline()
            while constant_line:
                line = constant_line.split('\t')
                quer = line[0].split('/')[1].split('_')[0]
                refer = line[1].split('/')[1].split('_')[0]
                if quer != refer:
                    quer_node = line[0].split('/')[2].rsplit('.', 1)[0]
                    reference_node = line[1].split('/')[2].rsplit('.', 1)[0]
                    ani = line[2]
                    orths = line[3]
                    total = line[4].strip()
                    if quer_node in bins[quer]:
                        quer_bin = bins[quer][quer_node]
                    else:
                        quer_bin = 'NoBin'
                    # Repeat for reference
                    if reference_node in bins[refer]:
                        ref_bin = bins[refer][reference_node]
                    else:
                        ref_bin = 'NoBin'
                    # Write all information to a file
                    new_info = '\t'.join([quer, refer, quer_node,
                                          reference_node, ani, orths,
                                          total, quer_bin, ref_bin]) + '\n'
                    o.write(new_info)
                else:
                    pass
                constant_line = af.readline()
    return writelist


if __name__ == "__main__":
    append_bins_to_ani(binlist, argument.ANI, argument.Output)
