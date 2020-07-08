'''
Program designed to take in the raw input from the fastANI snakemake
pipeline and standardize it with bin information from RefineM output.
Input requires a tab-delimited list of the binfiles as well.
Example useage:
$ python add_bins.py -a ani.txt -o out.txt -b bin1.txt bin2.txt
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
for value in argument.Bins:
    binlist.append(value)


def get_bindict(binfile):
    # Start with query bins
    bindict = {}
    with open(binfile) as qf:
        line = qf.readline()
        while line:
            node = line.split('\t')[0].split('_')[1]
            node = str(int(node))
            bin_name = line.split('\t')[1].strip()
            bindict[node] = bin_name
            line = qf.readline()
    return bindict


def get_all_bindicts(binlist):
    master = {}
    for bins in binlist:
        try:
            name = bins.rsplit('/',)[1].split('-')[0]
        except IndexError:
            name = bins.split('-')[0]
        master[name] = get_bindict(bins)
    return master


def append_bins_to_ani(binlist, anifile):
    writelist = []
    bins = get_all_bindicts(binlist)
    # Parse through ANI file
    with open(anifile) as af:
        constant_line = af.readline()
        while constant_line:
            line = constant_line.split('\t')
            quer_num = line[0].split('/')[2].split('_')[1]
            ref_num = line[1].split('/')[2].split('_')[1]
            query = line[0].split('/')[1]
            reference = line[1].split('/')[1]
            quer_node = line[0].split('/')[2].rsplit('.', 1)[0]
            reference_node = line[1].split('/')[2].rsplit('.', 1)[0]
            ani = line[2]
            orths = line[3]
            total = line[4].strip()
            # print(f"Query={query}, Ref={reference}, Nums={quer_num}\t{ref_num},\nNodes:\t{quer_node}\t{reference_node}\nANI:\t{ani}")
            # See if there's a bin for the nodes
            if quer_num in bins[query]:
                quer_bin = bins[query][quer_num]
            else:
                quer_bin = 'NoBin'
            # Repeat for reference
            if ref_num in bins[reference]:
                ref_bin = bins[reference][ref_num]
            else:
                ref_bin = 'NoBin'
            # Write all information to a file
            new_info = '\t'.join([query, reference, quer_node, reference_node,
                                  ani, orths, total, quer_bin, ref_bin])
            writelist.append(new_info)
            constant_line = af.readline()
    return writelist


def write_new_ani(binlist, anifile, output):
    writelist = append_bins_to_ani(binlist, anifile)
    with open(output, 'w') as o:
        for line in writelist:
            o.write(line + '\n')
    print('Done')


if __name__ == "__main__":
    write_new_ani(binlist, argument.ANI, argument.Output)
