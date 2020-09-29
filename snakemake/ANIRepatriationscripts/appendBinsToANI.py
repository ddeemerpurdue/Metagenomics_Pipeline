'''
Author: Dane
Program designed to take in the raw input from the fastANI snakemake
pipeline and standardize it with bin information.
Input requires a tab-delimited list of the binfiles as well.
Note: Each bin to node identification file (1 per sample) NEEDS to be
labeled in the format: <sample>.*
where <sample> matches the sample completeley in the ANI file.
For example, if in the ani file that sample names are 'particle' and
'supernatant', the bin identification files must have a prefix starting
with 'particle.' and 'supernatant.', such as 'particle.originalBins.txt'
and 'supernatant.originalBins.txt'.
Example useage:
$ python add_bins.py -a ani.txt -o out.txt -b bin1.txt bin2.txt
Again: bin files need to be in the specified format: sample.*
'''
import os
import time


def get_bindict(binfile):
    '''
    Open a bin identification file and output a dictionary in
    the form: dict[node] = bin_number
    '''
    bindict = {}
    with open(binfile) as qf:
        line = qf.readline()
        while line:
            bin_name = line.split('\t')[0]
            node = line.split('\t')[1].strip()
            node = node.strip('>')
            bindict[node] = bin_name
            line = qf.readline()
    return bindict


def get_all_bindict(binlist):
    '''
    Given a list of bin identification files, read each in a dictionary
    in the form of: dict[name][node] = bin_number
    This is a nested dictionary 1 level higher than get_bindict.
    '''
    all_dict = {}
    for file in binlist:
        name = os.path.basename(file)
        name = name.split('.')[0]
        dic = get_bindict(file)
        all_dict[name] = dic
    return all_dict


def append_bins_to_ani(binfile, anifile, output, logfile):
    '''
    Append bin identifiers to ANI's raw output.
    '''
    now = time.localtime(time.time())
    writelist = []
    bins = get_all_bindict(binfile)
    match_count = 0
    nonmatch_count = 0
    # Parse through ANI file
    with open(output, 'w') as o:
        with open(anifile) as af:
            constant_line = af.readline()
            while constant_line:
                line = constant_line.split('\t')
                # Below, grab sample name out of first 2 fields in raw ANI output
                quer = line[0].split('/')[5]  # With pipeline, is 5th entry
                refer = line[1].split('/')[5]
                # Do not need ANI results within samples (e.g., R1 vs. R1)
                if quer != refer:
                    quer_node = line[0].split('/')[6].rsplit('.', 1)[0]
                    reference_node = line[1].split('/')[6].rsplit('.', 1)[0]
                    ani = line[2]
                    orths = line[3]
                    total = line[4].strip()
                    if quer_node in bins[quer]:
                        quer_bin = bins[quer][quer_node]
                        match_count += 1
                    else:
                        quer_bin = 'NoBin'
                        nonmatch_count += 1
                    # Repeat for reference
                    if reference_node in bins[refer]:
                        ref_bin = bins[refer][reference_node]
                        match_count += 1
                    else:
                        ref_bin = 'NoBin'
                        nonmatch_count += 1
                    # Write all information to a file
                    new_info = '\t'.join([quer, refer, quer_node,
                                          reference_node, ani, orths,
                                          total, quer_bin, ref_bin]) + '\n'
                    o.write(new_info)
                else:
                    pass
                constant_line = af.readline()
    # Assert that some bins matches
    # assert match_count > 0, "No contigs matched to a bin!"
    # Finally, write the logfile:
    with open(logfile, 'w') as lg:
        lg.write(f"Time started: {time.asctime(now)}\n")
        for bin_dic in bins.keys():
            lg.write(
                f"Sample {bin_dic} contains {len(bins[bin_dic])} contigs.\n")
        lg.write(
            f"{match_count} contigs were identified and {nonmatch_count} were labeled NoBin.\n")
        now = time.localtime(time.time())
        lg.write(f"Time finished: {time.asctime(now)}\n")
    return writelist


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-b", "--Bins", required=True, nargs="*",
                        help="One or more files with contig name and bin info")
    parser.add_argument("-a", "--ANI", help="ANI file created from FASTANI",
                        required=True)
    parser.add_argument("-o", "--Output", help="Standardized ANI file to write to",
                        required=True)
    parser.add_argument("-l", "--Log", help="Log file to write to",
                        required=True)
    argument = parser.parse_args()
    append_bins_to_ani(argument.Bins, argument.ANI,
                       argument.Output, argument.Log)
