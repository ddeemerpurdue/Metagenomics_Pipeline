"""
Program designed to take an ANI input file annotated
with bin IDs and repatriate contigs

Input:
Tab-delimited file with fields:
Input from add_ani_bins.py!
Query, Reference, Qnode, Rnode, ANI, Orthos, Total, Qbin, Rbin

Example usage:
$ python aniContigRecycler.py -a ANIFILE.txt -t 90 -m 200 -o output90T200M.txt

"""
import pandas as pd
import argparse


def get_cois(file, threshold):
    '''
    First function: read in the ANI file and return 3 separate
    data frames depending on the conditions:
    i) Contig matches (CMs) where both contigs are not binned
    ii) CMs without bin in Query but with bin in Reference
    iii) CMs with bin matchs in both Query and Reference
    '''
    cnames = ['Query', 'Reference', 'Qnode', 'Rnode', 'ANI',
              'Orthos', 'Total', 'Qbin', 'Rbin']
    df = pd.read_csv(file, sep='\t',
                     names=cnames)
    df = df[df.ANI >= int(threshold)]
    nobin_nobin = df[(df.Qbin == 'NoBin') & (df.Rbin == 'NoBin')]
    nobin_bin = df[(df.Qbin == 'NoBin') & (df.Rbin != 'NoBin')]
    bin_bin = df[(df.Qbin != 'NoBin') & (df.Rbin != 'NoBin')]
    return nobin_nobin, nobin_bin, bin_bin


def write_cois(file, threshold):
    '''
    Optional function to write the subsetted dataframes
    to a tab-delimited file
    '''
    nb_nb, nb_b, b_b = get_cois(file, threshold)
    nb_nb.to_csv("Nobin_Nobin.txt", sep='\t', index=False)
    nb_b.to_csv("Nobin_Bin.txt", sep='\t', index=False)
    b_b.to_csv("Bin_Bin.txt", sep='\t', index=False)
    print('Done!')


def create_nobin_dictionary(file, threshold):
    '''
    Function designed to create a 3x-nested-dictionary in the
    form of: {Query{Reference{Contig[Rbin1, Rbin2, ...]}}}
    This function identifies where query contig nonbinners match
    up to bins in all other files.
    '''
    # df = pd.read_csv(file, sep='\t')
    _, df, _ = get_cois(file, threshold)
    large_dictionary = {}
    for quer in df.Query.unique():
        reference_dict = {}
        for refer in df.Reference.unique():
            # No need to compare to own sample
            if quer == refer:
                pass
            else:
                # Only look at subset of df corresponding to for loop
                a = df[(df.Query == quer) & (df.Reference == refer)]
                contig_dict = {}
                for i in range(a.shape[0]):
                    # Query contig is assigned 'qcontig'
                    qcontig = a.iloc[i, 2]
                    # Reference bin is assigned 'rbin'
                    rbin = a.iloc[i, 8]
                    # This loop is for multi-matchers
                    if (qcontig in contig_dict) and (rbin not in contig_dict[qcontig]):
                        contig_dict[qcontig].append(rbin)
                    else:
                        contig_dict[qcontig] = [rbin]
                reference_dict[refer] = contig_dict
        large_dictionary[quer] = reference_dict
    return large_dictionary


def create_bin_dictionary(file, threshold):
    '''
    Function designed to take as input the subset of ANI information
    where both query and reference-matched contigs are assigned bins.
    The data is captured in a 3x-nested dictionary in the form of:
    {Query{Reference{Qbin:[Rbin, # of matches]}}}. The Rbin(s) from
    this function can be compared to the create_nobin_dictionary
    function to assign non-binners a query bin.
    '''
    _, _, df = get_cois(file, threshold)
    large_dictionary = {}
    all_num_matches = []
    for quer in df.Query.unique():
        reference_dict = {}
        for refer in df.Reference.unique():
            if quer == refer:
                pass
            else:
                contig_dict = {}
                for bins in df.Qbin.unique():
                    a = df[(df.Query == quer) &
                           (df.Reference == refer) &
                           (df.Qbin == bins)]
                    for i in range(a.shape[0]):
                        qbin = a.iloc[i, 7]
                        rbin = a.iloc[i, 8]
                        # Make a dictionary where key == qbin and values
                        # are associated rbin
                        if qbin in contig_dict:
                            contig_dict[qbin].append(rbin)
                        else:
                            contig_dict[qbin] = [rbin]
        #         # Count occurences of each Rbin per Qbin
        # # At the end, contig_dict would contain all unique Qbins
        # # for Quer/Refer
                # Go through each Qbin key
                for key in contig_dict.keys():
                    cnt_list = []
                    # Loop through the set of values (unique)
                    for l in set(contig_dict[key]):
                        all_num_matches.append(contig_dict[key].count(l))
                        # If it's first value...
                        if cnt_list == []:
                            cnt_list.append(l)
                            cnt_list.append(contig_dict[key].count(l))
                            # Ex. ['Bin013', 13]
                        else:
                            # If multiple Rbins
                            if contig_dict[key].count(l) > cnt_list[1]:
                                cnt_list = []
                                cnt_list.append(l)
                                cnt_list.append(contig_dict[key].count(l))
                            # Ex. ['Bin013', 13, 'Bin011', 4]
                            else:
                                pass
                    # Now have key go from each occurence to counts
                    contig_dict[key] = cnt_list
                    # Ex. {'qBin002': ['rBin13', 13, 'rBin011', 4]}
                    # Proceed for all qBins
                reference_dict[refer] = contig_dict
        large_dictionary[quer] = reference_dict
    return large_dictionary, all_num_matches


def map_unbinned_contigs(file, threshold, matches):
    '''
    Function designed to pair up non-binners with a new home :)
    if the situation passes user-defined parameters
    '''
    lines = [f"query\treference\tcontig\tbin_to_match\tMATCH\n"]
    nonbins = create_nobin_dictionary(file, threshold)
    binners, _ = create_bin_dictionary(file, threshold)
    # 1 - Loop through all queries
    # These are actually the same for both dictionaries
    for query in nonbins.keys():
        # 2 - Loop through all references
        # Again, same for both dictionaries
        for reference in nonbins[query].keys():
            # 3 Go through all unbinned contigs
            for contig in nonbins[query][reference].keys():
            # 4 Capture their nearest bin match
                bins_to_match = nonbins[query][reference][contig]
            # Ex. {R1: {R2: {Node_10_blah: ['bin015', 'bin001']}}}
            # 5 Since sometimes multiple matches, loop through each match
                for bin_to_match in bins_to_match:
            # Ex. 'bin015' - Saying contig from R1 matched bin 15 from R2
            # 6 Search 'binners' to find match
            # Where do other contigs from R1 match R2 bins?
                    # Go through all query bins from bin-bin dict
                    for quer_bin in binners[query][reference].keys():
            # 7 If there's a match between Rbin and Qbin/Rbin
            # If reference from non-binner matches query:REFERENCE in binner AND
            # that match occurs at least N+ times
                        if (bin_to_match in binners[query][reference][quer_bin] and
                                binners[query][reference][quer_bin][1] > int(matches)):
            # This would match binners {R1: {R2: {Qbin:['Rbin']}}} -- 'Rbin'
                            bin_match = quer_bin
                            lines.append(f"{query}\t{reference}\t{contig}\t{bin_to_match}\t{bin_match}\n")
                            break
                        else:
                            bin_match = "No_Recycle"
                        
    return lines


def write_recycled_bins(file, threshold, matches, output):
    a = map_unbinned_contigs(file, threshold, matches)
    with open(output, 'w') as o:
        for line in a:
            o.write(line)
    return None


if __name__ == "__main__":
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-a", "--ANI", help="ANI master file", required=True)
    parser.add_argument("-t", "--Threshold", help="ANI threshold to consider \
    matches", default=90)
    parser.add_argument("-m", "--Matches", help="How many binning contigs must\
    match the reference bin", default=10, required=False)
    parser.add_argument("-o", "--Output", help="Output file to write to",
                        required=False)
    parser.add_argument("-w", "--Write", help="If provided, write files for \
                        binners, nonbinners, and binners/nonbinners",
                        required=False, action='store_true')
    argument = parser.parse_args()
    if argument.Write:
        write_cois(argument.ANI, argument.Threshold)
    else:
        write_recycled_bins(argument.ANI, argument.Threshold,
                            argument.Matches,argument.Output)
