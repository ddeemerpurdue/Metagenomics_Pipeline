'''
Author: Dane
Remove contamination from a metagenomic analysis by comparing bit-scores across contig
annotations vs. full bin annotations. Scores are obtained from CATBAT.
Requirements:
*.Bin2C.names.txt [Default CATBAT output]

Two threshold levels can be tuned to change confidence in taxonomy filtering:
--RemoveThreshold is the bit score required between the divergent contig annotation
and the overall bin annotation.
--AddThreshold is the bit score required between a non-binned contig and a potential
bin in order for the contig to join that bin.

Example usage:
$ python taxonFilter.py -i binFile.txt -c {sample}.C2C.names.txt -b {sample}.Bin2C.names.txt
-m 70 -a 50 -o outputFile.txt
'''


def read_taxonomy(taxonFile, bins=False):
    '''
    Function designed to read in *Bin2C.names.txt file
    or a *.C2C.names.txt file from running Bat or Cat
    (in CATBAT)
    '''
    bin_dic = {}
    with open(taxonFile) as i:
        line = i.readline()
        while line:
            if line.startswith('#'):
                line = i.readline()
            else:
                line = line.strip().split('\t')
                if bins is True:
                    name = line[0].split('.')[1]
                else:
                    name = line[0]
                binTaxon = line[5:]
                bin_dic[name] = binTaxon
                line = i.readline()
    return bin_dic


def read_bin_identifier(binfile):
    '''
    Read in a bin ID file in the format of:
    dict[bin_num] = node
    '''
    bin_dic = {}
    with open(binfile) as i:
        line = i.readline()
        while line:
            bin_num = line.split('\t')[0]
            contig = line.split('\t')[1].strip()
            bin_dic[contig] = bin_num
            line = i.readline()
    return bin_dic


def compare_taxonomies(contigTax, binTax, threshold):
    '''
    Function designed to compare taxonomies. The return value (True, False)
    can be used by another function to remove the contig from the bin and
    look for another bin to place the contig into.
    '''
    threshold = float(threshold)
    matches = 0
    total = 0
    '''
    Part A takes care of either the bin or the contig having higher levels
    of taxonomic annotations because zipping will only compare the first n
    values, where n is the length of the shorter list.
    '''
    if binTax == 'NoBin':
        return True
    for cval, bval in zip(contigTax, binTax):
        curcon = cval.split(':')
        curbin = bval.split(':')
        if curcon[0] == curbin[0]:
            matches += 1
        else:  # If the values don't match and are of high quality, flag it.
            try:  # Weird try loop for a few exceptions that don't follow pattern
                bitscore = float(curcon[1])
            except ValueError:
                bitscore = float(curcon[2])
            try:
                bin_bitscore = float(curbin[1])
            except ValueError:
                bin_bitscore = float(curbin[2])
            if (bitscore > threshold and bin_bitscore > threshold):
                # Here return true since they diverge significantly
                return True
            else:
                pass
        total += 1
    # If the contigTaxon is not 'NoBin' and does not diverge from the bin by
    # more than the threshold, no need to search for a new home.
    return False


def compare_new_bin_taxonomies(contigTaxon, bin_dic, threshold):
    '''
    Program that looks for a new home for a removed contig based
    on bit-scores. Return value is a list containing:
    contig repatriation score, bin confidence, and bin
    '''
    test_container = [0.0, 0.0, 'NoTaxon', 'NoTaxon']
    for binID in sorted(list(bin_dic.keys())):
        # Init counting variables
        bin_confidence = 0
        matches = 0.0
        score = 0.0
        binTax = bin_dic[binID]
        # Loop through taxonomies
        for cval, bval in zip(contigTaxon, binTax):
            curcon = cval.split(':')
            curbin = bval.split(':')
            if curbin[0] == 'NoTaxon':
                pass
            else:
                try:  # Weird try loop for a few exceptions that don't follow pattern
                    bitscore = float(curcon[1])
                except ValueError:
                    bitscore = float(curcon[2])
                try:
                    bin_bitscore = float(curbin[1])
                except ValueError:
                    bin_bitscore = float(curbin[2])
                if curcon[0] == curbin[0]:
                    # Each level
                    matches += 1.0
                    bin_confidence = (bitscore * float(matches))
                    # First round, max score ==  1
                    score += (bitscore * float(matches))
                    # Second round, max == (1 + 2)
                    # Third round, max = (1 + 2 + 3)
                else:  # Below, test if discrepancy too big to repair
                    # If both annotations are high level and don't match, then don't
                    # resolve this contig and make score = 0 so they don't collapse.
                    if (bitscore > threshold and bin_bitscore > threshold):
                        score = 0
                    # Otherwise, keep score where it is and potentially we can
                    # add this information
                    else:
                        pass
        # Test to see if we have a potential re-patriation candidate
        if score > 28.0:  # This accounts for first 7 levels of taxonomy
            if score > test_container[0]:
                test_container[0] = score
                test_container[1] = bin_confidence
                test_container[2] = binID
                test_container[3] = binTax
            elif score == test_container[0]:
                if bin_confidence > test_container[1]:
                    test_container[0] = score
                    test_container[1] = bin_confidence
                    test_container[2] = binID
                    test_container[3] = binTax
                else:
                    pass
        else:
            pass
    return test_container


def loop_contig_taxonomies(binTaxonFile, contigTaxonFile, binidfile,
                           removeThresh, addThresh, readme, output):
    '''
    Fx that loops through standardized contig database file and re-assigns,
    removes, and adds contigs to bins based on CAT/BAT taxonomy annotations
    with bit-scores.
    '''
    contigs_moved = {}
    addThresh = float(addThresh) / 100
    removeThresh = float(removeThresh) / 100
    bin_tax = read_taxonomy(binTaxonFile, bins=True)
    contig_tax = read_taxonomy(contigTaxonFile)
    bin_ident = read_bin_identifier(binidfile)
    remove_count = 0
    add_count = 0
    with open(readme, 'w') as o:
        # Step 1: Loop through all contig: taxonomy pairs
        for contig in contig_tax.keys():
            # Assign a variable for the contig's taxonomy
            contaxon = contig_tax[contig]
            try:
                # Get values for contig's bin and associated bin taxonomy
                mybin = bin_ident[contig]
                # Get taxon for bin based on bin_tax dict
                current_bintaxon = bin_tax[mybin]
            except KeyError:
                mybin = "NoBin"
                current_bintaxon = False
            # If the contig is in a bin, compare with the bin it's been placed in
            if current_bintaxon:
                # Below, return T or F, depending on if contig diverges from bin
                remove = compare_taxonomies(
                    contaxon, current_bintaxon, removeThresh)
                if remove is True:
                    remove_count += 1
                # If contig does not diverge, write it as it was
                else:
                    contigs_moved[contig] = mybin
            else:  # Otherwise, don't compare
                remove = True  # Want to search for a home
            if remove:
                new = compare_new_bin_taxonomies(
                    contaxon, bin_tax, addThresh)
                if new[0] != 0:
                    add_count += 1
                    writeline = f"Contig: {contig}\nOriginal Bin: {mybin}\nContig Taxonomy:\n\t{contaxon}\nOriginal Bin Taxonomy:\n\t{current_bintaxon}\nNew Bin: {new[2]}\nNew Bin score: {new[0]}\nNew Bin Taxonomy:\n\t{new[3]}\n\n"
                    o.write(writeline)
                    contigs_moved[contig] = new[2]
            else:
                pass
    print(f"Removed: {remove_count}\nAdded: {add_count}\n")
    # Finally, write a master binID file:
    with open(output, 'w') as out:
        for contig in contigs_moved.keys():
            out.write(f"{contigs_moved[contig]}\t{contig}\n")
        for contig in bin_ident.keys():
            if contig in contigs_moved:
                pass
            else:
                out.write(f"{bin_ident[contig]}\t{contig}\n")

    return remove_count


if __name__ == '__main__':
    from datetime import datetime
    now = datetime.now()
    date = now.strftime("%m/%d/%Y")
    import argparse
    ''' Init arguments '''
    parser = argparse.ArgumentParser(description='Parser')
    parser.add_argument('-b', '--Bat',
                        help='BAT file from CATBAT', required=True)
    parser.add_argument('-c', '--Cat',
                        help='CAT file from CATBAT',
                        required=True)
    parser.add_argument('-i', '--BinID',
                        help='Bin identification file',
                        required=True)
    parser.add_argument('-m', '--RemoveThreshold',
                        help='Threshold to remove already binned contigs',
                        required=True)
    parser.add_argument('-a', '--AddThreshold',
                        help='Threshold to add already binned contigs',
                        required=True)
    parser.add_argument('-r', '--Readme',
                        help='README file to keep track of steps',
                        required=False, default=f'README.txt')
    parser.add_argument('-o', '--Output', required=True,
                        help='Output file name')
    arg = parser.parse_args()
    loop_contig_taxonomies(arg.Bat, arg.Cat, arg.BinID,
                           arg.RemoveThreshold, arg.AddThreshold,
                           arg.Readme, arg.Output)
