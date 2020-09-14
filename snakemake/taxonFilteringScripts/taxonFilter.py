'''
Author: Dane
Remove contamination from a metagenomic analysis by
comparing bit-scores across contig annotations vs.
full bin annotations. Scores are obtained from a source
such as CATBAT.
Example usage:
$ python taxonFilter.py -b binFile.txt -c masterDB.txt -m 0.70
-a 0.50 -o outputFile.txt
'''


def read_bin_taxonomies(binCategoricalFile):
    '''
    Function designed to read in a categorical bin file and return a
    dictionary containing the bin number and taxonomy.
    '''
    bin_dic = {}
    with open(binCategoricalFile) as i:
        line = i.readline()  # Skip header
        line = i.readline()
        while line:
            line = line.split('\t')
            binName = line[0]
            binTaxon = line[1]
            bin_dic[binName] = binTaxon
            line = i.readline()
    return bin_dic


def compare_taxonomies(contigTaxon, binTaxon, threshold):
    '''
    Function designed to compare taxonomies. The return value (True, False)
    can be used by another function to remove the contig from the bin and
    look for another bin to place the contig into.
    '''
    remove = False
    threshold = float(threshold)
    contigTax = contigTaxon.split(',')
    binTax = binTaxon.split(',')
    matches = 0
    total = 0
    '''
    Part A takes care of either the bin or the contig having higher levels
    of taxonomic annotations because zipping will only compare the first n
    values, where n is the length of the shorter list.
    '''
    if binTaxon == 'NoTaxon':
        return remove
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
                remove = True
            else:
                pass
        total += 1
    return remove


def compare_new_bin_taxonomies(contigTaxon, bin_dic, threshold):
    '''
    Program that looks for a new home for a removed contig based
    on bit-scores. Return value is a list containing:
    contig repatriation score, bin confidence, and bin
    '''
    contigTax = contigTaxon.split(',')
    test_container = [0.0, 0.0, 'NoTaxon', 'NoTaxon']
    for binID in sorted(list(bin_dic.keys())):
        # Init counting variables
        bin_confidence = 0
        matches = 0.0
        score = 0.0
        binTax = bin_dic[binID].split(',')
        # Loop through taxonomies
        for cval, bval in zip(contigTax, binTax):
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
                    score += (bitscore * float(matches))
                    bin_confidence += (bitscore * float(matches))
                else:  # Below, test if discrepancy too big to repair
                    # If both annotations are high level and don't match, then don't
                    # resolve this contig and make score = 0 so they don't collapse.
                    if (bitscore > threshold and bin_bitscore > threshold):
                        score = 0
                        break
                    # Otherwise, keep score where it is and potentially we can
                    # add this information
                    else:
                        pass
        # Test to see if we have a potential re-patriation candidate
        if score > 10.0:  # This accounts for first 4 levels of taxonomy
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


def loop_contig_taxonomies(binCategoricalFile, contigFile, removeThresh, addThresh, readme, output):
    '''
    Fx that loops through standardized contig database file and re-assigns,
    removes, and adds contigs to bins based on CAT/BAT taxonomy annotations
    with bit-scores.
    '''
    bin_dic = read_bin_taxonomies(binCategoricalFile)
    count = 0
    with open(readme, 'w') as o:
        with open(output, 'w') as out:
            with open(contigFile) as i:
                line = i.readline()  # Skip header
                out.write(line)
                line = i.readline()
                while line:
                    try:
                        line = line.split('\t')
                        contig = line[0]
                        mybin = line[70]
                        contaxon = line[71]
                        bintaxon = line[73]
                        if (contaxon == 'NoTaxon' or contaxon == ''):
                            out.write('\t'.join(line))
                        else:
                            remove = compare_taxonomies(
                                contaxon, bintaxon, removeThresh)
                            if remove:
                                count += 1
                                new = compare_new_bin_taxonomies(
                                    contaxon, bin_dic, addThresh)
                                writeline = f"Contig: {contig}\nOriginal Bin: {mybin}\nContig Taxonomy:\n\t{contaxon}\nOriginal Bin Taxonomy:\n\t{bintaxon}\nNew Bin: {new[2]}\nNew Bin score: {new[0]}\nNew Bin Taxonomy:\n\t{new[3]}\n\n"
                                line[70] = new[2]
                                if line[73] == 'NoTaxon':
                                    pass
                                else:
                                    line[73] = ','.join(new[3])
                                o.write(writeline)
                            else:
                                pass
                            out.write('\t'.join(line))
                    except IndexError:
                        pass
                    line = i.readline()
    print(count)
    return count


if __name__ == '__main__':
    from datetime import datetime
    now = datetime.now()
    date = now.strftime("%m/%d/%Y")
    import argparse
    ''' Init arguments '''
    parser = argparse.ArgumentParser(description='Parser')
    parser.add_argument('-b', '--BinFile',
                        help='Categorical bin file', required=True)
    parser.add_argument('-c', '--ContigFile',
                        help='Master database file created from annotationFormmatting.py',
                        required=True)
    parser.add_argument('-m', '--RemoveThreshold',
                        help='Threshold to remove already binned contigs',
                        required=True)
    parser.add_argument('-a', '--AddThreshold',
                        help='Threshold to add already binned contigs',
                        required=True)
    parser.add_argument('-r', '--Readme',
                        help='README file to keep track of steps',
                        required=False, default=f'README-{date}.txt')
    parser.add_argument('o', '--Output', required=True,
                        help='Output file name')
    arg = parser.parse_args()
    loop_contig_taxonomies(arg.BinFile, arg.ContigFile,
                           arg.RemoveThreshold, arg.AddThreshold,
                           arg.Readme, arg.Output)
