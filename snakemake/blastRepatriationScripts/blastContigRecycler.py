'''
Author: Dane Deemer
Program to recycle contigs back into bins based on blast results using
non-binning contigs as the query and a reference associated with a bin
as the subject. Input for this script should come from running blastn and
be in format 6.
Input files include a blastn results file from running with a bin as query
and the closest assembly as reference, along with a file with non-binners
as query and the same assembly as
Thresholds to consider when repatriating contigs:
Percent identity, contig length, and length of alignment.
Example usage:
$ python contig_recycler.py -b <blastfile.txt> -n <nb_blastfile.txt>
  -t 0.80
'''


def calculate_similarity(cont_length, pident, length):
    '''
    Fx to calculate the level of matching a bin has with it's putative
    (nearest) assembly. This is the percent identity multiplied by the
    ratio of the length of the alignment / contig length.
    '''
    p = pident / 100
    assert cont_length != 0, "Length of contig should not be 0."
    alignment_ratio = length / cont_length
    return p * alignment_ratio


def read_bastn_reference(blast_file):
    '''
    Function that reads in a blastn file and either returns:
    i) The mean - weighted by nucleotide length - percent identity
    for each contig against the reference OR
    ii) A dictionary in the form of dict[binid] = contig AND
    a subsetted blastn file containing information on all repatriated
    contigs.
    '''
    similarities = []
    with open(blast_file) as i:
        line = i.readline().split('\t')
        assert len(line) == 12, "Not a blastn format 6 file"
        while line:
            try:
                contig_length = float(line[0].split('_')[3])
                pident = float(line[2])
                length = float(line[3])
                similarity = calculate_similarity(contig_length,
                                                  pident, length)
                similarities.append(similarity)
            except IndexError:
                "Error with converting column 1, 3 or 4 length into float!"
            line = i.readline().split()
        # Return the mean similarity value
        mean = sum(similarities) / len(similarities)
        return mean


def read_blastn_nobinners(blast_file, nb_blast_file, threshold):
    '''
    Below, standardize names for new output bin identification and
    log files. Prefix will be bin number followed by threshold, whereas
    endings will be .txt or .log. This may change in the future!
    '''
    # Initialize threshold to compare non-binners to
    mean_bin_value = read_bastn_reference(blast_file)
    print(mean_bin_value)
    minimal_match = mean_bin_value - (mean_bin_value * 0.05)
    threshold = float(threshold)

    ''' String name formatting '''
    bin_num = str(nb_blast_file).split('.')[0]
    new_file_prefix = f"{bin_num}.{str(threshold)}."
    binID_file = new_file_prefix + "txt"
    binLog_file = new_file_prefix + "log"
    ''' String name formatting '''

    # If the mean assembly value does not pass the threshold, log and return
    if minimal_match < threshold:
        with open(binLog_file, 'a') as log:
            writeline = f"Minimal match of {minimal_match} is less than {threshold}\n"
            log.write(writeline)
        return 0

    # Otherwise, proceed and check each contig result.
    # Initialize a list of contigs that pass our writing criteria for later:
    contig_write_list = []
    with open(nb_blast_file) as i:
        line = i.readline().split('\t')
        assert len(line) == 12, "Not a blastn format 6 file"
        while line:
            contig = line[0]
            try:
                contig_length = float(line[0].split('_')[3])
                pident = float(line[2])
                length = float(line[3])
                similarity = calculate_similarity(contig_length,
                                                  pident, length)
                if similarity >= minimal_match:
                    '''If we have a match, we want to write this output
                    to a binID file and log the reason. But first, see if we've
                    already encountered and wrote this contig:'''
                    if contig in contig_write_list:
                        pass
                    else:
                        with open(binLog_file, 'a') as log:
                            writeline = '\t'.join(line) + '\n'
                            log.write(writeline)
                        with open(binID_file, 'a') as o:
                            writeline = f"{contig}\t{bin_num}\n"
                            o.write(writeline)
                        contig_write_list.append(contig)
            except IndexError:
                "Error with converting column 1, 3 or 4 length into float!"
                return 1
            line = i.readline().split('\t')
    return 0


if __name__ == '__main__':
    import argparse
    ''' ARGUMENTS '''
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-b", "--Blastfile", help="Blast file for binners.",
                        required=True)
    parser.add_argument("-n", "--NonBinnerBlastFile", help="Blast file for non-binners",
                        required=True)
    parser.add_argument("-t", "--Threshold", help="Percent identity threshold.\
                        Default = 0.80", required=False, default=0.80)
    argument = parser.parse_args()
    ''' ARGUMENTS '''
    read_blastn_nobinners(argument.Blastfile, argument.NonBinnerBlastFile,
                          argument.Threshold)
