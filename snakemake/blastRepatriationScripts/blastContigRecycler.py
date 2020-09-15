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
import os


def calc_similarity(cont_length, pident, length):
    '''
    Fx to calculate the level of matching a bin has with it's putative
    (nearest) assembly. This is the percent identity multiplied by the
    ratio of the length of the alignment / contig length.
    '''
    p = pident / 100
    assert cont_length != 0, "Length of contig should not be 0."
    # alignment_ratio = length / cont_length
    # return p * alignment_ratio
    return p


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
    total_length = 0
    with open(blast_file) as i:
        line = i.readline().split('\t')
        assert len(line) == 12, "Not a blastn format 6 file"
        while line:
            try:
                contig_length = float(line[0].split('_')[3])
                pident = float(line[2])
                length = int(line[3])
                total_length += length
                similarity = calc_similarity(contig_length,
                                             pident, length)
                similarities.append(similarity)
            except IndexError:
                "Error with converting column 1, 3 or 4 length into float!"
            line = i.readline().split()
        # Return the mean similarity value
        mean = sum(similarities) / len(similarities)
        return mean, total_length


def format_name(nb_blast_file, threshold):
    ''' String name formatting '''
    nb_file = os.path.basename(str(nb_blast_file))
    bin_num = str(nb_file).split('.')[0]
    new_file_prefix = f"{bin_num}.{str(threshold)}."
    binID_file = new_file_prefix + "txt"
    binLog_file = new_file_prefix + "log"
    return binID_file, binLog_file, bin_num
    ''' String name formatting '''


def write_blastn_nobinners(blast_file, nb_blast_file, threshold):
    '''
    Below, standardize names for new output bin identification and
    log files. Prefix will be bin number followed by threshold, whereas
    endings will be .txt or .log. This may change in the future!
    '''
    # Initialize threshold to compare non-binners to
    mean_bin_value, total_length = read_bastn_reference(blast_file)
    minimal_match = mean_bin_value - (mean_bin_value * 0.05)
    minimal_match_length = 10000
    print(f"Minimal match: {minimal_match}")
    print(f"Total alignment length: {total_length}")
    threshold = float(threshold)
    print(f"Threshold: {threshold}")

    ''' String name formatting '''
    binID_file, binLog_file, bin_num = format_name(
        nb_blast_file, str(threshold).split('.')[1])
    ''' String name formatting '''

    # If the mean assembly value does not pass the threshold, log and return
    if (minimal_match < threshold or total_length < minimal_match_length):
        print('Bin does not match assembly above the listed threshold.')
        with open(binLog_file, 'a') as log:
            writeline = f"Minimal match of {minimal_match} is less than {threshold}\n"
            log.write(writeline)
        return {}

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
                similarity = calc_similarity(contig_length,
                                             pident, length)
                if (similarity >= minimal_match and length > 1500):
                    '''If we have a match, we want to write this output
                    to a binID file and log the reason. But first, see if we've
                    already encountered and wrote this contig:'''
                    if contig in contig_write_list:
                        pass
                    else:
                        with open(binLog_file, 'a') as log:
                            writeline = '\t'.join(line)
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


def read_blastn_nobinners(blast_file, nb_blast_file, threshold):
    '''
    Below, standardize names for new output bin identification and
    log files. Prefix will be bin number followed by threshold, whereas
    endings will be .txt or .log. This may change in the future!
    '''
    # First, test if file is empty.
    if os.stat(blast_file).st_size == 0:
        print(f"File - {blast_file} -  is empty!")
        return {}, []
    # Initialize threshold to compare non-binners to
    mean_bin_value, total_length = read_bastn_reference(blast_file)
    minimal_match = mean_bin_value - (mean_bin_value * 0.05)
    threshold = float(threshold) / 100

    ''' String name formatting '''
    binID_file, binLog_file, bin_num = format_name(nb_blast_file, threshold)
    ''' String name formatting '''
    # If the mean assembly value does not pass the threshold, log and return
    if minimal_match < threshold:
        print('Bin does not match assembly above the listed threshold.')
        return {}, []

    log_information = []
    loop_dict = {}  # Dictionary in the form: dict[contig]=[Bin, PID, len]
    # to keep track of best contig matches
    with open(nb_blast_file) as i:
        line = i.readline().split('\t')
        assert len(line) == 12, "Not a blastn format 6 file"
        while line != ['']:
            contig = line[0]
            try:
                contig_length = float(line[0].split('_')[3])
                pident = float(line[2])
                length = float(line[3])
                similarity = calc_similarity(contig_length,
                                             pident, length)
                if (similarity >= minimal_match and length > 1500):
                    '''If we have a match, we want to write this output
                    to a binID file and log the reason. But first, see if we've
                    already encountered and wrote this contig:'''
                    if contig in loop_dict:
                        if (pident * length) > (loop_dict[contig][1] * loop_dict[contig][2]):
                            loop_dict[contig] = [bin_num, length, pident]
                            log_information.append(line)
                        else:
                            pass
                    else:
                        loop_dict[contig] = [bin_num, length, pident]
                        log_information.append(line)
            except IndexError:
                "Error with converting column 1, 3 or 4 length into float!"
                return 1
            line = i.readline().split('\t')
    return loop_dict, log_information


def write_all_blastn(blast_files, nb_blast_files, threshold):
    '''
    Function that loops through a list of blast files and takes care of contigs
    that pass the threshold in multiple bins by placing contigs into the best
    hit bin
    '''
    master_log = []
    threshold = float(threshold) / 100
    final_blast_dic = {}
    blast_files = sorted(blast_files)
    nb_blast_files = sorted(nb_blast_files)
    # Loop through all files
    for blast_file, nb_blast_file in zip(blast_files, nb_blast_files):
        print(
            f"BF: {str(os.path.basename(blast_file))}\nNBF: {str(os.path.basename(nb_blast_file))}")
        # Make sure both files are the same
        assert str(os.path.basename(blast_file)) == str(
            os.path.basename(nb_blast_file))
        current_dic, current_log = read_blastn_nobinners(
            blast_file, nb_blast_file, threshold)
        master_log.extend(current_log)  # Add log information to list
        for contig in current_dic.keys():
            if contig in final_blast_dic:
                # Compare the pid*len values of each and keep the biggger
                current_stat = current_dic[contig][1] * current_dic[contig][2]
                _stat = final_blast_dic[contig][1] * final_blast_dic[contig][2]
                if _stat >= current_stat:
                    pass
                else:
                    final_blast_dic[contig] = current_dic[contig]
            else:
                final_blast_dic[contig] = current_dic[contig]

    # Now write the results to an output file:
    output = f"FinalBlastNResults.T{str(threshold).split('.')[1]}.txt"
    log_output = f"FinalBlastNResults.T{str(threshold).split('.')[1]}.log"
    with open(output, 'w') as o:
        for key, value in final_blast_dic.items():
            if key == 'NoTaxon':
                pass
            else:
                writeline = value[0] + '\t' + key + '\n'
                o.write(writeline)
    with open(log_output, 'w') as log:
        for line in master_log:
            writeline = '\t'.join([str(val) for val in line])
            log.write(writeline)
    return 0


if __name__ == '__main__':
    import argparse
    ''' ARGUMENTS '''
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-b", "--Blastfile", help="Blast file for binners.",
                        required=True, nargs='*')
    parser.add_argument("-n", "--NonBinnerBlastFile", help="Blast file for non-binners",
                        required=True, nargs='*')
    parser.add_argument("-t", "--Threshold", help="Percent identity threshold.\
                        Default = 80", required=False, default=80)
    argument = parser.parse_args()
    ''' ARGUMENTS '''
    if len(argument.Blastfile) == 1:
        blastfile = argument.Blastfile[0]
        nbblastfile = argument.NonBinnerBlastFile[0]
        write_blastn_nobinners(blastfile, nbblastfile, argument.Threshold)
    else:
        write_all_blastn(argument.Blastfile,
                         argument.NonBinnerBlastFile, argument.Threshold)
