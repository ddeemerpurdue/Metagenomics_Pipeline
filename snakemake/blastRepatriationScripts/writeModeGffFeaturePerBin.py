'''
Author: Dane
Note - this can be combined with the gff_mine program in the future!!!!!
Script designed to take in a file created from the gff_mine.py script and produce
the top feature for each bin. This script takes in as input the output from gff_mine.py
with the option --Top specified.
Example usage:
$ python writeTopGFFFeaturePerBin.py <inputfeature.txt> <output.txt>
'''
import time


def read_attribute_output(attfile, logfile):
    '''
    Read in GFF attribute output and return a dictionary
    in the form of dict[bin] = attribute, where attribute is
    the attribute occuring the most for the bin
    '''
    output_dic = {}
    final_dic = {}
    adequate = False
    '''
    Output dictionary is in the form of:
    dict[bins] = [attr1, attr2, etc.]...
    final_dic just subsets for highest occuring attribute
    '''
    with open(attfile) as f:
        line = f.readline()  # Skip header
        line = f.readline()
        while line:
            line = line.split('\t')
            assert len(line) == 4, "Error in number of delimiters."
            attrib = line[1]
            bins = line[3].strip()
            try:
                output_dic[bins].append(attrib)
            except KeyError:
                output_dic[bins] = [attrib]
            line = f.readline()
    # Now loop through dictionary we created and log
    with open(logfile, 'a') as lg:
        lg.write(
            f"Full tab-delimited record of all bins with assemblies with >=50% contig consensus:\n")
        for bins in output_dic.keys():
                                                # Variable to get how many contigs in bin
            bin_size = len(output_dic[bins])
            threshold = bin_size * 0.5          # Must be 50% consensus or more
            for assembly in set(output_dic[bins]):
                                                # Count each occurence of attrib
                count = output_dic[bins].count(assembly)
                if count > threshold:
                    lg.write(f"{bins}\t{assembly}\t{count}\n")
                    final_dic[bins] = [assembly, count]
                    adequate = True             # If flagged, don't write log below
            if not adequate:
                lg.write(
                    f"Bin {bins} does not have an adequate assembly match.\n")
    return final_dic


def write_new_bin_file(attfile, output, logfile):
    '''
    Function that takes in a dictionary in the form of:
    dict[bin] = [Assembly, N]
    '''
    # Dict[node] = bin
    # bin_dic = get_bin_dictionary(binfile)
    # Dict[bin] = Attribute
    att_dict = read_attribute_output(attfile, logfile)
    with open(output, 'w') as o:
        header = f"Bin\tAssembly\tN"
        for key in att_dict.keys():
            vals = '\t'.join([str(val) for val in att_dict[key]])
            writeline = key + '\t' + vals + '\n'
            o.write(writeline)
    return 0


if __name__ == "__main__":
    now = time.localtime(time.time())
    import argparse
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-a", "--AttributeFile",
                        help="Attribute file", required=True)
    parser.add_argument("-o", "--Output",
                        help="Output file to write to", required=True)
    parser.add_argument("-l", "--Log",
                        help="Log file to write to", required=True)
    argument = parser.parse_args()
    with open(argument.Log, 'a') as lg:
        lg.write(f"Time started: {time.asctime(now)}\n")
        write_new_bin_file(argument.AttributeFile, argument.Output,
                           argument.Log)
        now = time.localtime(time.time())
        lg.write(f"Time ended: {time.asctime(now)}\n")
