'''
Author: Dane
Script designed to convert raw counts (e.g., output from htseq-count) into CPMs.
If the contig node is specified vs. a particular feature, then the script will
have to be ran slightly different.
Example usage:
1.) When count file has contig node as first field and count as second:
python convertCPMs.py -g <gff.gff> -f 0 -c <htseqCountfile.txt>
2.) When count file has a feature as field field and count as second:
python convertCPMs.py -g <gff.gff> -f 8 -c <htseqCountfile.txt> -n tigrfam_acc
Note: field 8 of the gff file is where the feature attributes are located, so
here we specified we are looking to count a feature (-f 8) and the specific
type of feature is tigrfam_acc (-n tigrfam_acc)
'''


def get_GFF_length_dic(gff, featColumn, featName=None):
    '''
    Read in a GFF file and output a dictionary in the format:
    dict[feature]: length
    For contigs, this will be total contig length. For a specific
    feature predicted by prodigal/an annotator, this will be the CDS
    region length.
    '''
    length_dic = {}
    with open(gff) as g:
        line = g.readline()  # Skip header
        line = g.readline()
        while line:
            line = line.split('\t')
            feature = line[int(featColumn)]
            if featName is not None:
                values = feature.split(';')
                for value in values:
                    if str(featName) == value.split('=')[0]:
                        feature = value.split('=')[1]
                        break
                    else:
                        pass
            length = int(line[4]) - int(line[3])
            length_dic[feature] = length
            line = g.readline()
    print(length_dic)
    return length_dic


def calcRPK(raw, length):
    '''
    Convert raw read counts into reads per kilobase.
    Length is the length of the feature.
    '''
    raw = int(raw)
    length = int(length)
    kblength = length / 1000
    rpk = raw / kblength
    return rpk


def process_htseqcount_output(gffFile, featColumn, htseqFile,
                              featName=None):
    '''
    Function that processes htseq-count output and returns
    CPM values
    '''
    length_dic = get_GFF_length_dic(gffFile, featColumn, featName)
    rpk_dic = {}
    # Part 0: Read in the file
    output = f"{str(htseqFile).rsplit('.', 1)[0]}.processed.cpm.txt"
    with open(htseqFile) as i:
        line = i.readline()
        while line:
            cur_line = line.split('\t')
            feature = cur_line[0]
            raw_count = cur_line[1]
            # First, check if a CDS region multi-mapped to many features:
            if len(feature.split(',')) > 1:
                pass  # Do not want to include ambiguous feature annotations
            else:
                try:
                    length = length_dic[feature]
                    rpk = calcRPK(raw_count, length)
                    rpk_dic[feature] = rpk
                except KeyError:
                    pass
            line = i.readline()
    # Loop and calculate CPM
    total = sum(rpk_dic.values()) / 1000000
    tpm = {k: v / total for k, v in rpk_dic.items()}
    with open(output, 'w') as o:
        o.write(f"Feature\tTPM\n")
        for key, value in tpm.items():
            o.write(f"{key}\t{value}\n")
    return 0


if __name__ == '__main__':
    import argparse
    ''' Arguments '''
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-g", "--GFF", required=True,
                        help="GFF file used for counting.")
    parser.add_argument("-f", "--FeatureCol", help="Feature column",
                        required=True)
    parser.add_argument("-c", "--CountFile", help="Count file output. 2 fields (tab-delimited)",
                        required=True)
    parser.add_argument("-n", "--FeatureName", help="Feature column",
                        required=False, default=None)
    argument = parser.parse_args()
    ''' Arguments '''
    process_htseqcount_output(argument.GFF, argument.FeatureCol,
                              argument.CountFile,
                              featName=argument.FeatureName)
