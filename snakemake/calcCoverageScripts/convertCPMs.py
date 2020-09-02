import sys


def get_GFF_length_dic(gff, featColumn, featName=None):
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
    process_htseqcount_output(sys.argv[1], sys.argv[2], sys.argv[3],
                              featName=sys.argv[4])
