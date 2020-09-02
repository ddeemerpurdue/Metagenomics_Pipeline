import sys


def bin_dict(binfile):
    '''
    Function that reads a binID file and returns
    a dictionary
    '''
    bin_dic = {}
    with open(binfile) as i:
        line = i.readline()  # Skip header
        line = i.readline().strip()
        while line:
            line = line.split('\t')
            b = line[0]
            node = line[1]
            bin_dic[node] = b
            line = i.readline().strip()
    return bin_dic


def trim_gff_by_feature(binfile, gff, feat_name):
    '''
    Function that removes all feature entries in
    a .GFF(3) file that do not contain the feature
    indicated
    '''
    bin_dic = bin_dict(binfile)
    new_output = f"{str(gff).rsplit('.', 1)[0]}.{feat_name}.gff"
    with open(new_output, 'w') as o:
        with open(gff) as i:
            line = i.readline()  # Skip header
            o.write(line)
            line = i.readline().strip()
            while line:
                try:
                    cur_line = line.split('\t')
                    node = cur_line[0]
                    features = cur_line[8].split(';')
                    for feature in features:
                        attrib = feature.split('=')[0]
                        if (attrib == str(feat_name) and node in bin_dic):
                            o.write(f"{line}\n")
                            break
                        else:
                            pass
                except IndexError:
                    pass
                line = i.readline().strip()


if __name__ == "__main__":
    trim_gff_by_feature(sys.argv[1], sys.argv[2], sys.argv[3])
