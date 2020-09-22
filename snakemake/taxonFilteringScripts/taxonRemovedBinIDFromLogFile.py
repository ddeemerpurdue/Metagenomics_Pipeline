'''
Scripts that takes in a log file from taxonFiltering.py
and write a bin identification file of the unique contigs
moved during this processing stage.

Example usage:
$ python taxRemovedBinID.py taxonLogFile.log output.txt
'''


def write_bin_id(logfile, outfile):
    current_contig, current_bin = '', ''
    bin_dic = {}
    with open(logfile) as i:
        lines = i.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith('Contig:'):
                current_contig = line.split(':')[1]
            elif line.startswith('New Bin: '):
                current_bin = line.split('New Bin: ')[1]
                bin_dic[current_contig] = current_bin
                current_contig, current_bin = '', ''
            else:
                pass
    with open(outfile, 'w') as o:
        for contig in bin_dic.keys():
            writeline = f"{bin_dic[contig]}\t{contig}\n"
            o.write(writeline)


if __name__ == "__main__":
    import sys
    write_bin_id(sys.argv[1], sys.argv[2])
