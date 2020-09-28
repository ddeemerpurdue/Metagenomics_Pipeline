'''
Program to filter .FASTA files based on sequences
that are at least N nucleotides long
Example usage:
$ python filterSeqlength.py <input.fasta> <length_threshold> <output.fasta>
'''

from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser


def filter_fasta(file, threshold, output, log):
    """
    Open up a .fasta file and write only entries with
    sequences >= a specified length.
    """
    now = time.localtime(time.time())
    kept = 0
    removed = 0
    with open(output, 'w') as o:
        with open(file) as f:
            for values in SimpleFastaParser(f):
                defline = values[0]
                length = len(values[1])
                if length >= int(threshold):
                    o.write(f">{defline}\n{values[1]}\n")
                    kept += 1
                else:
                    removed += 1
    with open(log, 'w') as o:
        timeline = f"Time started: {time.asctime(now)}\n"
        o.write(timeline)
        keptline = f"{kept} entries were kept.\n"
        o.write(keptline)
        removedline = f"{removed} entries were removed.\n"
        o.write(removedline)
    return 0


if __name__ == "__main__":
    import time
    import argparse
    ''' Init arguments '''
    parser = argparse.ArgumentParser(description='Parser')
    parser.add_argument('-a', '--Assembly',
                        help='1 or more assembly files to filter',
                        required=True)
    parser.add_argument('-l', '--Length',
                        help='Length to filter the assembly by',
                        required=True)
    parser.add_argument('-o', '--Output',
                        help='Output file to write new assembly to',
                        required=True)
    parser.add_argument('-g', '--Log',
                        help='Log file to keep track of what was written',
                        required=True)
    arg = parser.parse_args()
    assert str(arg.Assembly).endswith(
        ('.fasta', '.fa')), "Invalid assembly file extension!"
    filter_fasta(arg.Assembly, arg.Length, arg.Output, arg.Log)
