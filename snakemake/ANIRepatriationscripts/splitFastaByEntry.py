'''
Program designed to accomplish 3 tasks:
1. Split a multi-fasta file into many files (1 per entry)
2. Write a list specifying the path to each entry
3. Split that list into N parts and format the names: list_{number}.txt
This is needed for the fastANI portion of the contig repatriation
pipeline.
'''
import os
from random import randrange
from Bio.SeqIO.FastaIO import SimpleFastaParser


def write_entry_files(assembly, outputDirectory, listName, parts, log):
    '''
    Function that takes in a multi-fasta file and outputs a new
    file for each entry with the name being 'defline.fasta'
    Fx also outputs a master list containing paths to each file and
    subsets of that master list (broken in 'parts' parts).
    '''
    # Step 0: If outputDirectory does not exist, create it
    if not os.path.exists(str(outputDirectory)):
        os.makedirs(str(outputDirectory))
    now = time.localtime(time.time())
    all_entries = 0
    entry_count = 0
    split_count = 0
    with open(assembly) as f:
        for values in SimpleFastaParser(f):
            all_entries += 1
            defline = values[0]
            sequence = values[1]
            # Part 1: Write each entry to appropriate file/location
            filename = f"{defline}.fasta"
            full_filename = os.path.join(outputDirectory, filename)
            with open(full_filename, 'w') as o:
                o.write(f">{defline}\n{sequence}\n")
            # Part 2: Append this to the master list
            with open(listName, 'a') as listfile:
                listfile.write(f"{full_filename}\n")
                entry_count += 1
            # Part 3: Write filepaths to split files
            # Loop based on N and write to a different file
            # Randomly write to one file
            rand_number = randrange(1, int(parts) + 1)
            # of all splits
            split_name = f"{listName.rsplit('.', 1)[0]}_{str(rand_number)}.txt"
            with open(split_name, 'a') as splitfile:
                splitfile.write(f"{full_filename}\n")
                split_count += 1
    # Lastly, capture all of the log statistics
    with open(log, 'w') as lg:
        timeline = f"Process started at: {time.asctime(now)}\n"
        lg.write(timeline)
        now = time.localtime(time.time())
        finishtime = f"Process finished at: {time.asctime(now)}\n"
        lg.write(finishtime)
        lg.write(
            f"Total entries found: {all_entries}\nTotal entries written: {entry_count}\n")
        lg.write(f"Total split file entries written: {split_count}\n")
        lg.write(f"Split files written to: {outputDirectory}")
    return 0


if __name__ == "__main__":
    import time
    import argparse
    ''' Init arguments '''
    parser = argparse.ArgumentParser(description='Parser')
    parser.add_argument('-a', '--Assembly',
                        help='Assembly file to parse',
                        required=True)
    parser.add_argument('-o', '--OutputDirectory',
                        help='Directory to write split.fasta files to',
                        required=True)
    parser.add_argument('-l', '--ListName',
                        help='List name to write to (full path)',
                        required=True)
    parser.add_argument('-n', '--Parts',
                        help='How many parts to split the list file into',
                        required=True)
    parser.add_argument('-g', '--Log',
                        help='Log file', required=True)
    arg = parser.parse_args()
    write_entry_files(arg.Assembly, arg.OutputDirectory,
                      arg.ListName, arg.Parts, arg.Log)
