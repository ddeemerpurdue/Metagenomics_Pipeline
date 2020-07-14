#!/usr/bin/env python3
import io
import argparse
import pandas as pd
import os


def dict_from_file(d, filename):
    with open(filename) as f:
        sample = filename[48:50]
        line = f.readline()  # Skip header
        for line in f:
            row = line.split()
            bin = sample + "_" + row[2]
            if row[0] in d:
                d[row[0]].append(bin)
            else:
                d[row[0]] = [bin]
        return d


if __name__ == "__main__":
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-f", "--Files", help="Txt file which contains the words:",
                        required=True)
    parser.add_argument("-d", "--Directory", required=True)
    argument = parser.parse_args()
    master = {}
    for filename in os.listdir(argument.Directory):
        if argument.Files in filename:
            dict_from_file(master, argument.Directory+filename)
    f = open("top_acc_merged.txt", "w")
    for key in master.keys():
        f.write(key + ":")
        f.write(','.join(map(str, master.get(key))))
        f.write("\n")
    f.close()