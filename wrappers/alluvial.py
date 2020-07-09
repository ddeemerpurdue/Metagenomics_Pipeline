#!/usr/bin/env python

# lousy script that works with the outputs of `anvi-export-collection`
# to reconcstruct the fate of contigs for a given algorithm and bin and
# spit out some text to be visualized on https://app.rawgraphs.io/

import sys

from collections import OrderedDict

import anvio.utils as u

from anvio.errors import ConfigError

G = lambda x: u.get_TAB_delimited_file_as_dictionary(x, no_header=True, column_names=["split_name", "bin_name"])

# crappy way to do this, indeed, but it will suffice for today.
# you need to put the output of anvi-export-collection for each
# binning algorithm into this dict with the matching filename:
algorithms = OrderedDict({'All_Def': G('collection-Red1_BT2Default_All_Bins.txt'),
                          'All_99F': G('collection-Red1_BT2Default_All99F_Bins.txt'),
                          'R1_99F': G('collection-Red1_BT2Default_R199F_Bins.txt')})

def main(args):
    algorithm = args.algorithm
    bin_name = args.bin

    if algorithm not in algorithms:
        raise ConfigError("That algorithm we don't know. This is what we know: %s." % ', '.join(algorithms))
        sys.exit(-1)

    split_names = set()
    for split_name in algorithms[algorithm]:
        if algorithms[algorithm][split_name]['bin_name'] == bin_name:
            split_names.add(split_name)

    print(','.join(['split'] + list(algorithms.keys())))
    for split_name in split_names:
        line = [split_name]
        for algorithm in algorithms:
            try:
                line.append(algorithms[algorithm][split_name]['bin_name'])
            except KeyError:
                line.append('No_Bin')
        print(','.join([val for val in line]))
        # try:
        #     print(','.join([split_name] + [algorithms[algorithm][split_name]['bin_name'] for algorithm in algorithms]))
        # except KeyError:
        #     print('Error')


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser("Give this guy one of the algorithms (it knows all about '%s') and a bin name,\
                                      get back an input data for an alluvial diagram that shows which bins contain\
                                      the splits in that bin..")

    parser.add_argument('--algorithm', default=None, help="Which algorithm?")
    parser.add_argument('--bin', default=None, help="Which bin name to start with?")

    args = parser.parse_args()

    try:
        main(args)
    except ConfigError as e:
        print(e)
        sys.exit(-1)