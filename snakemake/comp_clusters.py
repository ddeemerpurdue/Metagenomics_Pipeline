#!/usr/bin/env python3
import io
import argparse
import difflib


def files_to_df(gc, ac):
    genome = []
    ani = []
    with open(gc) as f:
        for line in f:
            line = line.rstrip().split(":",1)[1]
            if "," in line:
                genome.append(line)

    with open(ac) as f:
        for line in f:
            line = line.split("\"")[1]
            ani.append(line)

    comp(genome, ani)

def comp(gc, ac):
    f = open("comp.txt", "wt")
    f.write("genome cluster"+"\t"+"ani cluster"+"\t"+"common cluster"+"\n")
    for g in gc:
        glist= g.split(",")
        found = False
        for a in ac:
            alist = a.split(",")
            alist_set = set(alist)
            intersection = alist_set.intersection(glist)
            intersection_as_list = list(intersection)
            alist = set(alist) - set(intersection_as_list)
            glist = set(glist) - set(intersection_as_list)
            intersection = ','.join(intersection_as_list)
            a = ','.join(alist)
            g = ','.join(glist)
            if intersection_as_list:
                found = True
                if a == "":
                    a = "none"
                if g == "":
                    g = "none"
                f.write(g + "\t" + a + "\t" + intersection + "\n")
        if (g != "\t") & (~found):
            f.write(g + "\t" + "none" + "\t" + "none" + "\n")
    f.close()


if __name__ == "__main__":
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-gc", "--Genome", help="Txt file with genome clusters",
                        required=True)
    parser.add_argument("-ac", "--Ani", help="Txt file with ani clusters",
                        required=True)
    argument = parser.parse_args()
    files_to_df(argument.Genome, argument.Ani)
