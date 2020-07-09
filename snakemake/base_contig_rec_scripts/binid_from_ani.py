"""
Program that takes in an ANI file and outputs the
node and bin number

Example:
$ python binid_from_ani.py ANIFILE.txt
"""
import sys


def read_ani(anifile):
    ani_dic = {}
    with open(anifile) as a:
        line = a.readline()
        while line:
            sample = line.split('\t')[0]
            node = line.split('\t')[2]
            binname = line.split('\t')[4].strip()
            ani_dic[node] = [sample, binname]
            line = a.readline()
    return ani_dic


def write_binid(anifile):
    ani_dic = read_ani(anifile)
    for node in ani_dic.keys():
        output = str(ani_dic[node][0]) + "-Bins-Updated.txt"
        with open(output, "a") as o:
            outline = f"{node}\t{ani_dic[node][1]}\n"
            o.write(outline)


if __name__ == "__main__":
    write_binid(sys.argv[1])
