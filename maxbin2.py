import os
import subprocess
import argparse

# file = sys.argv[1]
""" So if we import this into another module, you still
want to get the information from the file, but not execute
the code """

''' Create arguments for MaxBin, that way a flag can be raised to tell
the program to proceed with refineM '''
parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-f", "--File", help="File with paths for MaxBin input",
                    required=True)
parser.add_argument("-r", "--Refinem", help="Tell the program to proceed with \
                    refining the bins using RefineM after completion",
                    required=False, action='store_true')
parser.add_argument("-b", "--Bams", help="This argument is required \
                    if the --Refinem flag is given!", required=False,
                    nargs='*')
parser.add_argument("-d", "--Basedir", help="Base directory to save \
                    refinement files to", required=False)
parser.add_argument("-db", "--Reference", help="Path to the directory \
                    containing the database information required for refinem",
                    required=False)

argument = parser.parse_args()


def read_sample_file(file):
    ''' Increase functionality to allow multiple lines be read!!! '''
    with open(file) as f:
        line = f.readline().split('\t')
        assert line[0].lower() == 'assembly', 'Make sure first column is \
                                               assembly!'
        assert line[1].lower() == 'abundance', 'Make sure second column \
                                                is abundance!'
        #assert line[2].lower() == 'output', 'Make sure third column is output'
        line = f.readline().split('\t')
        assembly = line[0].strip()
        abundance = line[1].strip()
        output = line[2].strip()
    return assembly, abundance, output


if __name__ == "__main__":
    os.system('module load bioinfo')
    os.system('module load MaxBin/2.2.3')
    assembly, abundance, output = read_sample_file(argument.File)
    subprocess.check_call(['run_MaxBin.pl', '-contig', assembly, '-abund_list',
                           abundance, '-out', output])
    if argument.Refinem and (argument.Bams is None
                             or argument.Basedir is None
                             or argument.Reference is None):
        parser.error('When --Refinem is used, --Bams, --Basedir and \
                      --Reference are all required as well.')
    elif argument.Refinem and (argument.Bams and argument.Basedir
                               and argument.Reference):
        from refinem import run_refinem
        if os.path.exists(argument.Basedir):
            pass
        else:
            os.mkdir(argument.Basedir)
        try:
            bins = output.split('/')[-2]
        except IndexError:
            bins = os.getcwd()
        ref = argument.Reference
        print(ref)
        run_refinem(assembly, output, argument.Bams, argument.Basedir, ref)
    else:
        pass
