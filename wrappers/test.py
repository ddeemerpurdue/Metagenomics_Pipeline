import argparse

parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-f", "--File", help="File with paths for MaxBin input",
                    required=True)
parser.add_argument("-r", "--Refinem", help="Tell the program to proceed with \
                    refining the bins using RefineM after completion",
                    required=False, action='store_true')
parser.add_argument("-b", "--Bams", help="This argument is required \
		                if the --Refinem flag is given!", required=False)

argument = parser.parse_args()

if __name__ == "__main__":
	if argument.Refinem and (argument.Bams is None):
		parser.error('When --Refinem required --Bams please!')
	print(argument.File)
	if argument.Refinem:
		print(argument.Refinem)
		print(argument.Bams)