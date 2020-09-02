'''
Author: Dane
Script to calculate the variance between forward and reverse
reads given htseq-count output for each file.
If the -o/--Output flag is given, the mean of F/R will be
computed and a new count file will be created with the mean
count.
Example usage:
$ python calcFRVariance.py -f <forward.txt> -r <reverse.txt> \
-o <optionally-write-output.txt>
'''


def calculate_FR_deviation(forward, reverse, output=False):
    '''
    Fx that reads in two files, calculated deviation and,
    if specified, writes a new mean output file.
    '''
    n = 0
    mean_values = {}
    mean_variances = []
    if output:
        with open(output, 'a') as o:
            o.write(f"Feature\tMeanCPM\n")
    with open(forward) as f:
        with open(reverse) as r:
            f_line, r_line = f.readline(), r.readline()
            f_line, r_line = f.readline().strip(), r.readline().strip()
            while (f_line and r_line):
                n += 1
                f_line = f.readline().split('\t')
                r_line = r.readline().split('\t')
                f_feat, r_feat = f_line[0], r_line[0]
                try:
                    f_val, r_val = float(f_line[1]), float(r_line[1])
                except IndexError:
                    break
                if f_feat != r_feat:
                    return f"Values don't match on line {n}."
                else:
                    mean = (f_val + r_val) / 2.0
                    var = (((f_val - mean) ** 2.0) +
                           ((r_val - mean) ** 2.0)) / 2.0
                    mean_variances.append(var)
                    if output:
                        with open(output, 'a') as o:
                            o.write(f"{f_feat}\t{mean}\n")
                f_line, r_line = f.readline().strip(), r.readline().strip()
    mean_variance = (sum(mean_variances) / len(mean_variances))
    if output:
        with open(output, 'a') as o:
            o.write(f"Mean_Variance\t{mean_variance}")
    print(f"The mean variance is: {mean_variance}")
    return 0


if __name__ == "__main__":
    import argparse
    ''' Arguments '''
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-f", "--Forward", required=True,
                        help="Forward file")
    parser.add_argument("-r", "--Reverse", required=True,
                        help="Reverse file")
    parser.add_argument("-o", "--Output", required=False,
                        help="Output file name")
    argument = parser.parse_args()
    ''' Arguments '''
    if argument.Output:
        calculate_FR_deviation(argument.Forward, argument.Reverse,
                               argument.Output)
    else:
        calculate_FR_deviation(argument.Forward, argument.Reverse)

