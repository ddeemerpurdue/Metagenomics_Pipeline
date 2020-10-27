'''
HMM Parser, since the .sh one given by HMMER is faulty.
For now, get all output.

Example usage:
$ python hmmparser.py -i hmm-file.dm -o hmm-file.dm.ps
'''
import sys


def read_hmm_ds(file):
    '''
    Read in a HMM file and output a tabular file
    '''
    header = f"Family-HMM\tLength\tQID\tQLen\tE-value\tHMM-Start\tHMM-End\tQStart\tQEnd\tCoverage\n"
    check_dic = {}
    with open(file) as i:
        line = i.readline()
        while line:
            if line.startswith('#'):
                pass
            else:
                # Read in the HMM information
                line = line.split()
                family = line[0]
                hmm_length = line[2]
                qid = line[3]
                qlen = line[5]
                e_val = line[6]
                hmm_start = line[15]
                hmm_end = line[16]
                q_start = line[17]
                q_end = line[18]
                cov = line[21]
                current_bin = str(file).split('.')[2]
                a_len = int(q_end) - int(q_start)
                all_values = [str(a_len), family, hmm_length, qid, qlen, e_val,
                              hmm_start, hmm_end, q_start, q_end, cov, current_bin]
                # Test if redundant
                if qid in check_dic:
                    # Test to see which is better
                    if float(e_val) < float(check_dic[qid][5]):
                        check_dic[qid] = all_values
                    # If they are of equal E-value, grab longest match
                    elif float(e_val) == float(check_dic[qid][5]):
                        if int(a_len) > int(check_dic[qid][0]):
                            check_dic[qid] = all_values
                        else:
                            pass
                else:
                    check_dic[qid] = all_values
            line = i.readline()
    return check_dic


def write_output(file, output):
    '''
    Write the results
    '''
    write_dic = read_hmm_ds(file)
    with open(output, 'w') as o:
        for family in write_dic:
            writeline = '\t'.join(write_dic[family][1:]) + '\n'
            o.write(writeline)
    return 0


write_output(sys.argv[1], sys.argv[2])
