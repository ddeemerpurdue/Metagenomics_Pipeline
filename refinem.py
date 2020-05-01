import os
import sys
import subprocess
print(os.getcwd())
from maxbin2 import read_sample_file
import argparse


''' SAMPLE TEXT '''
# If this is being called on its own, we need it to be loading int
# the necessary modules:
os.system('module load bioinfo')
os.system('module load RefineM/0.0.25')
os.system('module load blast/2.9.0+')
os.system('module load HMMER/3.2.1')


assembly, _, _ = read_sample_file(sys.argv[1])
mb_output = '/scratch/snyder/d/ddeemer/WhiteRed/MaxBin_Results/BWA_All/Bins'
bamfiles = '/scratch/snyder/d/ddeemer/WhiteRed/BAM/W1_Aligned.sorted.bam'
basedir = os.chdir('/scratch/snyder/d/ddeemer/RM_Test')
# Runthrough RefineM
subprocess.check_call(['refinem', 'scaffold_stats', '-c', '20',
                       '-x', 'fasta', assembly, mb_output,
                       'stats_output_dir', bamfiles])