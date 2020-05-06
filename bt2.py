from pathlib import Path
import os
import subprocess
import argparse

''' Arguments '''
parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-fq", "--Fastq", help="Fastq file directory",
                    required=True)
parser.add_argument("-i", "--Index", help="Directory of index file",
                    required=True)
parser.add_argument("-o", "--Output", help="Output directory to write \
                    files to", required=True)
argument = parser.parse_args()
''' Arguments '''

''' Load in modules '''
os.system('module load bowtie2/2.3.5.1')
os.system('module load samtools/1.4')


def run_bt2(fq_dir, index, save_directory):
    Path(save_directory).mkdir(parents=True, exist_ok=True)
    os.chdir(save_directory)
    forward, reverse = [], []
    for file in os.listdir(fq_dir):
        if (file.endswith('.fastq') or file.endswith('.fq')):
            if 'R1' in file:
                forward.append(os.path.join(fq_dir, file))
            elif 'R2' in file:
                reverse.append(os.path.join(fq_dir, file))
            else:
                print('Something goofy is going on...')
    forward, reverse = sorted(forward), sorted(reverse)
    # Run alignment
    for f, r in zip(forward, reverse):
        # Go through alignment process for each read-pair.
        samfile = ''.join([f_base.split('.')[0], '_BT2_Def.sam'])
        align = subprocess.Pope((
            'bowtie2', '-x', index,
            '-1', f, '-2', r, '-S', samfile)
        )
        f_base = os.path.basename(f)
        bamfile = ''.join([f_base.split('.')[0], '_BT2_Def.bam'])
        # Save alignment to bamfile
        subprocess.run([
            'samtools', 'view', '-b', samfile, '-o', bamfile],
        )
        align.wait()
        # Proceed to sort & index the BAM file
        bamfile_sorted = ''.join([bamfile.split('.')[0], '_Sorted.bam'])
        subprocess.check_call(
            ['samtools', 'sort', '-o',
             bamfile_sorted,
             bamfile]
        )
        subprocess.check_call(
            ['samtools', 'index',
             bamfile_sorted]
        )


if __name__ == "__main__":
    run_bt2(argument.Fastq, argument.Index, argument.Output)
