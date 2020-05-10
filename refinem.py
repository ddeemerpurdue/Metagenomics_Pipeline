import os
import sys
import subprocess


'''
Python wrapper to call RefineM through.

Variables needed:
Scaffold (From Maxbin), BINs, WDIR, BAMs, Reference, Base_Directory

'''
# If this is being called on its own, we need it to be loading int
# the necessary modules:

os.system('module load RefineM/0.0.25')
os.system('module load blast/2.9.0+')
os.system('module load HMMER/3.2.1')


# Read in information file
def read_refinem_info_file(file):
    with open(file) as f:
        # Read in header information
        line = f.readline().split('\t')
        assert line[0].lower() == 'assembly', "Make sure first column is \
            'assembly'"
        assert line[1].lower() == 'bins', "Make sure second column is 'bin'"
        assert line[2].lower() == 'bam', "Make sure third column is 'bam'"
        assert line[3].lower() == 'base', "Make sure the fourth column is \
            'base'"
        assert line[4].lower() == 'reference', "Make sure the fifth column \
            is 'reference"
        # Read in information now
        line = f.readline().split('\t')
        assembly = line[0]
        bins = line[1]
        bamfiles = line[2]
        basedir = line[3]
        reference = line[4]
        return assembly, bins, bamfiles, basedir, reference


# Runthrough RefineM
def run_refinem(assembly, bins, bams, basedir, reference):
    # Remove based on genomic properties
    rf = []
    rf.extend(['refinem', 'scaffold_stats', '-c', '20',
              '-x', 'fasta', assembly, bins, 'stats_output_dir'])
    rf.append(bams)
    subprocess.check_call(rf)
    sc_stats = os.path.join(basedir, 'stats_output_dir', 'scaffold_stats.tsv')
    subprocess.check_call(['refinem', 'outliers', '--no_plots', sc_stats,
                           'outlier_output'])
    outliers = os.path.join(basedir, 'outlier_output', 'outliers.tsv')
    subprocess.check_call(['refinem', 'filter_bins', '-x', 'fasta',
                           bins, outliers, 'genomic_filtered_output'])
    # Remove based on taxnomonic properties
    subprocess.check_call(['refinem', 'call_genes', '-c', '40', '-x',
                           'fasta', 'genomic_filtered_output', 'gene_output'])
    diam = os.path.join(basedir, reference,
                        'gtdb_r89_protein_db.2010-09-27.faa.dmnd')
    tax = os.path.join(basedir, reference, 'gtdb_r89_taxonomy.2019-09-27.tsv')
    subprocess.check_call(['refinem', 'taxon_profile', '-c', '40', '-x',
                           'fna', 'gene_output', sc_stats, diam, tax,
                           'taxon_profile_output_dir'])
    t_filter = os.path.join(basedir, 'taxon_filter.tsv')
    subprocess.check_call(['refinem', 'taxon_filter', '-c', '40',
                           'taxon_profile_outpir_dir', t_filter])
    subprocess.check_call(['refinem', 'filter_bins', '-x', 'fasta',
                           'genomic_filtered_output', t_filter,
                           'taxon_filtered_output_dir'])


if __name__ == "__main__":
    (assembly, bins, bamfiles,
        basedir, reference) = read_refinem_info_file(sys.argv[1])
    # Make the base directory below
    if os.path.exists(basedir):
        pass
    else:
        os.mkdir(basedir)
    os.chdir(basedir)
    run_refinem(assembly, bins, bamfiles, basedir, reference)
