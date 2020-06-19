import os
import subprocess
import argparse


parser = argparse.ArgumentParser(description="Parser")
parser.add_argument("-a", "--Assembly", help="Assembly file",
                    required=True)
parser.add_argument("-b", "--Bins", help="Bin location (directory)",
                    required=True)
parser.add_argument("-m", "--Bams", help="Bam files location - 1+",
                    required=True, nargs='*')
parser.add_argument("-s", "--Basedir", help="Base directory to save \
                    refinement files to", required=True)
parser.add_argument("-r", "--Reference", help="Path to the directory \
                    containing the database information required for refinem",
                    required=True)

argument = parser.parse_args()


# Runthrough RefineM
def run_refinem(assembly, bins, bams, basedir, reference):
    # Remove based on genomic properties
    stats_output_dir = os.path.join(basedir, 'stats_output_dir')
    rf = []
    rf.extend(['refinem', 'scaffold_stats', '-c', '20',
               '-x', 'fasta', assembly, bins, stats_output_dir])
    rf = rf + bams
    subprocess.check_call(rf)
    sc_stats = os.path.join(stats_output_dir, 'scaffold_stats.tsv')
    outlier_output = os.path.join(basedir, 'outlier_output')
    subprocess.check_call(['refinem', 'outliers', sc_stats,
                           outlier_output, '--no_plots'])
    outliers = os.path.join(basedir, 'outlier_output', 'outliers.tsv')
    genomic_filtered_output = os.path.join(basedir, 'genomic_filtered_output')
    subprocess.check_call(['refinem', 'filter_bins', '-x', 'fasta',
                           bins, outliers, genomic_filtered_output])
    # Remove based on taxnomonic properties
    gene_output = os.path.join(basedir, 'gene_output')
    subprocess.check_call(['refinem', 'call_genes', '-c', '40', '-x',
                           'fasta', genomic_filtered_output, gene_output])
    diam = os.path.join(reference,
                        'gtdb_r89_protein_db.2019-09-27.faa.dmnd')
    tax = os.path.join(reference, 'gtdb_r89_taxonomy.2019-09-27.tsv')
    taxon_profile_output_dir = os.path.join(basedir,
                                            'taxon_profile_output_dir')
    subprocess.check_call(['refinem', 'taxon_profile', '-c', '40', '-x',
                           'fna', gene_output, sc_stats, diam, tax,
                           taxon_profile_output_dir])
    t_filter = os.path.join(basedir, 'taxon_filter.tsv')
    subprocess.check_call(['refinem', 'taxon_filter', '-c', '40',
                           taxon_profile_output_dir, t_filter])
    taxon_filtered = os.path.join(basedir, 'taxon_filtered_output_dir')
    subprocess.check_call(['refinem', 'filter_bins', '-x', 'fasta',
                           genomic_filtered_output, t_filter,
                           taxon_filtered])


if __name__ == "__main__":
    if os.path.exists(argument.Basedir):
        pass
    else:
        os.mkdir(argument.Basedir)
    os.chdir(argument.Basedir)
    run_refinem(argument.Assembly, argument.Bins,
                argument.Bams, argument.Basedir,
                argument.Reference)
