'''
Snakemake file that starts with the RAW output from fastANI
and creates a list of predicted repatriated contigs along
with a directory containing fasta entries for each bin that
had contigs to be repatriated. These fasta files ONLY contain
the new contigs and must be merged with the original fasta files
downstream.

Input: raw ani file, original assembly for each sample analyzed,
and a bin identifier file (binID\tcontigName\n)
'''

configfile: "config/contigrecyc.yaml"

# Input requires:
# 1. fastANI output
# 2. A file in the format: bin\tcontig\n
# 3. Pass

# Script files located under ./contigRepatration/scripts/

rule all:
    input:
        expand("RepatFastas_{ani_filename}_T{t}_M{m}", ani_filename=config["anifile"],
                                                       t=config["threshold"],
                                                       m=config["matches"])

# Add bin information to raw ANI output from fastANI snakemake.
rule append_bins:
    input:
        ani="{ani_filename}.txt"
    output:
        outfile="Master.{ani_filename}.txt"
    shell:
        "python scripts/appendBinsToANI.py -a {input.ani} -b *.bins.txt -o {output.outfile}"

# Analyze contigs for repatriation and output a file in the format:
# f"query\treference\tcontig\tbin_to_match\tMATCH\n"
rule write_repatrated_contigs:
    input:
        ani="Master.{ani_filename}.txt"
    params:
        threshold="{t}",
        matches="{m}"
    output:
        outfile="Repat.{ani_filename}.T{t}.M{m}.txt"
    shell:
        "python scripts/aniContigRecycler.py -a {input.ani} -t {params.threshold} -m {params.matches} -o {output.outfile}"

# No need to write FASTA files and above script writes a new bin ID file
'''
rule write_new_fasta_files:
    input:
        repat="Repat.{ani_filename}.T{t}.M{m}.txt",
        assembly=config["assembly"]
    params:
        samplename=config['samplename']
    output:
        repatdir=directory("RepatFastas_{ani_filename}_T{t}_M{m}")
    shell:
        """
        mkdir {output.repatdir}
        python scripts/recycledBinFastaWriter.py -r {input.repat} -s {input.assembly} -q {params.samplename} -d {output.repatdir}
        """
'''
