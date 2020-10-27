'''
Snakemake rules to run .fasta files against a CAZyme database
using HMMER
'''

rule create_bin_id_file:
    input:
        bins = "FinalBinSets/{sample}/Bin.{number}.fasta"
    params:
        dm = "Results/{sample}.Bin.{number}.dm"
        db = "dbCAN-fam-HMMs.txt"
    output:
        "Results/{sample}.Bin.{number}.out"
    shell:
        """
        hmmscan --domtblout {params.dm} {params.db} {input.bins} > {output}
        """
