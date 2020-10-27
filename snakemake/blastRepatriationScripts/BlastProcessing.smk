'''
Author: Dane
Snakemake rules for BlastN bin processing/repatriation.
'''
p_nums, = glob_wildcards(
    "BlastBinners/particle/Full10000_OriginalTaxonRemovedA80R99ANIRepatT95M20/blastn.NoBin_{num}.txt")


s_nums, = glob_wildcards(
    "BlastBinners/supernatant/Full10000_OriginalTaxonRemovedA80R99ANIRepatT95M20/blastn.NoBin_{num}.txt")


# ~~~~~~~~~~ STEP 3: BlastN Processing ~~~~~~~~~~ #


# ~~~ Part A: Prepping for Repatriation ~~~ #
# Find the top genomedb_acc feature per contig
rule grab_contig_top_genomedb_acc:
    input:
        gff = "../input/GFF/{sample}/{sample}.All.gff",
        bin_id = "BinIdentification/{sample}.{processing}.txt"
    params:
        attribute = "genomedb_acc",
        bin_id = "BinIdentification/{sample}.{processing}.Full.txt"
    output:
        out_file = "GFFAnnotation/{sample}/{sample}.{processing}.TopContigGenomeDBAcc.txt"
    shell:
        """
        python scripts/gffMine.py -g {input.gff} -a {params.attribute} -b {input.bin_id} -o {output.out_file} --Top
        """


# Find the top genomedb_acc feature per bin
rule grab_bin_top_genomedb_acc:
    input:
        contig_annotations = "GFFAnnotation/{sample}/{sample}.{processing}.TopContigGenomeDBAcc.txt"
    log:
        "logs/BlastRepatration/{sample}.{processing}.TopBinGenomeDBAcc.log"
    output:
        bin_annotations = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.txt"
    shell:
        """
        python scripts/writeModeGffFeaturePerBin.py -a {input.contig_annotations} -o {output.bin_annotations} -l {log}
        """


# Download all genomedb_acc from rule above.
# Assemblies download in the form: {bin_num}.{assembly_#}.fasta
rule download_genomedb_acc:
    input:
        bin_annotations = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.txt"
    params:
        assembly_database = protected(directory(config['AssembliesDatabase']))
    log:
        "logs/BlastRepatriation/assemblyDownloading-{sample}-{processing}.log"
    output:
        "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.success.txt"
    shell:
        """
        python scripts/downloadAssemblies.py -g {input.bin_annotations}  -a {params.assembly_database}/{wildcards.sample}/{wildcards.processing} -l {log}
        """


# Could add a global wildcard to make sure same amount of bins are produced
rule write_Fasta_from_binID:
    input:
        binID = "BinIdentification/{sample}.{processing}.txt",
        assembly = "../input/Assembly/{sample}.original500.fasta"
    params:
        fasta_directory = directory("Fastas/{sample}/{processing}/"),
    output:
        fasta_bin = directory("Fastas/{sample}/{processing}/")
    shell:
        """
        python scripts/writeFastaFromBinID.py -b {input.binID} -a {input.assembly} -f {params.fasta_directory}
        """


# ~~~ Part B: Blasting and Repatriating ~~~ #
rule blast_binners:
    input:
        fasta_file = "Fastas/{sample}/{processing}/Bin.{number}.fasta",
        reference = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.success.txt"
    params:
        fmt = "6"
    log:
        "logs/BlastRepatriation/blastn_{sample}.{processing}.{number}.log"
    output:
        outfile = "BlastBinners/{sample}/{processing}/blastn.{number}.txt"
    shell:
        """
        python scripts/blastAssemblies.py -query {input.fasta_file} -f {params.fmt} -l {log}
        """


rule blast_no_binners:
    input:
        fasta_file = "Fastas/{sample}/{processing}/Bin.NoBin.fasta",
        reference = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.success.txt"
    params:
        fmt = "6"
    log:
        "logs/BlastRepatriation/blastn_{sample}.{processing}.NoBin.log"
    output:
        outfile = "BlastBinners/{sample}/{processing}/blastn.NoBins.success.txt"
    shell:
        """
        python scripts/blastAssemblies.py --NoBin -q {input.fasta_file} -f {params.fmt} -l {log}
        """


# Repatriate contigs based on reference assemblies
# Will have to change the NoBin_{num} to {num}.NoBin
rule blast_repatriation:
    input:
        nonbin_results = expand(
            "BlastBinners/{{sample}}/{{processing}}/blastn.NoBin_{num}.txt",
            num=[]),
        bin_results = expand(
            "BlastBinners/{{sample}}/{{processing}}/blastn.{num}.txt",
            num=[])
    params:
        threshold = "85"
    log:
        "logs/BlastRepatriation/{sample}.{processing}_BlastnT85L2000.log"
    output:
        binid = "BinIdentification/{sample}.{processing}_BlastnT85L2000.txt"
    shell:
        """
        python blastContigRecycler.py -n {input.nonbin_results} -b {input.bin_results} -t {params.threshold}
        """
