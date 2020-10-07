'''
Author: Dane
Purpose: Snakemake pipeline that processes multiple samples bin files
automatically. Various parameters can be specified in the config/config.yaml
file and metadata will automatically be logged to compare various parameters
effects on bin processing.
Starting input requires:
1. Bin files in 1+ directories (.fasta)
2. MetaErg annotation - {sample}.gff file
3. CatBat annotation files {sample}.C2C.names.txt & {sample}.Bin2C.names.txt
'''

configfile: "../config/config.yaml"

# A list of all possible taxon processing on original results.
all_taxon_processing = []
for one in config['TaxonAddThresh']:
    for two in config['TaxonAddThresh']:
        all_processing.append(f'OriginalTaxonRemovedA{one}R{two}')

gff_processing = []
for value in [90, 95, 99]:
    for value2 in [90, 95, 99]:
        for length in [2000, 5000]:
            new_value = f"Full.{length}.TaxonRemovedA{value}R{value2}.ANIRepatT95M20"
            gff_processing.append(new_value)

supernatant_bins = ["{0:03}".format(i) for i in range(1, 255)]
particle_bins = ["{0:03}".format(i) for i in range(1, 235)]


# Start of all of the rules #
rule all:
    input:
        taxonfilter = expand(
            "BinIdentification/{sample}.{processing}TaxonRemovedA{add}R{remove}.txt",
            sample=config['samples'],
            add=config['TaxonAddThresh'],
            remove=config['TaxonRemoveThresh'],
            processing=config['TaxonProcessing']
        ),
        filter = expand(
            "../input/Assembly/Filtered/{sample}.Assembly{length}.fasta",
            length=config['ANIAssemblyFilterSize'],
            sample=config['samples']
            ),
        fastani = expand(
            "FastANI/Filtered_{length}/Q{query}_R{reference}.{length}_{split}.txt",
            length=10000,
            query=config['samples'],
            reference=config['samples'],
            split=config['ANIAssemblySplits']),
        last_fastani = expand("BinIdentification/{sample}.Full{length}_{processing}ANIRepatT{thresh}M{match}.txt",
            sample=config['samples'],
            length=config['ANIAssemblyFilterSize'],
            processing=all_taxon_processing,
            thresh=config['ANIRepatIdentThreshold'],
            match=config['ANIRepatCountThreshold'])


# ~~~~~~~~~~ STEP 0: General Processing ~~~~~~~~~~ #


# Create a BinID file from list of .FASTA files
rule create_bin_id_file:
    input:
        bins = "../input/OriginalBins/{sample}/Bin.001.fasta"
    params:
        bins = directory("../input/OriginalBins/{sample}/")
    log:
        "logs/generalProcessing/{sample}.BinIDCreation.log"
    output:
        protected("BinIdentification/{sample}.Original.txt")
    shell:
        """
        python scripts/getContigBinIdentifier.py -f {params.bins}/*.fasta -o {output} -l {log}
        """


# ~~~~~~~~~~ STEP 1: Taxonomic Processing ~~~~~~~~~~ #


# Add and remove contigs based on taxonomies of contigs and respective bins
rule filter_taxonomy:
    input:
        bin_id = "BinIdentification/{sample}.{processing}.txt",
        cat = "../input/Cat/{sample}/{sample}.C2C.names.txt",
        bat = "../input/Bat/{sample}/{sample}.Bin2C.names.txt"
    params:
        addThresh = "{add}",
        removeThresh = "{remove}"
    wildcard_constraints:
        add = "\d+",
        remove = "\d+"
    output:
        new_bin_id = "BinIdentification/{sample}.{processing}TaxonRemovedA{add}R{remove}.txt"
    log:
        readme = "logs/taxonFiltering/{sample}.{processing}TaxonRemovedA{add}R{remove}.log"
    shell:
        """
        python scripts/taxonFilter.py -i {input.bin_id} -c {input.cat} -b {input.bat} \
        -m {params.removeThresh} -a {params.addThresh} -o {output.new_bin_id} -r {log.readme}
        """


# ~~~~~~~~~~ STEP 1: ANI-Based Processing ~~~~~~~~~~ #


# Filter the assembly to contain only contigs >= n
# This is prep for fastANI
rule filter_contigs:
    input:
        assembly = "../input/Assembly/{sample}.original500.fasta"
    params:
        length = "{length}"
    log:
        "logs/FastANI/filtering{sample}Assembly{length}.log"
    wildcard_constraints:
        length = "\d+"
    output:
        outputs = "../input/Assembly/Filtered/{sample}.Assembly{length}.fasta"
    shell:
        "python scripts/filterSeqLength.py -a {input.assembly} -l {params.length} -o {output.outputs} -g {log}"


# Prep for fastANI by splitting the assembly into multiple subsetted files
# Output consists of a .fasta file for EVERY entry in the assembly, a list
# of all the locations to said files, and n splits of that list file.
rule split_filtered_contigs:
    input:
        assembly = "../input/Assembly/Filtered/{sample}.Assembly{length}.fasta"
    params:
        parts = config['ANIAssemblySplitSize']
    log:
        "logs/generalProcessing/{sample}.{length}.FastaEntrySplitting.log"
    output:
        files = directory(
            "../input/Assembly/Filtered/Split-Files-{length}/{sample}/"),
        filelist = "../input/Assembly/Filtered/Split-Files-{length}/{sample}.AllContigsList{length}.txt",
        splitlist = expand("../input/Assembly/Filtered/Split-Files-{{length}}/{{sample}}.AllContigsList{{length}}_{split}.txt", split=range(
            1, int(config['ANIAssemblySplitSize']) + 1))
    shell:
        """
        python scripts/splitFastaByEntry.py -a {input.assembly} -o {output.files} -l {output.filelist} -n {params.parts} -g {log}
        """


# Run the fastANI program and get results for all splits.
rule run_fastani:
    input:
        full_list = "../input/Assembly/Filtered/Split-Files-{length}/{query}.AllContigsList{length}.txt",
        split_list = "../input/Assembly/Filtered/Split-Files-{length}/{reference}.AllContigsList{length}_{split}.txt"
    params:
        minfrac = config['FastANIMinFraction'],
        fraglen = config['FastANIFragLength']
    log:
        "logs/FastANI/Q{query}_R{reference}.{length}_{split}.log"
    output:
        outputs = "FastANI/Filtered_{length}/Q{query}_R{reference}.{length}_{split}.txt"
    shell:
        """
        fastANI -t 20 --minFraction {params.minfrac} --fragLen {params.fraglen} --ql {input.full_list} --rl {input.split_list} -o {output.outputs} &> {log}
        """


# Concatenate all of the parallelized fastANI results and remove all the split
# fasta entries (see temporary output file)
rule concatenate_output:
    input:
        files = expand("FastANI/Filtered_{{length}}/Q{query}_R{reference}.{{length}}_{split}.txt",
                       query=config['samples'], reference=config['samples'],
                       split=config['ANIAssemblySplits'])
    output:
        split_files = temp(directory(expand(
            "../input/Assembly/Filtered/Split-Files-{{length}}/{sample}/", sample=config['samples']))),
        outputs = "FastANI/Filtered_{length}/AllRawOriginalFastANIResults.{length}.txt"
    wildcard_constraints:
        length = "\d+"
    shell:
        """
        cat {input.files} > {output.outputs}
        """


# Append bins to the default output from fastANI
rule append_bins_to_ani:
    input:
        ani_file = "FastANI/Filtered_{length}/AllRawOriginalFastANIResults.{length}.txt",
        bin_id = expand(
            "BinIdentification/{sample}.{{processing}}.txt", sample=config['samples'])
    log:
        "logs/FastANI/BinAppending_{length}.{processing}.FastANI.log"
    output:
        new_ani = "FastANI/Filtered_{length}/AllProcessed.{processing}.FastANIResults.{length}.txt"
    shell:
        """
        python scripts/appendBinsToANI.py -a {input.ani_file} -b {input.bin_id} -o {output.new_ani} -l {log}
        """


# The the ANI-based repatriation script to add new contigs
rule ani_based_contig_repatriation:
    input:
        ani_file = "FastANI/Filtered_{length}/AllProcessed.{processing}.FastANIResults.{length}.txt"
    params:
        ident_thresh = "{thresh}",
        count_thresh = "{match}",
        bin_directory = "BinIdentification"
    wildcard_constraints:
        match = "\d+",
        thresh = "\d+"
    log:
        "logs/FastANI/ContigRepatriation_{length}_{processing}.ANIRepatT{thresh}M{match}.{length}.log"
    output:
        bin_files = "FastANI/Filtered_{length}/{processing}.ANIRepatT{thresh}M{match}.{length}.txt",
        out = expand(
            "BinIdentification/{sample}.{{length}}.{{processing}}.ANIRepatT{{thresh}}M{{match}}.txt", sample=config['samples'])
    shell:
        """
        python scripts/aniContigRecycler.py -a {input.ani_file} -t {params.ident_thresh} -m {params.count_thresh} -d {params.bin_directory} -o {output.bin_files}
        """


# This is where ANI gets merged with previous bin identification files, since
# the ANI-repat step only write a binID file with the new contigs.
rule concat_ani_bin_ident:
    input:
        bins_to_add_to = "BinIdentification/{sample}.{processing}.txt",
        ani_bins = "BinIdentification/{sample}.{length}.{processing}.ANIRepat{params}.txt"
    output:
        new_bins = "BinIdentification/{sample}.Full{length}_{processing}ANIRepat{params}.txt"
    shell:
        """
        cat {input.bins_to_add_to} {input.ani_bins} > {output.new_bins}
        """


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
        python scripts/blastAssemblies.py -NoBin -query {input.fasta_file} -f {params.fmt} -l {log}
        """

# Repatriate contigs based on reference assemblies
rule blast_repatriation:
    input:
        nonbin_results = expand(
            "blast_nobin_{sample}/blastnNoBins.{num}.txt", sample='supernatant', num=supernatant_bins),
        bin_results = expand("blast_bin_{sample}/blastn.{num}.txt",
        sample='supernatant', num=supernatant_bins)
    params:
        threshold = "85"
    output:
        log = "logs/BlastNResults.T85.txt",
        binid = "BinIdentification/BlastNResults.T85.txt"
    shell:
        """
        python blastContigRecycler.py -n {input.nonbin_results} {input.bin_results} -t {params.threshold}
<<<<<<< HEAD
        """
== == == =
        """
>>>>>>> origin/edits_dgd
