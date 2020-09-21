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
all_processing = []
for one in config['TaxonAddThresh']:
    for two in config['TaxonAddThresh']:
        all_processing.append(f'TaxonRemovedA{one}R{two}')

gff_processing = []
for value in [90, 95, 99]:
    for value2 in [90, 95, 99]:
        for length in [2000, 5000]:
            new_value = f"Full.{length}.TaxonRemovedA{value}R{value2}.ANIRepatT95M20"
            gff_processing.append(new_value)


rule all:
    input:
        all = expand(
            "BinIdentification/{sample}.TaxonRemovedA{add}R{remove}.txt",
            sample=config['samples'],
            add=config['TaxonAddThresh'],
            remove=config['TaxonRemoveThresh']
        ),
        ani = expand(
            "FastANI/Filtered_{length}/Q{query}_R{reference}.{length}_{split}.txt",
            length=config['ANIAssemblyFilterSize'],
            query=config['samples'],
            reference=config['samples'],
            split=config['ANIAssemblySplits']
            ),
        all_ani = expand(
            "FastANI/Filtered_{length}/AllRawOriginalFastANIResults.{length}.txt",
            length=config['ANIAssemblyFilterSize']),
        append_bins = expand(
            "FastANI/Filtered_{length}/AllProcessed.{processing}.FastANIResults.{length}.txt",
            length=config['ANIAssemblyFilterSize'],
            processing=all_processing
        ),
        final_ani = expand(
            "BinIdentification/{sample}.Full.{length}.{processing}.ANIRepatT{thresh}M{match}.txt",
            sample=config['samples'],
            length=config['ANIAssemblyFilterSize'],
            processing=all_processing,
            thresh=config['ANIRepatIdentThreshold'],
            match=config['ANIRepatCountThreshold']
            ),
        gff_results = expand(
            "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.txt",
            sample=config['samples'],
            processing=gff_processing
        )



# General Processing: Create a BinID file from list of .FASTA files
rule create_bin_id_file:
    input:
        bins = "../input/OriginalBins/{sample}/Bin.001.fasta"
    params:
        bins = directory("../input/OriginalBins/{sample}/")
    log:
        "logs/generalProcessing/{sample}.BinIDCreation.log"
    output:
        "BinIdentification/{sample}.OriginalBinID.txt"
    shell:
        """
        python scripts/getContigBinIdentifier.py -f {params.bins}/*.fasta -o {output} -l {log}
        """


# Filter contigs based on Cat/Bat taxonomic scores.
rule filter_taxonomy:
    input:
        bin_id = "BinIdentification/{sample}.OriginalBinID.txt",
        cat = "../input/Cat/{sample}/{sample}.C2C.names.txt",
        bat = "../input/Bat/{sample}/{sample}.Bin2C.names.txt"
    params:
        addThresh = "{add}",
        removeThresh = "{remove}"
    wildcard_constraints:
        add="\d+",
        remove="\d+"
    output:
        new_bin_id = "BinIdentification/{sample}.TaxonRemovedA{add}R{remove}.txt"
    log:
        readme = "logs/taxonFiltering/{sample}.TaxonRemovedA{add}R{remove}.log"
    shell:
        """
        python scripts/taxonFilter.py -i {input.bin_id} -c {input.cat} -b {input.bat} \
        -m {params.removeThresh} -a {params.addThresh} -o {output.new_bin_id} -r {log.readme}
        """


### ~~~~~ fastANI portion of the processing pipeline ~~~~~ ###


# Given an assembly file, filter to only contain contigs > 5kb (or whatever spec. in config file)
rule filter_contigs:
    input:
        assembly = "../input/Assembly/{sample}.Assembly500.fasta"
    params:
        length = "{length}"
    log:
        "logs/ANI/filtering{sample}Assembly{length}.log"
    output:
        outputs = "../input/Assembly/Filtered/{sample}.Assembly{length}.fasta"
    script:
        "scripts/filterContigsSm.py"

# Split up assembly file into many files, each corresponding to 1 fasta entry
rule split_filtered_contigs:
    input:
        assembly = "../input/Assembly/Filtered/{sample}.Assembly{length}.fasta"
    params:
        folder = "../input/Assembly/Filtered/Split-Files-{length}/{sample}"
    output:
        outputs = "../input/Assembly/Filtered/Split-Files-{length}/{sample}/{sample}.{length}.complete.tkn"
    shell:
        """
        scripts/split_mfa.sh {input.assembly} {params.folder}
        touch ../input/Assembly/Filtered/Split-Files-{wildcards.length}/{wildcards.sample}/{wildcards.sample}.{wildcards.length}.complete.tkn
        """

# From the list of contigs from split_filtered_contigs,
# create a list specifying their path.
rule make_contig_list:
    input:
        samples = "../input/Assembly/Filtered/Split-Files-{length}/{sample}/{sample}.{length}.complete.tkn"
    output:
        outputs = "../input/Assembly/Filtered/Split-Files-{length}/{sample}.{length}.AllContigsList.txt"
    script:
        "scripts/makelist.py"


# Split the list files into N smaller files specified by config['ANIAssemblySplitSize']
rule split_lists:
    input:
        lists = "../input/Assembly/Filtered/Split-Files-{length}/{sample}.{length}.AllContigsList.txt"
    params:
        split_size = config['ANIAssemblySplitSize'],
        directory = "../input/Assembly/Filtered/Split-Files-{length}"
    output:
        outputs = "../input/Assembly/Filtered/Split-Files-{length}/{sample}.{length}_{split}"
    script:
        "scripts/splitList.py"


# Run the fastANI program now! using full lists as query and split lists as points to actual files
rule run_fastani:
    input:
        full_lists = "../input/Assembly/Filtered/Split-Files-{length}/{query}.{length}.AllContigsList.txt",
        split_lists = "../input/Assembly/Filtered/Split-Files-{length}/{reference}.{length}_{split}"
    params:
        minfrac = config['FastANIMinFraction'],
        fraglen = config['FastANIFragLength']
    log:
        "logs/FastANI/Q{query}_R{reference}.{length}_{split}.log"
    output:
        outputs = "FastANI/Filtered_{length}/Q{query}_R{reference}.{length}_{split}.txt"
    shell:
        """
        fastANI -t 20 --minFraction {params.minfrac} --fragLen {params.fraglen} --ql {input.full_lists} --rl {input.split_lists} -o {output.outputs} &> {log}
        touch FastANI/Filtered_{wildcards.length}/AniComplete.tkn
        """


rule concatenate_output:
    input:
        files = expand("FastANI/Filtered_{{length}}/Q{query}_R{reference}.{{length}}_{split}.txt",
                       query=config['samples'], reference=config['samples'],
                       split=config['ANIAssemblySplits'])
    output:
        outputs = "FastANI/Filtered_{length}/AllRawOriginalFastANIResults.{length}.txt"
    wildcard_constraints:
        length="\d+"
    shell:
        """
        cat {input.files} > {output.outputs}
        """


### ~~~~~ fastANI REPATRIATION portion of the processing pipeline ~~~~~ ###


# Append bins to the default output from fastANI
rule append_bins_to_ani:
    input:
        ani_file = "FastANI/Filtered_{length}/AllRawOriginalFastANIResults.{length}.txt",
        bin_id = expand(
            "BinIdentification/{sample}.{{processing}}.txt", sample=config['samples'])
    output:
        new_ani = "FastANI/Filtered_{length}/AllProcessed.{processing}.FastANIResults.{length}.txt"
    shell:
        """
        python scripts/appendBinsToANI.py -a {input.ani_file} -b {input.bin_id} -o {output.new_ani}
        """


# Run aniContigRecycler.py on the results and output
rule ani_based_contig_repatriation:
    input:
        ani_file = "FastANI/Filtered_{length}/AllProcessed.{processing}.FastANIResults.{length}.txt"
    params:
        ident_thresh = "{thresh}",
        count_thresh = "{match}",
        bin_directory = "BinIdentification"
#    wildcard_constraints:
#        match = "\d+",
#        thresh = "\d+"
    output:
        bin_files = "FastANI/Filtered_{length}/{processing}.ANIRepatT{thresh}M{match}.{length}.txt",
        out = expand("BinIdentification/{sample}.{{length}}.{{processing}}.ANIRepatT{{thresh}}M{{match}}.txt", sample=config['samples'])
    shell:
        """
        python scripts/aniContigRecycler.py -a {input.ani_file} -t {params.ident_thresh} -m {params.count_thresh} -d {params.bin_directory} -o {output.bin_files}
        """


# This is where we tie in the bin identification to the ANI processing
rule concat_ani_bin_ident:
    input:
        bins_to_add_to = "BinIdentification/{sample}.{processing}.txt",
        ani_bins = "BinIdentification/{sample}.{length}.{processing}.ANIRepat{params}.txt"
    output:
        new_bins = "BinIdentification/{sample}.Full.{length}.{processing}.ANIRepat{params}.txt"
    shell:
        """
        cat {input.bins_to_add_to} {input.ani_bins} > {output.new_bins}
        """


### ~~~~~ blastn REPATRIATION portion of the processing pipeline ~~~~~ ###


# Find the top genomedb_acc feature per contig
rule grab_contig_top_genomedb_acc:
    input:
        gff = "../input/GFF/{sample}/{sample}.All.gff",
        bin_id = "BinIdentification/{sample}.{processing}.txt"
    params:
        attribute = "genomedb_acc",
#        bin_id = "BinIdentification/{sample}.{processing}.Full.txt"
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
    output:
        bin_annotations = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.txt"
    shell:
        """
        python scripts/writeModeGffFeaturePerBin.py {input.contig_annotations} {output.bin_annotations}
        """


# Download all genomedb_acc from rule above.
rule download_genomedb_acc:
    input:
        bin_annotations = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.txt"
    params:
        directory = directory(
            "GFFAnnotation/AssemblyFiles/{sample}_{processing}/")
    output:
        out_tkn = "GFFAnnotation/AssemblyFiles/{sample}_{processing}/NCBI_Assembly_Download.tkn"
    shell:
        """
        sh scripts/download_acc_ncbi.bash {input.bin_annotations} {params.directory}
        """

rule download_genomedb_acc:
    input:
        bin_annotations = "GFFAnnotation/{sample}/{sample}.{processing}.TopBinGenomeDBAcc.txt"
    params:
        directory_name = "Assemblies"
    output:
        outdir=directory("Assemblies")
    shell:
        """
        ./datasets download assembly GCA_003269275.1,GCF_000243215.1,GCF_001314995.1 --filename assemblydownloads.zip
        unzip assemblydownloads.zip -d {params.directory_name}
        """
