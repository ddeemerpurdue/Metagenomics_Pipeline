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

supernatant_bins = ["{0:03}".format(i) for i in range(1, 255)]
particle_bins = ["{0:03}".format(i) for i in range(1, 235)]

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
        ),
        # For now, needing 2 separate expand functions (1 per sample)
        blast_bin_p = expand(
            "blast_bin_{sample}/blastn.{number}.txt",
            sample='particle', number=particle_bins),
        blast_bin_s = expand(
            "blast_bin_{sample}/blastn.{number}.txt",
            sample='supernatant', number=supernatant_bins),
        blast_nobin_p = expand(
            "blast_nobin_{sample}/blastnNoBins.{number}.txt",
            sample='particle', number=particle_bins),
        blast_nobin_s = expand(
            "blast_nobin_{sample}/blastnNoBins.{number}.txt",
            sample='supernatant', number=supernatant_bins)


# General Processing: Create a BinID file from list of .FASTA files
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
        python scripts/getContigBinIdentifier.py -f {params.bins}*.fasta -o {output} -l {log}
        """


# Filter contigs based on Cat/Bat taxonomic scores.
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


### ~~~~~ fastANI portion of the processing pipeline ~~~~~ ###


# Given an assembly file, filter to only contain contigs > 5kb (or whatever spec. in config file)
rule filter_contigs:
    input:
        assembly = "../input/Assembly/{sample}.original500.fasta"
    params:
        length = "{length}"
    log:
        "logs/generalProcessing/filtering{sample}Assembly{length}.log"
    wildcard_constraints:
        length = "\d+"
    output:
        outputs = "../input/Assembly/Filtered/{sample}.Assembly{length}.fasta"
    shell:
        "python scripts/filterSeqLength.py -a {input.assembly} -l {params.length} -o {output.outputs} -g {log}"


# Split up assembly file into many files, each corresponding to 1 fasta entry
# and write all output to a list, along with N subsetted lists.
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

# Thought - test out global wildcards
# split_fastas=glob_wildcards("../input/Assembly/Filtered/Split-Files-{length}/{sample}/{file}.fasta")
# Then make this temporary input in run_fastani
# Make a wildcard for all files in the above 'files' output:
split_fasta_files, = glob_wildcards(
    "../input/Assembly/Filtered/Split-Files-2000/{name}.fasta")
# Run the fastANI program now! using full lists as query and split lists as points to actual files
rule run_fastani:
    input:
        queryfiles = temp(split_fasta_files),
        reffiles = temp(split_fasta_files),
        full_lists = "../input/Assembly/Filtered/Split-Files-{length}/{query}.AllContigsList{length}.txt",
        split_lists = "../input/Assembly/Filtered/Split-Files-{length}/{reference}.AllContigsList{length}_{split}.txt"
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
        length = "\d+"
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
        out = expand(
            "BinIdentification/{sample}.{{length}}.{{processing}}.ANIRepatT{{thresh}}M{{match}}.txt", sample=config['samples'])
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


'''
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
'''


rule write_Fasta_from_binID:
    input:
        binID = "BinIdentification/{sample}.{processing}.txt"
    output:

    shell:
        """
        """

# Blast all binned contigs against their reference assembly
rule blast_binners:
    input:
        bin_file = "Fastas/{sample}TRA90R99ANIT95M20/Bin.{number}.fasta",
        assembly = "ReferenceFiles/{sample}/{number}.fasta"
    params:
        fmt = "6"
    output:
        outfile = "blast_bin_{sample}/blastn.{number}.txt"
    shell:
        """
        blastn -query {input.bin_file} -subject {input.assembly} -outfmt {params.fmt} -out {output.outfile}
        """

# Blast non-binning contigs against all reference assemblies
rule blast_nonbinners:
    input:
        nonbinner = "Fastas/{sample}TRA90R99ANIT95M20/Bin.NoBin.fasta",
        assembly = "ReferenceFiles/{sample}/{number}.fasta"
    params:
        fmt = "6"
    output:
        outfile = "blast_nobin_{sample}/blastnNoBins.{number}.txt"
    shell:
        """
        blastn -query {input.nonbinner} -subject {input.assembly} -outfmt {params.fmt} -out {output.outfile}
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
