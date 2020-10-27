'''
Author: Dane
Snakemake rules for ANI bin processing/repatriation.
'''

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