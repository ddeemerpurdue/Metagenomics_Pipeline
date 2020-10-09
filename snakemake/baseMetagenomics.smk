'''
Authors: Dane Deemer and Renee Oles
Purpose: Snakemake pipeline that takes fastq DNA reads and bins
using a reference assebmly 
Input: 
- Samples 
- Index (sample identifiers)
- Threads 
- Threshold
- Readsize
results:
- readsize
- threshold
- maxbin_abunds
- refinem_bams
'''
configfile: "../config/config.yaml"


rule all:
    input:
        expand("AssemblyIndex/{index}.1.bt2", index=config['index']),
        expand("Bam/Original/{sample}.sorted.bam.bai", sample=config['samples'])


# ~~~~~~~~~~ STEP 0: Essentials ~~~~~~~~~~ #
# [Assembly Indexing - Alignment - Bam Sorting - Bam Indexing]

# Filter fastq files
'''
rule filter_fastq:
    input:
        reads=expand("FastqFiltered/{sample}_{num}.fastq", sample=config['samples'], num=["R1","R2"])
    params:
        threads=config["threads"]
    conda:
        "envs/fastqc.yaml"
    output:
        directory("FastqFiltered")
    log:
        "logs/filter_fastq.{sample}_{num}.log"
    shell:
        "fastqc -t {params.threads}  --outdir {output} {input.reads} >> {log}"
'''

# Index the assembly (or any fasta file)
rule index_assembly:
    input:
        #needed=directory("FastqFiltered"),
        assembly="../input/Assembly/{index}.fasta"
    params:
        threads=config["threads"],
        index="AssemblyIndex/{index}"
    conda:
        "envs/alignment.yaml"
    output:
        "AssemblyIndex/{index}.1.bt2"
    log:
        "logs/assemblyIndexing.{index}.log"
    shell:
        "bowtie2-build -f --threads {params.threads} {input.assembly} {params.index} >> {log}"

# Align fastq files to indexed assembly from above rule
rule bowtie2_alignment:
    input:
        reads=expand("../input/Fastq/{{sample}}_{num}.fastq", num=["R1","R2"]),
        i=f"AssemblyIndex/{config['index']}.1.bt2"
    params:
        index=f"AssemblyIndex/{config['index']}",
        threads=config['threads']
    conda:
        "envs/alignment.yaml"
    output:
    	"Bam/Original/{sample}.bam"
    log:
        "logs/bowtie2Alignment.{sample}.log"
    shell:
        "(bowtie2 --threads {params.threads} -x {params.index} -1 {input.reads[0]} -2 {input.reads[1]} | samtools view -b -o {output}) &> {log}"

# Sort bam file created from above rule
rule bam_sorting:
    input:
        "Bam/Original/{sample}.bam"
    params:
        prefix="Bam/Original/pref.{sample}",
        threads=config['threads']
    conda:
        "envs/alignment.yaml"
    output:
        "Bam/Original/{sample}.sorted.bam"
    log:
        "logs/samtoolsSort.{sample}.log"
    shell:
        "samtools sort -@{params.threads} -T {params.prefix} -o {output} {input} &> {log}"

# bam_indexing: Index the bam file (required for downstream programs)
rule bam_indexing:
    input:
        "Bam/Original/{sample}.sorted.bam"
    conda:
        "envs/alignment.yaml"
    output:
        "Bam/Original/{sample}.sorted.bam.bai"
    log:
        "logs/samtoolsIndex.{sample}.log"
    shell:
        "samtools index -@20 {input} &> {log}"


# ~~~~~~~~~~ STEP 1: General Processing ~~~~~~~~~~ #
# [Bam Filtering - Calculating Coverage - Abundance Filter - Abundance File Merge]

# Filter BAM files to only contain reads mapping above a certain threshold of percent identity.
rule filter_bams:
    input:
        "Bam/Original/{sample}.sorted.bam"
    params:
        readsize=config['readsize'],
        threshold=config['threshold']
    conda:
        "envs/alignment.yaml"
    output:
        "Bam/Filtered/{sample}.{params.threshold}F.sorted.bam"
    log:
        "logs/bamFiltering.{sample}.log"
    shell:
        "(samtools view -h {input} | python scripts/sam_threshold_filter.py -s {params.readsize} -t {params.threshold} | samtools view -b -o {output} && samtools index -@20 {output}) &> {log}"

# Calculate coverage stats per contig
rule bbmap_stats:
    input:
        samples=lambda wildcards: expand(f"Bam/{wildcards.path}/{{sample}}.sorted.bam")
    conda:
        "envs/bbmap.yaml"
    output:
        "Abundances/{path}/{sample}.coverage.txt"
    log:
        "log/bbmap-abunds.{path}.{sample}.log"
    shell:
        "pileup.sh in={input} out={output} &> {log}"

rule get_cov_file:
    input:
        "Abundances/{path}/{sample}.coverage.txt"
    output:
        "Abundances/{path}/{sample}.abundance.txt"
    shell:
        "cut -f 1,5 {input} | grep -v '^#' > {output} 2> {log}"

rule get_abund_list:
    input:
        samples=lambda wildcards: expand(f"Abundances/{wildcards.path}/{{sample}}.abundance.txt")
    output:
        "Abundances/{path}/abundance_list.txt"
    script:
        "scripts/make_abund_list.py"

rule run_maxbin:
    input:
        abund="Abundances/{path}/abundance_list.txt",
        scaff="../input/Assembly/{index}.fasta"
    conda:
        "envs/binning.yaml"
    output:
        directory("Maxbin/{path}/MBALL")
    log:
        log1="logs/maxbin2_{path}_.log"
    shell:
        """
        module load MaxBin/2.2.3
        run_MaxBin.pl -contig {input.scaff} -out {output} -abund_list {input.abund} &> {log}
        touch {output}
        """

# Estimate completion and contamination for original binning
rule run_checkm:
    input:
        "Maxbin/{path}/{path}.log"
    conda:
        "envs/binning.yaml"
    params:
        threads=config['threads'],
        out="CheckM/{path}/"
    output:
        "CheckM/{path}/lineage.ms"
    log:
        "logs/checkm_{path}.log"
    shell:
        """
        checkm lineage_wf -t {params.threads} -x fasta {input} {params.out} &> {log}
        checkm qa -t {params.threads} {output} {params.out} &> {log}
        """

# Further refine bins using refineM
rule run_refinem:
    input:
        bam=expand("Bam/Original/{sample}.sorted.bam", sample=config['samples']),
    params:
        assembly="../input/Aseembly/{index}.fasta",
        reference=config["reference"],
        mb=directory("Maxbin/")
    conda:
        "envs/refinem.yaml"
    output:
        directory("RefineM")
    log:
        "logs/refinem.log"
    shell:
        "python scripts/refinem_snakemake.py -a {params.assembly} -b ../{input.mb} -m ../{input.bam} -s RefineM -r {params.reference} &> {log}"

'''
# Change the below variable!!!
#rm_path = "Default"
rule run_refinem:
    input:
        bam=expand("Bam/{mypath}/{sample}.sorted.bam", mypath=rm_path, sample=config['samples']),
        mb=directory("Maxbin/{path}/")
    params:
        assembly="../input/Aseembly/{index}.fasta",
        reference=config["reference"]
    conda:
        "envs/refinem.yaml"
    output:
        directory("RefineM/{path}")
    log:
        "logs/refinem_{path}.log"
    shell:
        "python scripts/refinem_snakemake.py -a {params.assembly} -b ../{input.mb} -m ../{input.bam} -s RefineM -r {params.reference} &> {log}"

 Complete this last part later 
# Estimate completion and contamination after refinement
rule run_checkm:
    input:
        "RefineM/{path}/{path}.log"
    conda:
        "envs/binning.yaml"
    params:
        threads=config['threads'],
        out="CheckM/{path}/"
    output:
        "CheckM/{path}/lineage.ms"
    log:
        "logs/checkm_{path}.log"
    shell:
        """
        checkm lineage_wf -t {params.threads} -x fasta {input} {params.out} &> {log}
        checkm qa -t {params.threads} {output} {params.out} &> {log}
        """
Complete this last part later '''
