configfile: "config/config.yaml"
'''
configfile should have the following attributes:
Step 0:
- assembly
- index  # Unique identifier per sample
- threads
- readpath
Step 1:
- readsize
- threshold
- maxbin_abunds
- refinem_bams
'''


rule all:
    input:
        expand("AssemblyIndex/{index}.1.bt2", index=config['index'])
        expand("Bam/{sample}.bam", sample=config['samples'])


# PART 0: Essentials
# [Assembly Indexing - Alignment - Bam Sorting - Bam Indexing]
# Speed: **
''' Essentials Start '''

# index_assembly: index the assembly (or any fasta file)
rule index_assembly:
    input:
        config["assembly"]
    params:
        index=config["index"],
        threads=config["threads"]
    conda:
        "envs/alignment.yaml"
    output:
        f"AssemblyIndex/{config['index']}.1.bt2"
    log:
        f"logs/assemblyIndexing.{config['index']}.log"
    shell:
        """
        bowtie2-build -f --threads {params.threads} {input} AssemblyIndex/{params.index} >> {log}
        """

# bowtie2_alignment: align fastq files to indexed assembly from above rule
rule bowtie2_alignment:
    input:
        reads=expand("{path}{{sample}}_{num}_10000.fastq", path=config['readpath'], num=['R1','R2']),
        i=f"AssemblyIndex/{config['index']}.1.bt2"
    params:
        index=f"AssemblyIndex/{config['index']}",
        threads=config['threads']
    conda:
        "envs/alignment.yaml"
    output:
        temp("Bam/{sample}.bam")
    log:
        "logs/bowtie2Alignment.{sample}.log"
    shell:
        "(bowtie2 --threads {params.threads} -x {params.index} -1 {input.reads[0]} -2 {input.reads[1]} | samtools view -b -o {output}) &> {log}"

# bam_sorting: sort bam file created from above rule
rule bam_sorting:
    input:
        "Bam/{sample}.bam"
    params:
        prefix="Bam/pref.{sample}",
        threads=config['threads']
    conda:
        "envs/alignment.yaml"
    output:
        protected("Bam/{sample}.sorted.bam")
    log:
        "logs/samtoolsSort.{sample}.log"
    shell:
        "samtools sort -@{params.threads} -T {params.prefix} -o {output} {input} &> {log}"

# bam_indexing: Index the bam file (required for downstream programs)
rule bam_indexing:
    input:
        "Bam/{sample}.sorted.bam"
    conda:
        "envs/alignment.yaml"
    output:
        "Bam/{sample}.sorted.bam.bai"
    log:
        "logs/samtoolsIndex.{sample}.log"
    shell:
        "samtools index -@20 {input} &> {log}"

''' Essentials Complete '''
''' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '''

# PART 1: Base Processing
# [Bam Filtering - Calculating Coverage - Abundance Filter - Abundance File Merge]
# Speed: ****
''' Base Processing Start '''

# filter_bams: filter BAM files to only contain reads mapping above a certain threshold of percent identity.
rule filter_bams:
    input:
        "Bam/{sample}.sorted.bam"
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

# bbmap_stats: Calculate coverage stats per contig
rule bbmap_stats:
    input:  # Wildcard path used for filtered or normal: for normal, path == "."
        "Bam/{path}/{sample}.sorted.bam"
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
        samples=lambda wildcards: expand(f"Abundances/{wildcards.path}/{{sample}}.abundance.txt", sample=config['maxbin_abunds2'])
    output:
        "Abundances/{path}/abundance_list.txt"
    script:
        "scripts/make_abund_list.py"

rule run_maxbin:
    input:
        abund="Abundances/{path}/abundance_list.txt",
        scaff=config["assembly"]
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

rule run_checkm:
    input:
        directory("Maxbin/{path}/MBALL")
    output:
        directory("{input}/CHECKM")
    shell:
        "checkm lineage_wf -t 8 -x fasta {input} {output}"

rule run_refinem:
    input:
#        checkm="CheckM/All",
        abunds=expand(f"Abundances/Default/{{files}}.abundance.txt", files=config['maxbin_abunds']),
        mb=directory("Maxbin/{path}"),
    params:
        assembly=config["assembly"],
        bins=directory("Maxbin/{path}"),
        bams=lambda wildcards: expand(f"Bam/{wildcards.path}/{{ref_bams}}.sorted.bam", ref_bams=config['refinem_bams']),
        rm_basedir=config["rm_basedir"],
        reference=config["reference"]
    conda:
        "envs/refinem.yaml"
    output:
        directory("RefineM/{path}")
    shell:
        "python scripts/refinem_snakemake.py -a {params.assembly} -b {params.bins} -m {params.bams} -s {params.rm_basedir} -r {params.reference}"
