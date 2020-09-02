'''
Author: Dane Deemer
Purpose: To take in 1 or more .Bam files and get a count value for each feature in the .BAM file.
Directories and file patterns will need to be changed, but the skeleton should suffice.
Note: A .GFF file must be provided. In the case of counting reads per each contig, a .GFF file will
need to be created that has entries corresponding to each contig. If a different feature in going
to be used to calculate counts on, a .GFF file must be used where EVERY ENTRY has the field specified.
It might be necessary to filter this .GFF file for only entries containing the feature of choice if
annotations are not uniform across all contigs (for example, tigrfam_acc may be annotated in some contigs
and not in others due to no matches).
Also note that htseq-count - when using "module load HTSeqcount" on Snyder - is outdated and defaults
back to python2, a version incompatible with snakemake. Make sure to run htseq-count via installation
on pip (pip install HTSeq) or bioconda (conda install -c bioconda HTSeq) for compatibility with python3.
'''

# Below, pattern matching to get a list of wildcards needed for "rule all"
allsamples = []
for i, j in zip(range(685, 720), range(1, 36)):
    one = f"030{str(i)}_{str(j)}_S{str(j)}_R1"
    two = f"030{str(i)}_{str(j)}_S{str(j)}_R2"
    allsamples.extend([one, two])
allsamples.remove("030691_7_S7_R1")


rule all:
    input:
        expand("Counts/Supernatant/{sample}.count", sample=allsamples)
        expand("Counts/SupernatantGenes/{sample}.processed.cpm.txt", sample=allsamples)

# Run counts on contigs.
rule run_htseqcount:
    input:
        cursample="supernatant_normal_contigs_500_maxbin/{sample}_filtered_sort.bam"
    output:
        curoutput="Counts/Supernatant/{sample}.count"
    shell:
        """
        htseq-count -f bam -r pos -t contig --stranded=no {input.cursample} supernatant.gff > {output.curoutput}
        """

# Run counts on specific feature.
# Note: If feature is not present in all entries, must filter first.
rule run_gene_htseqcount:
    input:
        cursample="supernatant_normal_contigs_500_maxbin/{sample}_filtered_sort.bam"
    conda:
        "envs/py27.yml"
    params:
        attrib="tigrfam_acc",
        threads=20
    output:
        curoutput="Counts/SupernatantGenes/{sample}.count"
    shell:
        """
        htseq-count -f bam -r pos -t CDS -i {params.attrib} -n {params.threads} --stranded=no {input.cursample} GFFs/supernatantAll.tigrfam_acc.gff > {output.curoutput}
        """

# From raw counts, calculate CPM values. 
rule calc_TPM:
    input:
        cursample="Counts/SupernatantGenes/{sample}.count"
    params:
        gff="GFFs/supernatantAll.tigrfam_acc.gff",
        fcol="8",
        fname="tigrfam_acc"
    output:
        out="Counts/SupernatantGenes/{sample}.processed.cpm.txt"
    shell:
        """
        python convert_cpms.py {params.gff} {params.fcol} {input.cursample} {params.fname}
        """
        
rule calc_mean_CPM:
    input:
        forward="Counts/SupernatantGenes/{sample}_R1.processed.cpm.txt",
        reverse="Counts/SupernatantGenes/{sample}_R2.processed.cpm.txt"
    output:
        "Counts/SupernatantGenes/Means/{sample}.processed.cpm.mean.txt"
    shell:
        """
        python calcMeanCounts.py -f {input.forward} -r {input.reverse} -o {output}
        """
