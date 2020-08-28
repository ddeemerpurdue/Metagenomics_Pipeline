
'''
Snakemake plugin that runs NucDiff on all possible bin comparisons between two samples.
The two samples can be different, or can be the same sample to for intra-comparisons.
This is the same basic skeleton as the DnaDiff.smk file.
'''

# Change bins to whatever appropriate values are needed #
refbins = ["{0:03}".format(i) for i in range(1, 255)]
querbins = ["{0:03}".format(i) for i in range(1, 235)]
# Change bins to whatever appropriate values are needed #

rule all:
    input:
        expand("NucDiffOutput/R{r}.{rbin}_Q{q}.{qbin}.report", r="Supernatant", q="Particle", rbin=refbins, qbin=querbins)


rule run_dnadiff:
    input:
        reference="{path}/{bin}.fasta",
        query="{path2}/{bin2}.fasta"
    params:
        prefix="R{path}{bin}_Q{path2}{bin2}"
    output:
        directory("NucDiffOutput")
    shell:
        "nucdiff {input.reference} {input.query} {output} {params.prefix}"
