'''
Snakemake plugin that runs nucmer on all possible bin comparisons between two samples.
The two samples can be different, or can be the same sample to for intra-comparisons
'''

refbins = ["{0:03}".format(i) for i in range(1, 255)]
querbins = ["{0:03}".format(i) for i in range(1, 235)]

rule all:
    input:
        expand("DnaDiffOutput/R{r}.{rbin}_Q{q}.{qbin}.report", r="Supernatant", q="Particle", rbin=refbins, qbin=querbins)


rule run_dnadiff:
    input:
        reference="{path}/{bin}.fasta",
        query="{path2}/{bin2}.fasta"
    params:
        prefix="DnaDiffOutput/R{path}.{bin}_Q{path2}.{bin2}"
    output:
        "DnaDiffOutput/R{path}.{bin}_Q{path2}.{bin2}.report"
    shell:
        "dnadiff {input.reference} {input.query} --prefix {params.prefix}"
