'''
Generic snakemake file for MetaErg
'''

'''
SANDBOX
Bin path: /depot/lindems/data/Dane/Shared/AshProject/20Aug20Merging/Particle
Modules to load:
ml bioinfo
ml perl/5.26.1
ml java/11.0.2
ml minced/0.2.0
ml ARAGORN/1.2.38
ml HMMER/3.2.1
ml blast/2.9.0+
ml diamond/0.9.26
ml prodigal/2.60
ml MinPath/1.2
ml signalp/4.1c
ml tmhmm/2.0c

perl ../../../local/metaerg/metaerg/bin/metaerg.pl --dbdir $RCAC_SCRATCH/local/metaerg/metaerg/db --outdir $d  --prefix Supernatant $f

'''


rule run_metaerg:
    input:
        assembly="{sample_num}.fasta"
    params:
        executable="/scratch/snyder/d/ddeemer/local/metaerg/metaerg/bin/metaerg.pl",
        db="$RCAC_SCRATCH/local/metaerg/metaerg/db",
        outdir="{sample_num}_MetaErg",
        prefix="MetaErg_{sample_num}"
    log:
        logs/{sample}.metaerg.log
    output:
        {sample_num}_MetaErg/{params.prefix}_MetaErg.fna
    shell:
        """
        ml bioinfo perl/5.26.1 java/11.0.2 minced/0.2.0 ARAGORN/1.2.38 HMMER/3.2.1 blast/2.9.0+ diamond/0.9.26 prodigal/2.60 MinPath/1.2 signalp/4.1c tmhmm/2.0c
        source ~/.sm_base
        perl {params.executable} --dbdir {params.db} --outdir {params.outdir} --prefix {params.prefix} {input.assembly} > {log}
        """
