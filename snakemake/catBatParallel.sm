'''
Author: Dane
Purpose:
To parallelize the process of running CATBAT and split up each job by bin. For example,
if an experiment has 200 bins, the user can run each bin as a separate job, with the only
change to this script being the "all_samples" pattern, along with directory information.
'''


nSamples = 234
allsamples = ["{0:03}".format(i) for i in range(1, len(nSamples) + 1)]

rule all:
    input:
        expand("PartSM.Bin{number}.MAG.ORF2LCA.txt", number=allsamples)

rule run_bat:
    input:
        samples="particle_bins/Bin.{num}.fasta"
    output:
        out="PartSM.Bin{num}.MAG.ORF2LCA.txt"
    shell:
        "/scratch/snyder/d/ddeemer/WhiteRed/scgs/catbat/CAT/CAT_pack/CAT bin -b {input.samples} \
        -d /scratch/snyder/d/ddeemer/WhiteRed/scgs/catbat/2020-05-22_CAT_database/ \ 
        -t /scratch/snyder/d/ddeemer/WhiteRed/scgs/catbat/2020-05-22_taxonomy -o PartSM.Bin{wildcards.num}.MAG"
