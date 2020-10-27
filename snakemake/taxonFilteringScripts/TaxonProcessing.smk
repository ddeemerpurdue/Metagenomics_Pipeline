'''
Author: Dane Deemer
Rules to filter bins based on Cat/Bat taxonomic annotations.
'''


# Add and remove contigs based on taxonomies of contigs and respective bins
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