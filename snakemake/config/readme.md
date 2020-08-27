    Folder containing configuration files for snakemake pipeline.

config.baseMetagenomics.yaml contains a running list of all config attributes used.

Guidelines:

- Any 'soft-code' that is intended to be changed on a per-run basis should
be at the top of each Step

- Follow pattern-matching schemes to make the config file dynamic

- Avoid redundancies
