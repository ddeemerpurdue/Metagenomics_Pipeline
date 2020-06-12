    Folder containing conda environment files to be used for local rules.
    
Guidelines:

1.) Each env file should be specific to one rule unless multiple rules share the
same dependencies.

2.) Environment files should only contain the required packages to run that rule
in order to reduce waiting times for loading multiple unneccessary packages and to
avoid potential package disagreements.

3.) If the rule only requires base packages, omit the conda section of snakemake
altogether and rule without loading any extra packages.
