    Folder containing conda environment files to be used for local rules.

## allEnvironments.yaml contains a running list of all programs and version to choose from
- If a new program and/or version needs to be used, please consider the following:
1. If it is available via conda, add it to the allEnvironments.yaml file so all collaborators
are using the same version.
2. If it is not available via conda, look for availability on Purdue's Snyder and add it to the
moduleList.txt file.
3. If the program is not available on either platform, provide a link to the installation instructions
in this folder or provide your own installation instructions. The next step would be to reach out to
Purdue's bioinformatics department to have them install the module within Snyder.

Guidelines:

1.) Each env file should be specific to one rule unless multiple rules share the
same dependencies.

2.) Environment files should only contain the required packages to run that rule
in order to reduce waiting times for loading multiple unneccessary packages and to
avoid potential package disagreements.

3.) If the rule only requires base packages, omit the conda section of snakemake
altogether and rule without loading any extra packages.
