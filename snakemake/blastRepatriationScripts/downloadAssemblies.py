'''
Author: Dane
Date: 21Sep20
Purpose: Given the output from writeModeGffFeaturePerBin.py, a file with 3
tab-delimited fields: [bin_number, assembly_id, and count], download each
assembly using the ./datasets command line tool provided by ncbi.
To make this compatible with snakemake and reduce the amount of downloading
needed, this will first check if an assembly is not already in a database
that we've previously created, and will download if necessary. The file used
that links bins to assemblies will be used to create a dictionary in the form:
dict[bin_number] = assembly.

Example usage:
$ python downloadAssemblies.py
NOTE:
Unfortunately, ./datasets cannot find all assemblies at this time
'''
import os
import sys

print(os.getcwd())


def read_mode_assembly_file(modeGffFile):
    '''
    File to read in the output from writeModeGffFeaturePerBin.py into a
    dictionary in the form:
    dict[bin_number] = assembly
    '''
    bin_dic = {}
    with open(modeGffFile) as i:
        line = i.readline().strip()
        assert len(line.split('\t')) == 3, 'Invalid number of fields!'
        while line:
            line = line.split('\t')
            bin_number = line[0]
            assembly = line[1]
            bin_dic[bin_number] = assembly
            line = i.readline().strip()
    return bin_dic


def check_assembly_downloads(
        modeGffFile,
        assemblyLocation='Assemblies/ncbi_dataset/data/'
):
    '''
    Function to check if assembly is already downloaded, and if not then
    it downloads it and adds to the standard directory
    '''
    assemblies_to_download = []
    bin_dic = read_mode_assembly_file(modeGffFile)
    for bin_num, assembly in bin_dic.items():
        if assembly in os.listdir(assemblyLocation):
            pass
        else:
            assemblies_to_download.append(assembly)
    print(f"There are {len(assemblies_to_download)} assemblies to download.")

    # Feed all assemblies into a comma-delimited string to prep for command
    assembly_string = ','.join(assemblies_to_download)
    download_command = f"./datasets download assembly {assembly_string} --filename tmpAssemblies.zip"
    unzip_command = f"unzip tmpAssemblies.zip -d tmpAssemblies"
    move_command = f"mv tmpAssemblies/ncbi_dataset/data/G* {assemblyLocation}"
    os.system(download_command)
    os.system(unzip_command)
    os.system(unzip_command)
    return 0


if __name__ == "__main__":
    check_assembly_downloads(sys.argv[1])
