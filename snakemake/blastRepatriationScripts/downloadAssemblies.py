'''
Author: Dane
Date: 01Oct20
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
from subprocess import call
import subprocess


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
        logfile,
        assemblyLocation='/depot/lindems/data/Dane/Assemblies/'
):
    '''
    Function to check if assembly is already downloaded, and if not then
    it downloads it and adds to the standard directory
    '''
    assemblies_to_download = {}
    bin_dic = read_mode_assembly_file(modeGffFile)
    try:
        os.makedirs(assemblyLocation)
    except FileExistsError:
        pass
    with open(logfile, 'a') as log:
        # No need to download if we've stumbled across this before
        for bin_num, assembly in bin_dic.items():
            assembly_match = f"{bin_num}.{assembly.replace('.', '_')}.fasta"
            # If the assembly already exists in the database folder
            if (assembly_match in os.listdir(assemblyLocation)):
                full_path = os.path.join(assemblyLocation, assembly_match)
                # And if that file is not empty (0 bytes)
                if os.stat(full_path).st_size != 0:
                    log.write(
                        f"Bin {bin_num}'s assembly <{assembly}> is already present.\n")
                else:
                    log.write(
                        f"Bin {bin_num}'s assembly <{assembly}> is empty and needs downloading.\n")
                    assemblies_to_download[bin_num] = assembly
            else:  # Otherwise, add it to a list of assemblies to download
                assemblies_to_download[bin_num] = assembly
                log.write(
                    f"Bin {bin_num}'s assembly <{assembly}> needs downloading.\n")
        log.write(
            f"{len(assemblies_to_download.keys())} assemblies need to be downloaded.\n\n")

        # Loop assemblies and download what is not present or empty.
        successes = 0
        for bin_num, assembly in assemblies_to_download.items():
            # Don't like having ".1.fasta", change to "_1.fasta"
            write_assembly = f"{bin_num}.{assembly.replace('.', '_')}"
            full_assembly = os.path.join(assemblyLocation, write_assembly)
            if assembly.startswith("GCF"):
                gca = False
                download_command = f"esearch -db nucleotide -query {assembly} | efetch -format fasta > {full_assembly}.fasta > /dev/null 2>&1"
            elif assembly.startswith("GCA"):
                gca = True
                command = f"esearch -db assembly -query GCA_900066305.1 | esummary | xtract -pattern DocumentSummary -element FtpPath_GenBank"
                get_download_link = subprocess.Popen(
                    command, stdout=subprocess.PIPE, shell=True)
                ftp_link = get_download_link.stdout.read().strip().decode("utf-8")
                file_link = f"{str(os.path.basename(ftp_link))}_genomic.fna.gz"
                full_link = f"{ftp_link}/{file_link}"
                print(full_assembly)
                download_command = f"wget {full_link} -O {full_assembly}.fasta.gz && {full_assembly}.fasta.gz > /dev/null 2>&1"
            log.write(
                f"Command to download {assembly}:\n{download_command}\n\n")
            try:
                # Execute the command
                retcode = call(download_command, shell=True)
                if retcode > 0:  # Any error, let's catch it
                    # For errors like this, want it sent to the stderr too
                    print("Process was terminated by signal",
                          retcode, file=sys.stderr)
                    log.write(f"Process was terminated by signal {retcode}\n")
                    return 1  # If return 1, we won't write the token
                    # output file, halting the pipeline.
                else:
                    pass
                    log.write(f"Download of {assembly} returned {retcode}\n")
                    successes += 1

            except OSError as e:
                # Want a message in stdout
                print("Execution failed: ", e, file=sys.stderr)
                log.write(f"Download of {assembly} triggered {e}.\n")
                return 1
        log.write(f"{successes} assemblies successfully downloaded.\n")
    return 0


if __name__ == "__main__":
    import argparse
    import time
    now = time.localtime(time.time())
    """ Arguments """
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-g", "--GffMode",
                        help="Mode GFF file",
                        required=True)
    parser.add_argument("-a", "--AssemblyLocation",
                        help="Location containing all assembly files",
                        required=False,
                        default='/depot/lindems/data/Dane/Assemblies/')
    parser.add_argument("-l", "--Log",
                        help="Log file to write to",
                        required=True)
    argument = parser.parse_args()
    with open(argument.Log, 'a') as log:
        log.write(f"Time started: {time.asctime(now)}\n")
    run = check_assembly_downloads(argument.GffMode, argument.Log,
                                   argument.AssemblyLocation)
    if run == 0:  # If everything runs correctly, create a token file.
        command = f"touch {str(argument.GffMode).replace('.txt', '.success.txt')}"
        os.system(command)
    now = time.localtime(time.time())
    with open(argument.Log, 'a') as log:
        log.write(f"\nTime ended: {time.asctime(now)}\n")
