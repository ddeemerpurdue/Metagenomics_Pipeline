'''
Script that integrates into the binProcessing.smk pipeline. This works with blasting
both binners and non-binners.

Example usage:
$ python blastAssemblies.py -q queryFastaFile.fasta -f 6
'''
import os
import sys
from subprocess import call
import time


def run_blast(output, full_reference, query, log, fmt):
    '''
    Description
    '''
    make_directory_command = f"mkdir -p {output.split('/blastn')[0]}"
    os.system(make_directory_command)
    download_command = f"blastn -query {query} -subject {full_reference} -outfmt {str(fmt)} -out {output} -max_target_seqs 3"
    print(download_command)
    log.write("Attempting to blast...\n")
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
            log.write(
                f"BlastN of {query} vs {full_reference} returned {retcode}\n")
    except OSError as e:
        # Want a message in stdout
        print("Execution failed: ", e, file=sys.stderr)
        log.write(
            f"BlastN of {query} vs {full_reference} triggered {e}.\n")
        return 1
    return 0


def blast_bin_files(query, fmt, log):
    '''
    Description
    '''
    now = time.localtime(time.time())
    blast = False
    sample = query.split('/')[1]
    processing = query.split('/')[2]
    bin_number = os.path.basename(query).split('.')[1]
    reference_path = f"References/{sample}/{processing}/"
    with open(log, 'w') as lg:
        lg.write(f"Time started: {time.asctime(now)}\n")
        lg.write(f"Attempting to blast {query} against a reference.\n")
        # Look to see if there is a query file:
        for file in os.listdir(reference_path):
            if file.startswith(f"{bin_number}."):
                blast = True
                lg.write(f"File {file} matches query bin id.\n")
                # accession = file.split('.')[1]
                full_reference = os.path.join(reference_path, file)
                break
            else:
                pass
        # At this point, if there was a reference then we have all information
        # we need to run a blast search
        output = f"BlastBinners/{sample}/{processing}/blastn.{bin_number}.txt"
        print(f"{output}\n{full_reference}\n{query}\n{lg.name}")
        if blast:
            if run_blast(output, full_reference, query, lg, fmt) == 1:
                return 1
        else:
            lg.write(
                f"No reference found for {query} - no blasting performed.\n")
            # We still want to write the output file so snakemake knows
            touch_command = f"touch {output}"
            os.system(touch_command)
    # Return a clean exit status
        now = time.localtime(time.time())
        lg.write(f"Time ended: {time.asctime(now)}\n")
    return 0


def blast_nobin_files(query, fmt, log):
    '''
    Description
    '''
    now = time.localtime(time.time())
    blast = False
    sample = query.split('/')[1]
    processing = query.split('/')[2]
    bin_number = os.path.basename(query).split('.')[1]
    assert bin_number == "NoBin", "Input wrong fasta file!"
    reference_path = f"References/{sample}/{processing}/"
    with open(log, 'w') as lg:
        lg.write(f"Time started: {time.asctime(now)}\n")
        lg.write(f"Attempting to blast {query} against all references.\n")
        # Look to see if there is a query file:
        for file in os.listdir(reference_path):
            print(file)
            if (file.endswith(('.fasta', '.fa')) and not file.startswith('NoBin')):
                blast = True
                lg.write(f"File {file} will be blasted against {file}.\n")
                ref_number = file.split('.')[0]
                # accession = file.split('.')[1]
                full_reference = os.path.join(reference_path, file)
            else:
                blast = False
            # At this point, loop and blast
            output = f"BlastBinners/{sample}/{processing}/blastn.NoBin_{ref_number}.txt"
            print(output)
            if blast:
                if run_blast(output, full_reference, query, lg, fmt) == 1:
                    return 1
            else:
                lg.write(
                    f"No reference found for {query} - no blasting performed.\n")
                # We still want to write the output file so snakemake knows
                touch_command = f"touch {output}"
                os.system(touch_command)
    # Return a clean exit status
        now = time.localtime(time.time())
        lg.write(f"Time ended: {time.asctime(now)}\n")
    # Since everything went correct, let's touch a token file:
    touch_command = f"touch BlastBinners/{sample}/{processing}/blastn.NoBins.success.txt"
    os.system(touch_command)
    return 0


if __name__ == "__main__":
    import argparse
    ''' ARGUMENTS '''
    parser = argparse.ArgumentParser(description="Parser")
    parser.add_argument("-q", "--Query", help="Blast query file",
                        required=True)
    parser.add_argument("-f", "--Format", help="Blast file format",
                        required=False, default=6)
    parser.add_argument("-l", "--Log", help="Log file",
                        required=True)
    parser.add_argument("-n", "--NoBin", help="If used, blast non-binners",
                        required=False, default=False, action='store_true')
    argument = parser.parse_args()
    if argument.NoBin:
        blast_nobin_files(argument.Query, argument.Format, argument.Log)
    else:
        blast_bin_files(argument.Query, argument.Format, argument.Log)
