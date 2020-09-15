#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --time=04:00:00
#SBATCH --ntasks=20
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=roles@purdue.edu

#This job submission takes the txt output from the python file consensus_add_gff.py 
#The input file should have three columns: accession, number, and bin number
#The function is to download the genomedb for each accession number
#Then to run mummer on each one 
# Need to change directories before each use.

# Step 1: Load required modules:
module load bioinfo
module load edirect

cd /scratch/snyder/r/roles/ash

#Step 2: Define variables
{
read
while IFS= read -r line
	do
		BIN=($(echo "$line" | awk -F "\t" '{ print $1 }'))
		ACC=($(echo "$line" | awk -F "\t" '{ print $2 }'))
		esearch -db nucleotide -query $ACC < /dev/null | efetch -format fasta > particle/$BIN.fasta
		echo "$BIN"
		echo "$ACC"
    
    done 
}< "/scratch/snyder/r/roles/ash/particleAssemblies_GCF.txt"

mapfile -t ACCS < <(awk '{print $2}' /scratch/snyder/r/roles/ash/particleAssemblies_GCA.txt)
mapfile -t BINS < <(awk '{print $1}' /scratch/snyder/r/roles/ash/particleAssemblies_GCA.txt)

#while IFS= read -r line; do
for i in "${!ACCS[@]}"; do
		#BIN=($(echo "$line" | awk -F "\t" '{ print $1 }'))
		#ACC=($(echo "$line" | awk -F "\t" '{ print $2 }'))
		ACC="${ACCS[i]}"
		BIN="${BINS[i]}"
		esearch -db assembly -query $ACC \
			| esummary \
		    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
		    | while read -r line ; 
			do
			fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
			wget -O $BIN.fasta.gz "$line/$fname";
			gunzip $BIN.fasta.gz
			done
		echo "$BIN"
		echo "$ACC"
done
