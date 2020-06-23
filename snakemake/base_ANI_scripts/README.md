# Location for all the script files that are associated with base_ANI.sm

The purpose of base_ANI.sm is to start with multiple multi-fasta files
and do a pairwise comparison across all files and all sequences in each
file in an efficient manner. The process works by parallelizing each
sample-sample comparison, and then further breaking down each of those
comparisons into 10 parts to be run in parallel. The product can be submitted
to a cluster to greatly increase the efficiency of obtaining ANI's for each
sequence.
