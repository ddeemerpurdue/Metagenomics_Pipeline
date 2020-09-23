# 3 Step Bin Processing

### Based on:
1. Taxonomy
2. Average nucleotide identity
3. Reference assembly

## Introduction
This 3-phase processing pipeline is designed to both add and remove contigs from bins based on different reference annotations. Three main programs are utilized during this pipeline and correspond to a step: i.) CatBat ii.) fastANI and iii.) MetaErg. Step 1 both removes and add contigs, whereas steps 2 and 3 are limited to only adding contigs into bins for now. The end results are bin identification files and fasta files. The purpose of this is to be able to automatically run the processing steps with various combinations of parameters to see how the final bin product compares the to original results.

### Initial Setup
A specific directory structure is required to make sure the pipeline runs correctly. First, create a directory for the entire project.
$ mkdir MyProject
$ cd MyProject
Note: From now on, reference to '~' corresponds to the /MyProject/ folder.
Three main folders are required at the top level:
1. config
- This contains a config.yaml file, which is specifies the project-specific variables.
2. input
- All files required for this analysis.
3. workflow
- This is where all results will be saved for the project.

The initial tree structure should look as follows:

./MyProject
+-- config
+-- input
+-- workflow

#### Configuration directory:
Inside this directory there should be 2 files:
1. config.yaml
This species multiple variables needed for the pipeline. Below is an example:
---
\#config/config.yaml  

\# General attributes:  
samples: ['particle', 'supernatant']  

\# (1) Taxonomy Based Processing  
TaxonAddThresh: [90, 95, 99]  
TaxonRemoveThresh: [90, 95, 99]  

\# (2) ANI Based Processing  
ANIAssemblyFilterSize: [2000, 5000]  
ANIAssemblySplitSize: 9  
ANIAssemblySplits: ['aa','ab','ac','ad','ae','af','ag','ah','ai','aj']  
Note: split size will normally split file into N + 1 files (if modulo != 0)  
FastANIMinFraction: 0.2  
FastANIFragLength: 1000  
ANIRepatIdentThreshold: 95  
ANIRepatCountThreshold: 20  
---

2. cluster.json
---
\#config/cluster.json  

{  
    "\__default__":  
    {  
        "account": "standby",  
        "mem": "20G",  
        "time": "04:00:00",  
        "cpus": 20,  
        "ntasks-per-node": 20,  
        "nodes": 1  
    }  
}  
---


#### Input directory:
./input  
+-- Assembly  
    +-- sample1.assembly.fasta  
    +-- sample2.assembly.fasta  
+-- OriginalBins  
    +-- {sample1}  
        +-- Bin.{number}.fasta  
        +-- Bin.{number}.fasta  
        +-- etc.
    +-- {sample2) 
        +-- Bin.{number}.fasta  
        +-- Bin.{number}.fasta  
        +-- etc. 
+-- Cat  
    +-- {sample1}/{sample1}.C2C.names.txt  
    +-- {sample2}/{sample1}.C2C.names.txt  
+-- Bat  
    +-- {sample1}/{sample1}.Bin2C.names.txt  
    +-- {sample2}/{sample1}.Bin2C.names.txt  
+-- GFF  
    +-- {sample1}/{sample1}.All.gff  
    +-- {sample2}/{sample1}.All.gff  

Note: The Cat and Bat directories correspond to output files from the program CatBat. When copying into these folders, most likely the names will need to change to comply with this pipelines rules.  
For the GFF file, this file must contain an attribute with the name *genomedb_acc* in order for this to work. MetaErg provides a gff file with this annotation.  

***