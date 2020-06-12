configfile: "config.yaml"

rule tipp:
	input:
		samples=lambda wildcards: f"{config['samplepath']}{wildcards.samples}.contigs.fasta"
	output:
    		"Metaphlan/profiled_metagenome.{samples}.txt"
	log:
    		"log/{samples}.log"
	conda:
		"mpa.yaml"
	params:
    		input_type="fasta"
	shell:
    		"metaphlan {input.samples} --nproc 20 --input_type {params.input_type} -o {output}"