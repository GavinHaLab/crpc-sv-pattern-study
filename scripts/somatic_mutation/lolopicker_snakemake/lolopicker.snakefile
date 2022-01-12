#lolopicker.snakefile

#Ha Lab
#Fred Hutchinson Cancer Research Center

"""
#before running snakemake, do in tmux terminal:
ml snakemake/5.8.1-foss-2016b-Python-3.7.4
ml Python/3.7.4-foss-2016b-fh1
ml Pysamstats/1.1.2-foss-2016b-Python-3.7.4
ml BEDOPS/2.4.35-foss-2016b

#command to run snakemake (remove -np end when done validating):
snakemake -s lolopicker.snakefile --latency-wait 60 --keep-going --restart-times 2 --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" -j 50 -np
"""

configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
    input:
        expand("results/bed_files/{tumors}/whole_genome_no_indels.bed", tumors = config["samples"]),
        expand("results/bed_files/{tumors}/whole_genome_snvs_only.bed", tumors = config["samples"]),
        expand("results/bed_files/{tumors}/{chromosomes}_snvs.bed", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/{chromosomes}/raw_somatic_variants.txt", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/{chromosomes}/control_stats.txt", tumors = config["samples"], chromosomes = config["chromosomes"]),
        expand("results/{tumors}/merged_control_stats.txt", tumors = config["samples"]),
        expand("results/{tumors}/reject_calls.txt", tumors = config["samples"]),
        expand("results/{tumors}/stats_calls.txt", tumors = config["samples"])

"""
rule format_bed_files:
    input:
        expand(config["mutect_results_path"] + "/{tumors}/filtered_all.vcf.gz", tumors = config["samples"])
    output:
        whole_genome_no_indels_bed = protected("results/bed_files/{tumors}/whole_genome_no_indels.bed"),
        whole_genome_snvs_only_bed = protected("results/bed_files/{tumors}/whole_genome_snvs_only.bed"),
        chr_bed_files = protected("results/bed_files/{tumors}/{chromosomes}_snvs.bed")
    params:
        mutect_results_path = config["mutect_results_path"],
        bed_formatting_script = config["bed_formatting_script"]
    shell:
        "zcat {input} | vcf2bed --snvs > {output.whole_genome_no_indels_bed}"
        "python {params.bed_formatting_script} results/bed_files"
        "bedextract {chromosomes} {output.whole_genome_snvs_only_bed} > results/bed_files/{tumors}/{chromsomes}_snvs.bed"
"""

rule mutect2_vcf_to_no_indel_bed:
    input:
        config["mutect_results_path"] + "/{tumors}/filtered_all.vcf.gz"
    output:
        protected("results/bed_files/{tumors}/whole_genome_no_indels.bed")
    shell:
        "zcat {input} | grep 'PASS' | vcf2bed --snvs > {output}"

rule make_snv_only_bed:
    input:
        "results/bed_files/{tumors}/whole_genome_no_indels.bed"
    output:
        protected("results/bed_files/{tumors}/whole_genome_snvs_only.bed")
    run:
        no_indels_bed_file = input[0]
        #make filename of new bed file
        snvs_only_bed = no_indels_bed_file.replace("_no_indels.bed", "_snvs_only.bed")
        with open(no_indels_bed_file, "r") as input_bed:
            with open(snvs_only_bed, "w") as output_bed:
                #for every line in the input bed
                for line in input_bed.readlines():
                    #if the alt allele (col 7) is only one letter (aka event is SNV)
                    if len(line.split("\t")[6]) == 1:
                        #write that line to the output bed
                        output_bed.write(line)

rule split_bed_by_chr:
    input:
        "results/bed_files/{tumors}/whole_genome_snvs_only.bed"
    output:
        protected("results/bed_files/{tumors}/{chromosomes}_snvs.bed")
    shell:
        "bedextract {wildcards.chromosomes} {input} > {output}"

rule call_variants:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumors][0],
        normal_filepath = lambda wildcards: config["samples"][wildcards.tumors][2],
        bed = "results/bed_files/{tumors}/{chromosomes}_snvs.bed"
    output:
        protected("results/{tumors}/{chromosomes}/raw_somatic_variants.txt")
    params:
        reference_genome = config["reference_genome"],
        somatic_script = config["somatic_script"],
        phasing_mode = config["phasing_mode"]
    log:
        "logs/call_variants/{tumors}_{chromosomes}_call_variants.txt"
    shell:
        "(python {params.somatic_script} \
        -t {input.tumor_filepath} \
        -n {input.normal_filepath} \
        -r {params.reference_genome} \
        -b {input.bed} \
        --phasing_mode {params.phasing_mode} \
        -o results/{wildcards.tumors}/{wildcards.chromosomes}) &> {log}"

rule inspect_controls:
    input:
        "results/{tumors}/{chromosomes}/raw_somatic_variants.txt"
    output:
        protected("results/{tumors}/{chromosomes}/control_stats.txt")
    params:
        control_script = config["control_script"],
        control_tsv_file = config["control_tsv_file"],
        reference_genome = config["reference_genome"]
    log:
        "logs/inspect_controls/{tumors}_{chromosomes}_inspect_controls.txt"
    shell:
        "(python {params.control_script} \
        -l {params.control_tsv_file} \
        -r {params.reference_genome} \
        -n 4 \
        -o results/{wildcards.tumors}/{wildcards.chromosomes}) &> {log}"

rule merge_control_stats:
    input:
        chr1_control_stats = "results/{tumors}/chr1/control_stats.txt",
        chr2_control_stats = "results/{tumors}/chr2/control_stats.txt",
        chr3_control_stats = "results/{tumors}/chr3/control_stats.txt",
        chr4_control_stats = "results/{tumors}/chr4/control_stats.txt",
        chr5_control_stats = "results/{tumors}/chr5/control_stats.txt",
        chr6_control_stats = "results/{tumors}/chr6/control_stats.txt",
        chr7_control_stats = "results/{tumors}/chr7/control_stats.txt",
        chr8_control_stats = "results/{tumors}/chr8/control_stats.txt",
        chr9_control_stats = "results/{tumors}/chr9/control_stats.txt",
        chr10_control_stats = "results/{tumors}/chr10/control_stats.txt",
        chr11_control_stats = "results/{tumors}/chr11/control_stats.txt",
        chr12_control_stats = "results/{tumors}/chr12/control_stats.txt",
        chr13_control_stats = "results/{tumors}/chr13/control_stats.txt",
        chr14_control_stats = "results/{tumors}/chr14/control_stats.txt",
        chr15_control_stats = "results/{tumors}/chr15/control_stats.txt",
        chr16_control_stats = "results/{tumors}/chr16/control_stats.txt",
        chr17_control_stats = "results/{tumors}/chr17/control_stats.txt",
        chr18_control_stats = "results/{tumors}/chr18/control_stats.txt",
        chr19_control_stats = "results/{tumors}/chr19/control_stats.txt",
        chr20_control_stats = "results/{tumors}/chr20/control_stats.txt",
        chr21_control_stats = "results/{tumors}/chr21/control_stats.txt",
        chr22_control_stats = "results/{tumors}/chr22/control_stats.txt",
        chrX_control_stats = "results/{tumors}/chrX/control_stats.txt",
        chrY_control_stats = "results/{tumors}/chrY/control_stats.txt"
    output:
        protected("results/{tumors}/merged_control_stats.txt")
    log:
        "logs/merge_control_stats/{tumors}_merge_control_stats.txt"
    shell:
        "(find results/{wildcards.tumors}/* -name control_stats.txt | xargs -n 1 tail -n +2 > {output}) 2> {log}"

rule core_stats:
    input:
        "results/{tumors}/merged_control_stats.txt"
    output:
        protected("results/{tumors}/reject_calls.txt"),
        protected("results/{tumors}/stats_calls.txt")
    log:
        "logs/core_stats/{tumors}_core_stats.txt"
    params:
        stats_script = config["stats_script"]
    shell:
        "(python {params.stats_script} \
        -o results/{wildcards.tumors} \
        --method FDR) &> {log}"
