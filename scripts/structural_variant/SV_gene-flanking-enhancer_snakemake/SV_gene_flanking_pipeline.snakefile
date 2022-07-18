# modules to load
# 
# ml Python/3.7.4-foss-2019b-fh1
# ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
# ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4
# ml snakemake/5.19.2-foss-2019b-Python-3.7.4
# ml R/3.6.2-foss-2019b-fh1
# ml BEDTools
#
#To run this snakemake on the cluster (remove -np): 
#snakemake -s SV_gene_flanking_pipeline.snakefile --latency-wait 60 --restart-times 3 --keep-going --cluster-config config/cluster_slurm.yaml --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output} --constraint=gizmok" -j 10 -np

configfile: "config/config.yaml"

rule all:
    input: 
        #expand("results/prepare_bed_and_intervals/{bedpe_prefix}_with_svID.bedpe", bedpe_prefix=config["bedpe_prefix"]),
        "results/prepare_bed_and_intervals/all_sv_events_breakpoints_only.txt",
    	"results/intersect_files/intersect_sv_and_gene_body.txt",
    	"results/intersect_files/short_sv_minus_gene_body_transection.txt",
    	"results/intersect_files/intersect_sv_and_flanking_up.txt",
    	"results/intersect_files/intersect_sv_and_flanking_dn.txt",
    	"results/comut_files/comut_intersect_sv_and_gene_body.csv",
    	"results/comut_files/comut_intersect_sv_and_flanking_up.csv",
    	"results/comut_files/comut_intersect_sv_and_flanking_dn.csv",
    	expand("results/longSV_analysis/{outPrefix}_longSV_summary_gene_info_single.txt", outPrefix=config["outPrefix"]),
    	expand("results/shortSV_enh_analysis/{outPrefix}_enhUpMatrix_H3K27ac_single.txt", outPrefix=config["outPrefix"]+"_shortSV"),
    	expand("results/shortSV_enh_analysis/{outPrefix}_enhDownMatrix_H3K27ac_single.txt", outPrefix=config["outPrefix"]+"_shortSV"),
    	expand("results/longSV_analysis/{outPrefix}_enhUpMatrix_H3K27ac_single.txt", outPrefix=config["outPrefix"]+"_longSV"),
    	expand("results/longSV_analysis/{outPrefix}_enhDownMatrix_H3K27ac_single.txt", outPrefix=config["outPrefix"]+"_longSV"),
    	expand("results/longSV_analysis/{outPrefix}_summary_H3K27ac_single.bedpe", outPrefix=config["outPrefix"]+"_longSV"),
    	"results/longSV_analysis/all_longSV_flanking_events_with_enhancer_single_summary_full_version.txt",
    	"results/shortSV_enh_analysis/shortSV_enh_analysis_complete.txt",
    	"results/comut_files/all_comut_SV_combined_single.csv",

#1. Prepare sv (breakpoints only in bed format) file and gene/flanking interval file
# This requires sv bedpe file, sample_list files and gene lists with coordinates information.

rule prepare_bed_and_intervals:
	input:
		bedpe=expand("{bedpe_prefix}.bedpe", bedpe_prefix=config["bedpe_prefix"]),
		samples=lambda wildcards: config["sample_file"],
		gene_list=lambda wildcards: config["gene_file"],	
	output: 
		"results/prepare_bed_and_intervals/gene_body_intervals.txt",
		"results/prepare_bed_and_intervals/genes_upstream_intervals.txt",
		"results/prepare_bed_and_intervals/genes_dnstream_intervals.txt",
		expand("results/prepare_bed_and_intervals/{bedpe_prefix}_with_svID.bedpe", bedpe_prefix=config["bedpe_prefix"]),
		"results/prepare_bed_and_intervals/all_sv_events_breakpoints_only.txt"
	params:
		script="script/prepare_bed_and_gene.R",
	shell:
		"Rscript {params.script} {input.bedpe} {input.samples} {input.gene_list}"

#2. Intersect between "all SV events" and gene body (gene body transection events) 

rule intersect_SV_geneBody:
	input:
		sv_bed="results/prepare_bed_and_intervals/all_sv_events_breakpoints_only.txt",
		gene_body_intervals="results/prepare_bed_and_intervals/gene_body_intervals.txt",
	output:
		"results/intersect_files/intersect_sv_and_gene_body.txt"
	shell:
		"intersectBed -a {input.sv_bed} -b {input.gene_body_intervals} -wb > {output}"

#3. Retrieve all SVs "not" contained in "intersect_sv_and_gene_body"
# Then create separate sv bedpe (short events, long events) based on SPAN size : 
# 1) For intra short events, produce two files: (1) short version for bedtools and (2) long version for follow-up H3K27ac enhancer analysis.
# 2) For intra long events (SPAN >=10MB) and translocation events, produce one file.

rule retrieve_SV_not_in_geneBody:
	input:
		bedpe_with_svID=expand("results/prepare_bed_and_intervals/{bedpe_prefix}_with_svID.bedpe", bedpe_prefix=config["bedpe_prefix"]),
		intersect_sv_and_gene_body="results/intersect_files/intersect_sv_and_gene_body.txt",
	output:
		"results/intersect_files/short_sv_minus_gene_body_transection.txt",
		"results/intersect_files/short_sv_minus_gene_body_transection_full.txt",
		"results/intersect_files/long_sv_minus_gene_body_transection_full.txt",
	params:
		script="script/filter_gene_body_SVs.R",
	shell:
		"Rscript {params.script} {input.bedpe_with_svID} {input.intersect_sv_and_gene_body}"

#4. Map the remaining short SVs (SPAN <= 10mb) to flanking region (1mb upstream and dnstream) of genes.

rule intersect_shortSV_in_flanking:
	input:
		short_sv_minus_gene_body="results/intersect_files/short_sv_minus_gene_body_transection.txt",
		upstream_intervals="results/prepare_bed_and_intervals/genes_upstream_intervals.txt",	
		dnstream_interval="results/prepare_bed_and_intervals/genes_dnstream_intervals.txt",
	output:
		intersect_sv_and_flanking_up="results/intersect_files/intersect_sv_and_flanking_up.txt",
		intersect_sv_and_flanking_dn="results/intersect_files/intersect_sv_and_flanking_dn.txt",
	shell:
		"""
			intersectBed -a {input.short_sv_minus_gene_body} -b {input.upstream_intervals} -wb > {output.intersect_sv_and_flanking_up}
			intersectBed -a {input.short_sv_minus_gene_body} -b {input.dnstream_interval} -wb > {output.intersect_sv_and_flanking_dn}
		"""

#5. Create comut file for gene body and shortSV flanking events

rule create_comut_file_for_gene_flanking:
	input:
		"results/intersect_files/intersect_sv_and_flanking_up.txt",
		"results/intersect_files/intersect_sv_and_flanking_dn.txt",
		"results/intersect_files/intersect_sv_and_gene_body.txt",
	params:
		script="script/create_comut_from_bedintersect_files.R",
		intersect_bed_dir="results/intersect_files/",
	output:
		"results/comut_files/comut_intersect_sv_and_gene_body.csv",
		"results/comut_files/comut_intersect_sv_and_flanking_up.csv",
		"results/comut_files/comut_intersect_sv_and_flanking_dn.csv"
	shell:
		"Rscript {params.script} {params.intersect_bed_dir}"

#6. Create flanking region mapping table for long SVs (also record their partner gene information)

rule create_mappingTable_longSVevents:
	input:
		longSV_info="results/intersect_files/long_sv_minus_gene_body_transection_full.txt",
		gene_list=lambda wildcards: config["gene_file"],
		annotation=lambda wildcards: config["annotation_file"],
		chrEnd=lambda wildcards: config["chrEnd_rds"],
	output:
		"results/longSV_analysis/{outPrefix}_longSV_summary_gene_info_single.txt"
	params:
		script="script/create_mappingTable_longSVevents.R",
		outPrefix=lambda wildcards: config["outPrefix"],
	shell:
		"Rscript {params.script} {input.longSV_info} {input.gene_list} {input.annotation} {input.chrEnd} {params.outPrefix}"

# 7-1. Short SV "Enhancer" events -- confirm if there is H3K27ac peaks contained in flanking region (SV should fully contain H3k27ac peaks)

rule makeSVmatrixByGenes_shortSV_Enhancer:
	input:
		shortSV_info="results/intersect_files/short_sv_minus_gene_body_transection_full.txt",
		gene_list=lambda wildcards: config["gene_file"],
		chrEnd=lambda wildcards: config["chrEnd_rds"],
		H3K27ac_peaks=lambda wildcards: config["H3K27ac_bed"],
	output:
		"results/shortSV_enh_analysis/{outPrefix}_enhUpMatrix_H3K27ac_single.txt",
		"results/shortSV_enh_analysis/{outPrefix}_enhDownMatrix_H3K27ac_single.txt",
	params:
		script="script/makeSVmatrixByGenes_with_H3K27ac_enhancer_for_intraSV.R",
		outPrefix=lambda wildcards: config["outPrefix"],
	shell:
		"Rscript {params.script} {input.shortSV_info} {input.gene_list} {input.chrEnd} {input.H3K27ac_peaks} {params.outPrefix}_shortSV"

# 7-2. Long SV "Enhancer" events  -- same rule applies as above
rule makeSVmatrixByGenes_longSV_Enhancer:
	input:
		longSV_info="results/intersect_files/long_sv_minus_gene_body_transection_full.txt",
		gene_list=lambda wildcards: config["gene_file"],
		chrEnd=lambda wildcards: config["chrEnd_rds"],
		H3K27ac_peaks=lambda wildcards: config["H3K27ac_bed"],
	output:
		"results/longSV_analysis/{outPrefix}_enhUpMatrix_H3K27ac_single.txt",
		"results/longSV_analysis/{outPrefix}_enhDownMatrix_H3K27ac_single.txt",
		"results/longSV_analysis/{outPrefix}_summary_H3K27ac_single.bedpe"
	params:
		script="script/makeSVmatrixByGenes_with_H3K27ac_enhancer_for_interSV.R",
		outPrefix=lambda wildcards: config["outPrefix"],
	shell:
		"Rscript {params.script} {input.longSV_info} {input.gene_list} {input.chrEnd} {input.H3K27ac_peaks} {params.outPrefix}_longSV"

# 8. Create summary table for long SV events (1) flanking and (2) gain of enhancer from partner (enhancer hijacking)
# This script also generates comut files for long SV with flanking + H3K27ac enhancer information.

rule create_summary_for_longSV_flanking_Enhancer:
	input:
		longSV_summary="results/longSV_analysis/"+config["outPrefix"]+"_longSV_summary_gene_info_single.txt",
		longSV_H3K27ac_info="results/longSV_analysis/"+config["outPrefix"]+"_longSV_summary_H3K27ac_single.bedpe",
		gene_list=lambda wildcards: config["gene_file"],
		chrEnd=lambda wildcards: config["chrEnd_rds"],
		H3K27ac_peaks=lambda wildcards: config["H3K27ac_bed"],
	output:
		"results/longSV_analysis/all_longSV_flanking_events_with_enhancer_single_summary_full_version.txt",
	params:
		script="script/create_final_annotation_table_for_longSVevents.R",
	shell:
		"Rscript {params.script} {input.longSV_summary} {input.longSV_H3K27ac_info} {input.gene_list} {input.chrEnd} {input.H3K27ac_peaks}"

#Create comut file for Short SV H3k27ac enhancer events (from step 7-1)

rule create_comut_file_for_shortSV_enh:
	input:
		"results/shortSV_enh_analysis/"+config["outPrefix"]+"_shortSV_enhUpMatrix_H3K27ac_single.txt",
		"results/shortSV_enh_analysis/"+config["outPrefix"]+"_shortSV_enhDownMatrix_H3K27ac_single.txt",
		gene_list=lambda wildcards: config["gene_file"],
	output:
		# "results/comut_files/comut_shortSV_enhUpMatrix_H3K27ac_single.csv",
		# "results/comut_files/comut_shortSV_enhDownMatrix_H3K27ac_single.csv",
		"results/shortSV_enh_analysis/shortSV_enh_analysis_complete.txt",
	params:
		script="script/create_comut_from_SVmatrixByGenes_with_H3K27ac.R",
		prefix=lambda wildcards: config["outPrefix"]+"_shortSV",
		out_dir="results/comut_files/"
	shell:
		"Rscript {params.script} {params.prefix} {input.gene_list} {params.out_dir}"


# Merge all comut files into one
rule combine_all:
	input:
		"results/comut_files/comut_intersect_sv_and_gene_body.csv",
		"results/comut_files/comut_intersect_sv_and_flanking_up.csv",
		"results/comut_files/comut_intersect_sv_and_flanking_dn.csv",
		"results/longSV_analysis/all_longSV_flanking_events_with_enhancer_single_summary_full_version.txt",
		"results/shortSV_enh_analysis/shortSV_enh_analysis_complete.txt",
	output:
		"results/comut_files/all_comut_SV_combined_single.csv",
	params:
		script="script/combine_all_SV_comut.R",
	shell:
		"Rscript {params.script}"


