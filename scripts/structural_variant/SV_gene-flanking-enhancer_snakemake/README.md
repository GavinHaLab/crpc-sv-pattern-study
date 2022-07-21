# Overview of Analysis Workflow

* SV calls assignment of each event to gene/flanking/enhancer by prioritizing SV events.

<img src="https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/SV-gene-flanking-pipeline.png" title="Analysis Workflow" height="300" width="700">

## Definition of SV events
### For short SV events (SPAN < 10mb)

<img src="https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/short_SV_events_definition.png" title="Analysis Workflow" height="400" width="600">

### For translocation events including intra-chromosomal long SV events (SPAN >= 10mb), flankingGene Left

<img src="https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/long_SV_events_definition_left.png" title="Analysis Workflow" height="400" width="800">

### For translocation events including intra-chromosomal long SV events (SPAN >= 10mb), flankingGene Right

<img src="https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/long_SV_events_definition_right.png" title="Analysis Workflow" height="400" width="800">

For translocation events including intra-chromosomal long SV, the direction of the SV should go towards the gene body region.
Long flanking with enhancer (gain of enhancer) is interesting (Enhancer-Associated Rearrangements); hypothesized that enhancer hijacking translocation events that bring potent enhancers to driver gene.

## Prioritize SV events 

* Gene transecting events are prioritized over flanking events for SV calls assignments.

<img src="https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/SV_gene_flanking_double_count_example.png" title="Analysis Workflow" height="400" width="800">

## Snakemake

* This workflow takes SV bepe (bedpe_prefix), sample list (sample_file) and gene list (gene_file). Please update this to your config.yaml. Current gene annotation information is in hg38 (data folder).

* SV bedpe format can be found from https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/CRPC10X42-PoNToolFilter_manual_edit_WCDT101_svaba-titan.min1kb.bedpe

* sample file format can be found from https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/sampleList_allCases.txt

* gene file format can be found from https://github.com/GavinHaLab/crpc-sv-pattern-study/blob/main/metadata/combined_gene_list_159.txt


