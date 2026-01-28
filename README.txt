Project Genome-wide transcript end mapping reveals novel processing events in plant organelles

Investigators: Leni Campbell-Clause, Ian Small, Catherine Colas des Francs-Small

This project contains

Julia scripts circus2.jl and circus3.jl for the mapping of specialised short-read RNAseq used to capture primary and processed transcript ends in Arabidopsis thaliana chloroplasts and mitochondria

The directory structure of this project is as follows:

>chloroplast
#julia scripts, notebooks, and output files are stored here
	>circus2
		>circularisation_checks
		#mappings output from circus - calculation of circularisation rates
		>circus2_deseq2_comp
		#jupyter notebook and files for comparison - beta distribution and deseq2
		>circus2.jl
		#circus2 mapping script
		>intron_splicing
		#jupyter notebook and known splice site comparison
	>circus3
		>AP000423.gff
		>beta_distribution.ipynb
		#jupyter notebook for beta distribution
		>circus3.jl
		#circus3 mapping script
		>mappings3
		#mappings3 output from circus
		>primary
		#jupyter notebook and files for comparisons
		>processed
		#jupyter notebook and files for comparisons
>mitochondria
#raw rna sequence data is stored here
	>processed
		>beta_distribution_output
		#output files from beta distribution
		>ends_BK010421.1.csv
		#known transcript ends with positions
		>ends_with_nearest_sig.csv
		#output file from overlaps_processed_RF_ipynb
		>joined_3prime
		#output file from overlaps_processed_RF_ipynb
		>joined_5prime
		#output file from overlaps_processed_RF_ipynb
		>overlaps_processed_RF_ipynb
		#jupyter notebook for comparison - beta distribution and known ends
		>TABLE5.xlsx
		#comparison table - beta distribution and known ends
	>RF
		>beta_distribution_RF.ipynb
		#jupyter notebook for beta distribution of RF libraries		>BK010421.gff3
		#annotation file for Arabidopsis thaliana mitochondrial genome		>circus3.jl
		#circus3 mapping script		>comparison_filemaker.ipynb
		#conducts log2 comparison analysis		>counts3
		#circus3.jl output files		>RF_plotting.ipynb
		#plotting of full genome and individual comparisons		>samples.tsv
		#library and sample names		>sig_comparisons
		#output from comparison_filemaker.ipynb containing log2 comparisons
