Overview

This pipeline is designed for inferring species trees and network evolution history from gene trees, as well as performing ancestral state reconstruction.

Scripts and Tools:

ABBA_BABA_analysis_NEW.py
Performs the ABBA-BABA test to detect gene flow and introgression between species.

ancestral_reconstruction.r
Uses the phytools package to perform ancestral state reconstruction, inferring ancestral traits for each internal node of the phylogenetic tree.

ancestral_state.R
This R script performs ancestral state reconstruction and evaluates Blomberg’s K and P values to assess the phylogenetic signal in traits. Additionally, it uses RRphylo to detect evolutionary shifts in trait evolution along the phylogeny.

mcmc_compare.R
Compares MCMC results under different root age constraints. 

plot_ABBA.R
Visualizes ABBA-BABA results using heatmaps to show the distribution of gene flow and introgression patterns across species/individuals.

process_species.ABBA_avg.py
Processes species data and computes the average ABBA-BABA values for species groups, summarizing gene flow and introgression patterns.

process_species.ABBA_avg_multiInd.py
Similar to the above script but processes multiple individuals per species, suitable for population-level data analysis.

retrieve_triple_from_tree.MultiInd.py
Extracts specific relationships (triplets) from phylogenetic trees with multiple individuals per species.

run_tree.sh
A shell script that filters low bootstrap branches using newick-utils, infers the species tree with ASTRAL, and performs network evolution analysis using PhyloNet to identify hybridization and gene flow events.
