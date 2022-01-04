#!/bin/bash
set -euo pipefail

#wget -O data/problematic_sites_sarsCov2.vcf https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf

#cat ../genomes/*.fasta seqs/root_seqs.fasta seqs/open_context_genomes.fasta > seqs/all_consensus.fasta

mkdir -p augur
# mask problematic sites and start/end from defaults in ncov
augur mask --mask-sites $(awk -F $'\t' '$1!~/^#/ {printf $2; printf " "}' < data/problematic_sites_sarsCov2.vcf) --sequences seqs/all_consensus.fasta --mask-from-beginning 100 --mask-from-end 50 --output augur/masked_seqs.fasta

# align with mafft
augur align --sequences augur/masked_seqs.fasta --nthreads 8 --method mafft --reference-sequence data/MN908947.3.fasta --fill-gaps --output augur/aligned.fasta 

# build tree with iqtree matching ncov settings
augur tree --alignment augur/aligned.fasta --tree-builder-args '-ninit 10 -n 4' --output augur/tree.nwk --nthreads 8

## refine tree matching ncov settings 
augur refine --tree augur/tree.nwk --output-tree augur/refined_tree.nwk --output-node-data augur/branch_lengths.json --root Wuhan/Hu-1/2019 Wuhan/WH01/2019 
#
# infer ancestral states and mutations
augur ancestral --tree augur/refined_tree.nwk --alignment augur/aligned.fasta --output-node-data augur/nt_muts.json --inference joint --infer-ambiguous

# translate nt to aa
augur translate --tree augur/refined_tree.nwk --ancestral-sequences augur/nt_muts.json --reference-sequence data/MN908947.3.gbk --output-node-data augur/aa_muts.json

# add clade info
augur clades --tree augur/refined_tree.nwk --mutations augur/nt_muts.json augur/aa_muts.json --clades data/clades.tsv --output-node-data augur/clades.json

# pangolin assigned
pangolin --all_versions > pango_version
pangolin --outfile pangolin_assignments.tsv seqs/all_consensus.fasta
