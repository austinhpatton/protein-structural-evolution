# README 
# This directory contains orthologs (both orthogroups [OGs] and hierarchical orthogroups [HOGs]) 
# as inferred using "The Comparative Set" of 196 species of Eukaryotes from the EukProt Dataset,
# and orthofinder. 

# To cluster protein sequences into initial orthogroups, MCL clustering was applied to 
# phylogenetically corrected similarity scores from All-v-All diamond comparisons. 
# For each, MAFFT L-INS-I was applied to infer multiple sequence alignments.
# RAxML-NG was used to infer gene-family/orthogroup trees (coordinated/parallelized using ParGenes).
# These gene family trees will then be used within the SpeciesRax pipeline to infer a species tree, 
# under a model of duplication, transfer and loss, and subsequently use this species tree for
# gene-tree species-tree reconciliation. 

# Directories contained herein (to be updated) includes:
# 1) ./OGs/ : "Raw" orthogroups inferred from MCL clustering of similarity scores
# 2) ./OGs/Core/ : "Core orthogroups" that include at least ten species, with a mean 
	                        per-species copy number of 3. These are used for species tree inference.
# 3) ./OGs/Core/MSAs : "Core orthogroup" multiple sequence alignments, trimmed using ClipKit.
# 4) ./OGs/Core/Trees : "Core orthogroup" gene-family trees
# 5) ./OGs/Remainder/ : Every other orthogroup. These are particularly large, or small. 
# 6) ./OGs/Remainder/MSAs : "Remainder" multiple sequence alignments
# 7) ./OGs/Remainder/Trees : "Remainder" gene-family trees
