library(tidyverse)
library(ape)
library(phytools)
# Pull out principal components for each protein within OG0003341,
# and orthogroup with is comprised almost entirely of ARPC5 protein 
# sequences.

# This will ultimately be done in two ways:
# 1) using the PCs calculated on shapemers across all proteins ("the 
# global set")
# 2) or, using PCs calculated from shapemers of just the proteins in 
# this orthogroup. 

# Get the uniprot accessions for these proteins
prots <- read_tsv('~/environment/analyses/orthofinder/TCS/OrthogroupAnnotations/GenericOrthogroups/og_uniprot_accessions/OG0003341_ProteinIDs.tsv')
unipAccs <- prots$Entry[!is.na(prots$Entry)]

# And remove any species not in this set
prots <- prots[which(prots$Entry %in% unipAccs),]

# Great, now read in the shapemer PCs (this is huge to begin with, but we'll clean up after).
# Load the python packages for use in R and to read in these data
#Initiate reticulate
library(reticulate)
np <- import("numpy")
pd <- import("pandas")

#Load TFIDF
proteins = np$load("../../data/shapemers/eukprot/eukprot_comparative_set_shapemer_tfidf_top5000_components.npy")
#Load protein names (connect with rownames of TFIDF)
n = np$load("../../data/shapemers/eukprot/eukprot_comparative_set_shapemer_tfidf_top5000_components_proteins.npy")

# identify which of these proteins are in our ARPC5 orthogroup/gene family
arpc5.index <- which(str_detect(n, paste(unipAccs, collapse = "|")))

# And pull these out
n.arpc5 <- n[arpc5.index]
proteins.arpc5 <- proteins[arpc5.index,]

# Clean up
rm(proteins)

#Add names to TFIDF matrix
rownames(proteins.arpc5) = unlist(lapply(strsplit(n.arpc5, "-"), function(v){v[2]}))

# How many proteins per species remain?
arpc5.prots <- prots[which(prots$Entry %in% rownames(proteins.arpc5)),]
spp.counts <- summary(as.factor(arpc5.prots$Taxon))
spp.counts

# Bummer. Still lots of duplicates. Will need to think about this. 
# Why don't we conduct a phylogenetic PCA on all the structural predictions, 
# using this to test, for each protein, the distance from the centroid of PC
# space. Then, for each species for which multiple protein structures remain, 
# we retain the one closest to the centroid. Conservative, as we'll lose the 
# unusual protein structures, even for proteins with presumably the same 
# function, but it's better to be conservative than to include proteins of 
# divergent function for this analysis. 

# Conduct a phylogenetic PCA using these shapemers. 
# Read in the species tree 
# Read in the gene tree
gtree <- read.tree('../pargenes/EukProt_TCS_ExtremeCoreOGs/mlsearch_run/results/OG0003341-clipkit_fa/OG0003341-clipkit_fa.raxml.bestTree')

# Reduce down to the tips for which we have structural predictions
gtree <- drop.tip(gtree, gtree$tip.label[-which(gtree$tip.label %in% arpc5.prots$EukProt_ID)])

# And rename tips with their accessions
gtree$tip.label <- arpc5.prots$Entry[match(gtree$tip.label, arpc5.prots$EukProt_ID)]

# Pull out the terminal branch lengths - perhaps those proteins that have been 
# already identified as ARPC5 homologs are characterized by shallower terminal
# branch lengths, meaning we can just retain the protein per species with the 
# shortest terminal branch-lengths
# Start by pulling out all terminal edge lengths and naming them according to 
# the respective protein sequence. 
arpc5.edgelens <- 
setNames(gtree$edge.length[sapply(1:length(gtree$tip.label),
         function(x,y) which (y==x),y=gtree$edge[,2])],gtree$tip.label)
         
# What is the length for annotated ARPC5 homologs? For the rest?
summary(arpc5.edgelens[which(names(arpc5.edgelens) %in% arpc5.prots$Entry[which(arpc5.prots$`Protein names` == "Actin-related protein 2/3 complex subunit 5")])])
summary(arpc5.edgelens[which(names(arpc5.edgelens) %in% arpc5.prots$Entry[-which(arpc5.prots$`Protein names` == "Actin-related protein 2/3 complex subunit 5")])])


# get the acessions of proteins that have previously been identified as ARPC5 homologs. 
arpc5.accs <- arpc5.prots$Entry[which(arpc5.prots$`Protein names` == "Actin-related protein 2/3 complex subunit 5")]

# calculate pairwise distance using the raw shapemers
arpc5.dist <- dist(proteins.arpc5)

# Get distance of each protein from the centroid
centroid.dists <- usedist::dist_to_centroids(d = arpc5.dist, g = rep(as.factor('A'), nrow(proteins.arpc5)))

# And summarize
summary(centroid.dists$CentroidDistance[which(centroid.dists$Item %in% arpc5.accs)])
t.test(centroid.dists$CentroidDistance[which(centroid.dists$Item %in% arpc5.accs)], 
       centroid.dists$CentroidDistance[-which(centroid.dists$Item %in% arpc5.accs)],
       alternative = 'less')

# Alright. Not a wildy robust test, but it does seem that the homologs externally 
# annotated as ARPC5 are closer to the centroid of the distance matrix calculated 
# from shapemers. So, for all species that have multiple orthologs, retain only 
# the one that is closest to the centroid of that distance matrix. 

# pull out the ones we don't have to worry about first (single og per spp). 
arpc5.single.prots <- arpc5.prots[which(arpc5.prots$Taxon %in% names(spp.counts[which(spp.counts == 1)])),]

# Now go through and pull out one protein per species with multiple orthologs, keeping the one closest to the centroid
arpc5.multi.prots <- arpc5.prots[-which(arpc5.prots$Taxon %in% names(spp.counts[which(spp.counts == 1)])),]
for(i in 1:length(unique(arpc5.multi.prots$Taxon))){
    spp <- unique(arpc5.multi.prots$Taxon)[i]
    multi <- arpc5.multi.prots[which(arpc5.multi.prots$Taxon == spp),]
    dists <- centroid.dists[which(centroid.dists$Item %in% 
                            arpc5.multi.prots$Entry[which(arpc5.multi.prots$Taxon == spp)]),]
    mindist <- which(dists$CentroidDistance == min(dists$CentroidDistance))
    keep <- multi[mindist,]
    arpc5.single.prots <- rbind(arpc5.single.prots, keep)
}

# Fantastic. Now we can conduct the phylogenetic PCA to summarize the 
# protein structural predictions!
# First read in the species tree and rename the tips to just Genus_species and 
# prune down. We'll then rename both the proteins and tips as "Genus_species_UniProtAccession"
spp.tree <- read.tree('~/environment/projects/protein-structural-evolution/speciesrax/TCS-SpeciesRax-Congruified-WP2011.nwk')
spp.tree <- drop.tip(spp.tree, spp.tree$tip.label[-which(spp.tree$tip.label %in% arpc5.single.prots$Taxon)])

# Pull out the shapemers for the single protein per species that we retained
sc.arpc5 <- proteins.arpc5[which(rownames(proteins.arpc5) %in% arpc5.single.prots$Entry),]

# Rename tips
arpc5.single.prots <- arpc5.single.prots[match(spp.tree$tip.label, arpc5.single.prots$Taxon),] 
spp.tree$tip.label <- paste(arpc5.single.prots$Taxon, arpc5.single.prots$Entry, sep = "_")
# and proteins
arpc5.single.prots <- arpc5.single.prots[match(rownames(sc.arpc5), arpc5.single.prots$Entry),]
rownames(sc.arpc5) <- paste(arpc5.single.prots$Taxon, arpc5.single.prots$Entry, sep = "_")

# And sort the shapemers to be in order of the tips in the tree.
sc.arpc5 <- sc.arpc5[match(spp.tree$tip.label, rownames(sc.arpc5)),]

# And run the pPCA!
arpc5.ppca <- phyl.pca(spp.tree, sc.arpc5, method="BM", mode="corr")

# Write the shapemers and results of pPCA out to file
write.table(sc.arpc5, file = 'SingleCopy-ARPC5-Shapemers.tsv', 
            sep = '\t', quote = F, row.names = T, col.names = T)
saveRDS(arpc5.ppca, file = "SingleCopy-ARPC5-Shapemer-phyloPCA.RDS") 

# Also just run a normal PCA to explore what's driving the pattern of PCs 
# primarily explaining differences between closely related pairs of species. 
