library(ape)
library(geiger)
library(phytools)
library(ggtree)
library(deeptime)

setwd('~/Documents/ArcadiaScience/Projects/Protein-SeqStruct-Evolution/TCS/')
tree <- read.tree('FastME-TCS-Rooted-SpeciesTree.nwk')
tree$tip.label <- sub(".*?_", "", tree$tip.label)
tree <- ladderize(tree)

# Figure out which species are even in UniProt and reduce down to this subset. 
up.tax <- read.table('./UniProt-Queryable-Taxa.txt', header = T, sep ='\t')
tree <- 
  drop.tip(tree, 
           tree$tip.label[-which(tree$tip.label %in% 
                                   up.tax$Name_to_Use)])
# spp.ids <- read.table('../EukProt/EukProt-TaxonSets/the-comparative-set-spp-names.tsv-spp-ids.txt', sep = ' ', header = F)
# spp.ids$V2 <- sub(".*?_", "", spp.ids$V2)
#spp.ids <-spp.ids[as.numeric(tree$tip.label),] 

# for(i in 1:nrow(spp.ids)){
#   tree$tip.label[which(tree$tip.label == i-1)] <- spp.ids$V2[i]
# }

plot.phylo(tree, cex = 0.55)

timetrees.table <- 
  read.table('~/Documents/ArcadiaScience/Projects/Protein-SeqStruct-Evolution/papers/Wegener-Parfrey-et-al-2011_EukaryoteTImeCal-Phylo_PNAS-TableS1-Trees.tsv', 
             sep = '\t', header = T)
timetrees <- list()
for(i in 1:16){
  timetrees[[i]] <- ladderize(
    castor::read_tree(string = timetrees.table$Newick.tree.string.with.divergence.times[i])
)
}
names(timetrees) <- timetrees.table$Name[1:16]

# Now, we're going to congruify our tree using all 16 time-calibrated phylogenies
# from Parfrey et al., 2011. 
for(tr in 1:16){
  if(tr == 1){
    congruified.trees <- list()
    genera.tt <- list()
    tt.df <- list()
    timetrees.clean <- list()
    
    genera.ml <- unlist(lapply(strsplit(tree$tip.label, '_'), `[[`, 1))
  }
  # We need to do some cleanup - keep only tips for genus that is found in the focal tree
  genera.tt[[tr]] <- unlist(lapply(strsplit(timetrees[[tr]]$tip.label, '_'), `[[`, 1))
  tt.df[[tr]] <- 
    data.frame(Genus = genera.tt[[tr]], 
               Species = timetrees[[tr]]$tip.label,
               Keep = 0,
               InML = 0)
  for(i in 1:nrow(tt.df[[tr]])){
    if(tt.df[[tr]]$Genus[i] %in% genera.ml){tt.df[[tr]]$Keep[i] <- 1}
    if(tt.df[[tr]]$Species[i] %in% tree$tip.label){tt.df[[tr]]$InML[i] <- 1}
  }
  
  # Drop the odd tip that isn't in our ML tree
  tt.df[[tr]] <- tt.df[[tr]][-which(tt.df[[tr]]$Keep == 0),]
  # Keep the first of replicate samples in each genus
  tt.df[[tr]] <- tt.df[[tr]][!duplicated(tt.df[[tr]]$Genus),]
  
  # And then for the remaining species names that don't match, replace with the 
  # species name for that genus in the ML tree
  tt.df[[tr]]$SpeciesNew <- tt.df$Species
  for(i in 1:nrow(tt.df[[tr]])){
    tt.df[[tr]]$SpeciesNew[i] <- 
      tree$tip.label[which(genera.ml == tt.df[[tr]]$Genus[i])][1]
  }
  
  # Clean up the old tree
  timetrees.clean[[tr]] <- 
    drop.tip(timetrees[[tr]], 
             timetrees[[tr]]$tip.label[-which(timetrees[[tr]]$tip.label %in%
                                                tt.df[[tr]]$Species)])
  
  # Now rename the tip labels with these values to make usable with congruify
  for(i in 1:nrow(tt.df[[tr]])){
    timetrees.clean[[tr]]$tip.label[i] <- tt.df[[tr]]$SpeciesNew[i]
  }
  
  plot(ladderize(timetrees.clean[[tr]]))
  timetrees.clean[[tr]] <- force.ultrametric(timetrees.clean[[tr]], method = 'extend')
  
  congruified.trees[[tr]] <- 
    congruify.phylo(reference = timetrees.clean[[tr]], target = tree, 
                    ncores = 6, scale = 'treePL')
  
  congruified.trees[[tr]]$phy <- ladderize(congruified.trees[[tr]]$phy)
}

# Get a sense of how similar each time-calibrated tree is from the tree inferred
# from the EukProt dataset
for(i in 1:16){
  obs.tips <- tree$tip.label
  tt.tips <- timetrees.clean[[i]]$tip.label
  
  tr1 <- 
    drop.tip(tree, 
             obs.tips[-which(obs.tips %in% tt.tips)])
  tr2 <- 
    drop.tip(timetrees.clean[[i]], 
             tt.tips[-which(tt.tips %in% tree$tip.label)])
  if(i == 1){
    rf.dists <- castor::tree_distance(tr1, tr2)
  }else{
    rf.dists[i] <- castor::tree_distance(tr1, tr2)
  }
}

plot(rf.dists)
# tree 4 is closest to the tree we estimated using the Eukprot, by a decent
# amount - move forward with this one. 

#congruified.arch$phy <- root(congruified.arch$phy, interactive = T)
plot(congruified.trees[[4]]$phy, cex = 0.5)

cphy <- cophylo(congruified.trees[[4]]$phy, timetrees[[4]])
plot.cophylo(cophylo(tree, timetrees[[4]]), fsize = 0.5)






# tcs.timetree <- 
#   read.tree('../EukProt/EukProt-TaxonSets/the-comparative-set-time-tree.nwk')
# 
# # The group containing giardia seems problematic for dating... remove
# # for now. 
# # also, 'Naegleria gruberi' seems problematic
# # And a couple problem taxa
# prob.tax <- 
#   c('Carpediemonas_membranifera', 'Giardia_intestinalis', 
#     'Paratrimastix', 'Naegleria_gruberi', 
#     'Mastigamoeba_balamuthi', 'Entamoeba_histolytica',
#     'Acanthamoeba_castellanii', 'Amoeba_proteus', 'Arcella_intermedia')
# # and the amoeboids
# 
# 
# tcs.timetree <-
#   drop.tip(tcs.timetree, prob.tax)

# # There is a differences in species synonym that need to be resolved. 
# # Fix this first.
# tcs.timetree$tip.label[which(tcs.timetree$tip.label == 'Nonionellina_labradorica')] <-
#   'Nonionella_labradorica'
#   
# # We need to do some cleanup - keep only tips for genus that is found in the focal tree
# genera.tt <- unlist(lapply(strsplit(tcs.timetree$tip.label, '_'), `[[`, 1))
# genera.ml <- unlist(lapply(strsplit(tree$tip.label, '_'), `[[`, 1))
# tt.df <- 
#   data.frame(Genus = genera.tt, 
#              Species = tcs.timetree$tip.label,
#              Keep = 0,
#              InML = 0)
# for(i in 1:nrow(tt.df)){
#   if(tt.df$Genus[i] %in% genera.ml){tt.df$Keep[i] <- 1}
#   if(tt.df$Species[i] %in% tree$tip.label){tt.df$InML[i] <- 1}
# }
# 
# # Drop the odd tip that isn't in our ML tree
# tt.df <- tt.df[-which(tt.df$Keep == 0),]
# # Keep the first of replicate samples in each genus
# tt.df <- tt.df[!duplicated(tt.df$Genus),]
# 
# # And then for the remaining species names that don't match, replace with the 
# # species name for that genus in the ML tree
# tt.df$SpeciesNew <- tt.df$Species
# for(i in 1:nrow(tt.df)){
#   tt.df$SpeciesNew[i] <- 
#     tree$tip.label[which(genera.ml == tt.df$Genus[i])][1]
# }
# 
# # Clean up the old tree
# tcs.timetree <- 
#   drop.tip(tcs.timetree, 
#            tcs.timetree$tip.label[-which(tcs.timetree$tip.label %in%
#                                            tt.df$Species)])
# 
# # Now rename the tip labels with these values to make usable with congruify
# for(i in 1:nrow(tt.df)){
#   tcs.timetree$tip.label[i] <- tt.df$SpeciesNew[i]
# }
# 
# plot(ladderize(tcs.timetree))
# tcs.timetree <- force.ultrametric(tcs.timetree, method = 'extend')

# tree <- root(tree, interactive = T)
# Now use congruify
# Think I need to change "smooth_scion" in the internal functions when using "scale".
# fixInNamespace(ns  = 'geiger', congruify.phylo)
congruified.tcs <- 
  congruify.phylo(reference = tcs.timetree, target = tree, 
                  ncores = 6, scale = 'treePL')

congruified.tcs$phy <- ladderize(congruified.tcs$phy)
#congruified.arch$phy <- root(congruified.arch$phy, interactive = T)
plot(congruified.tcs$phy, cex = 0.5)
cphy <- cophylo(tree, tcs.timetree)
plot(cphy, fsize = c(0.01, 0.5))
metadat <- 
  read.table('../EukProt/EukProt-Taxonomic-Metadata.tsv', header = T, sep = '\t')
tax.counts <- summary(as.factor(metadat$Supergroup_UniEuk))
core.groups <- names(tax.counts[which(tax.counts >5)])

sup.grp.trees <- list()
sup.grp.tree.plts <- list()
tips <- congruified.tcs$phy$tip.label
for(i in 1:length(core.groups)){
  taxon <- core.groups[i]
  spps <- metadat$Name_to_Use[which(metadat$Supergroup_UniEuk == taxon)]
  sup.grp.trees[[i]] <- 
    drop.tip(congruified.tcs$phy, 
             tip = tips[-which(tips %in% spps)])
  
  p <- 
    ggtree(sup.grp.trees[[i]]) +
    theme_tree2() + 
    ggtitle(paste0(taxon))
  p <- revts(p)
  sup.grp.tree.plts[[i]] <- 
    gggeo_scale(p, neg = T)
}

p <- p %<+% tip_metadata + geom_tippoint(aes(color=age_group), size=3)

tcs.p <- 
  ggtree(tree) + 
  theme_tree2()
tcs.p <- revts(tcs.p)

tip.metadat <- metadat[,c(2,5)]
colnames(tip.metadat)[1] <- 'taxa'
tcs.p %<+% 
  tip.metadat + 
  geom_tippoint(aes(color = Supergroup_UniEuk), size = 3) + 
  geom_tiplab(size = 2.5)
names(sup.grp.trees) <- core.groups
plot(sup.grp.trees$Alveolata)



