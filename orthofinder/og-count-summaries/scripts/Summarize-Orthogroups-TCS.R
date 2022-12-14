# Summaries of orthogroups and comparative genomic statistics for The Comparative SeT
# as inferred using OrthoFinder
library(tidyverse)
library(reshape)
library(ape)
library(vroom)
library(data.table)
#library(M3C)
library(MetBrewer)
library(ggExtra)
library(cowplot)
library(foreach)
library(umap)
library(ClusterR)
library(reticulate)

# If on laptop:
basedir <- '~/Documents/ArcadiaScience/Projects/Protein-SeqStruct-Evolution'

# If on home desktop:
basedir <- '~/ArcadiaScience/Protein-SeqStruct-Evolution/orthofinder/TCS/'
gitdir <- '~/ArcadiaScience/github/protein-structural-evolution/orthofinder/og-count-summaries/'

setwd(basedir)

# Read in metadata for samples in the comparative set
#metadat <- 
#  read.table('~/Arcadia Science Dropbox/Austin Patton/Protein-SeqStruct-Evolution/EukProt-V3-Data/EukProt-Taxonomic-Metadata.tsv', 
#             header = T, sep = '\t')

# Read in the gene counts for each orthogroup, for each sample/species in TCS. 
counts.mat <- vroom('Orthogroups.GeneCount.tsv', delim = '\t', col_names = T)
og.no <- vroom('./Results_Aug10_1/Phylogenetic_Hierarchical_Orthogroups/N0.tsv')
n0.ogs <- og.no[which(og.no$`Gene Tree Parent Clade` == 'n0'),]

# # Summarize gene couts of each N0.HOG for each species. 
# for(hog in 1:length(og.no$HOG)){
#   og.no
# }

# A function to summarize gene counts per orthogroup (e.g. mean, median, SD, SE) for a set of species. 
summarize_counts <- 
  function(counts.mat, species = NULL){
    if(is.null(species)){
      species <- colnames(counts.mat[-ncol(counts.mat)])
    }
    og.counts <- 
      counts.mat[,which(colnames(counts.mat) %in% species)]
    
    og.summs <- 
      data.frame(Orthogroup = counts.mat$Orthogroup,
                 Mean_CountPerSpp = NA,
                 Median_CountPerSpp = NA,
                 Variance_CountPerSpp = NA,
                 SD_CountPerSpp = NA,
                 SE_CountPerSpp = NA)
    
    st.dev <- 
      function(x){
        sqrt(var(x))
      }
    st.err <- 
      function(x){
        st.dev(x)/sqrt(length(x))
      }
    
    print('Calculating mean OG count across species')
    og.summs$Mean_CountPerSpp <- 
      apply(og.counts[-1], 1, mean)
    
    print('Calculating median OG count across species')
    og.summs$Median_CountPerSpp <- 
      apply(og.counts[-1], 1, median)
    
    print('Calculating variance in OG count across species')
    og.summs$Variance_CountPerSpp <- 
      apply(og.counts[-1], 1, var)
    
    print('Calculating standard deviation of OG count across species')
    og.summs$SD_CountPerSpp <- 
      apply(og.counts[-1], 1, st.dev)

    print('Calculating standard error of OG count across species')
    og.summs$SE_CountPerSpp <- 
      apply(og.counts[-1], 1, st.err)
    
    return(og.summs)
  }

og.count.summs <- 
  summarize_counts(counts.mat)
og.count.summs.mlt <- reshape::melt(og.count.summs, id.vars = 'Orthogroup')
# plot them
mean.count.plt <- 
  ggplot(data = og.count.summs, aes(x = Mean_CountPerSpp)) + 
  geom_histogram(fill = 'grey', color = 'black', bins = 50) + 
  scale_x_log10() + 
  annotation_logticks(sides = 'b')+
  theme_bw(base_size = 14) + 
  xlab('Mean Per-Species OG Count')

median.count.plt <- 
  ggplot(data = og.count.summs, aes(x = Median_CountPerSpp)) + 
  geom_histogram(fill = 'grey', color = 'black', bins = 50) + 
  scale_x_log10() + 
  annotation_logticks(sides = 'b') +
  theme_bw(base_size = 14) + 
  xlab('Median Per-Species OG Count')

var.count.plt <- 
  ggplot(data = og.count.summs, aes(x = Variance_CountPerSpp)) + 
  geom_histogram(fill = 'grey', color = 'black', bins = 50) + 
  scale_x_log10() + 
  annotation_logticks(sides = 'b') +
  theme_bw(base_size = 14) + 
  xlab('Variance of Per-Species OG Count')

se.count.plt <- 
  ggplot(data = og.count.summs, aes(x = SE_CountPerSpp)) + 
  geom_histogram(fill = 'grey', color = 'black', bins = 50) + 
  scale_x_log10() + 
  annotation_logticks(sides = 'b') +
  theme_bw(base_size = 14) + 
  xlab('Std. Error of Per-Species OG Count')

sd.count.plt <- 
  ggplot(data = og.count.summs, aes(x = SD_CountPerSpp)) + 
  geom_histogram(fill = 'grey', color = 'black', bins = 50) + 
  scale_x_log10() + 
  annotation_logticks(sides = 'b') +
  theme_bw(base_size = 14) + 
  xlab('Std. Deviation of Per-Species OG Count')

count.summ.plt <- 
  plot_grid(mean.count.plt, median.count.plt,
            sd.count.plt, se.count.plt, nrow = 2, 
            labels = c('a.', 'b.', 'c.', 'd.'))

ggsave(count.summ.plt, file = 'TCS-OG-CountSummaries.pdf', 
       width = 10, height = 10)

# Lets summarize these data using UMAP
# Start by pulling out the subset of our data that includes only those species for which we may correspond proteins to UniProt/AlphaFold.
up.tax <- read.table('../../../github/protein-structural-evolution/UniProt-Queryable-Taxa.txt', header = T, sep ='\t')
tax <- colnames(counts.mat)[1:196]
spps <- sub(tax, pattern = ".*?_", replacement = "")
tax <- data.frame(CombinedID = tax, Species = spps)
up.tax <- tax[which(tax$Species %in% up.tax[,1]),]

# Read in taxonomic metadata. 
taxonomy <- 
  read.table('~/ArcadiaScience/Protein-SeqStruct-Evolution/EukProt-V3-Data/EukProt-Taxonomic-Metadata.tsv',
             header = T)
taxonomy <- taxonomy[which(taxonomy$Name_to_Use %in% up.tax$Species),]
taxonomy$CombinedID <- paste(taxonomy$EukProt_ID, taxonomy$Name_to_Use, sep = "_")

# Pull out the names of supergroups, removing any with only a single species. This will introduce issues downstream. 
sup.grps <- names(which(summary(factor(taxonomy$Supergroup_UniEuk)) >= 2))

# Only include the count data for these supergroups - otherwise things'll turn out wonky
keepers <- taxonomy$CombinedID[which(taxonomy$Supergroup_UniEuk %in% sup.grps)]

up.og.counts.mat <- counts.mat[,c(1, which(colnames(counts.mat) %in% keepers))]
rm(counts.mat)


# And pull out the set of species in each supergroup
tax.sets <- list()
for(grp in 1:length(sup.grps)){
  tax.sets[[grp]] <-
    taxonomy[which(taxonomy$Supergroup_UniEuk == sup.grps[grp]),]
}
names(tax.sets) <- sup.grps

voom <- 
  function(counts){
    # Transform each count as a sample-specific transformation
    for(samp in 1:ncol(counts)){
      samp.counts <- unlist(c(counts[,samp]))
      voom <- log2((samp.counts+0.5) / (sum(samp.counts)+1) * 1000000)
      counts[,samp] <- voom
    }
    return(counts)
  }

# Set up a function to pull out and format the gene count data for any taxonomic group:
format.counts <- 
  function(taxonomic.group, count.data, all.tax = T){
    print(paste0("Pulling out counts for ", taxonomic.group))
    if(taxonomic.group != "all"){
      # pull out the counts for the focal taxonomic group. 
      counts <- 
        count.data[,c(1,which(colnames(count.data) %in% 
                                tax.sets[[taxonomic.group]]$CombinedID))]
    }else{
      counts <- count.data
    }
    
    if(all.tax == T){
      # And remove any orthogroup that is present in fewer than two species. 
      # Only do this if formatting conts for the full set - we want the data
      # to be as consistent across taxonomic groups when making predictions.
      counts <- counts[which(rowSums(counts[,-1] > 0) > 1),]
    }

    # Prepare for the voom transformation
    voom.counts <- counts
    voom.counts[,-1] <- voom(voom.counts[,-1])
  
    return(CountData = voom.counts)
  }

run.umap <- function(count.set, umap.config){
  umap.res <-
    umap(count.set[,-1],
         config = umap.config,
         method = 'umap-learn')
  
  # Plot it.
  umap.df <- as.data.frame(umap.res$layout)
  colnames(umap.df) <- c('UMAP_Axis1', 'UMAP_Axis2')
  
  # Calculate variances of rows, and use to color when plotting
  umap.df$CountVariance <- apply(count.set[,-1], var, MARGIN = 1)
  umap.df$Log10CountVariance <- log10(umap.df$CountVariance)
  umap.df$MeanCount <- apply(count.set[,-1], mean, MARGIN = 1)
  umap.df$MinCount <- apply(count.set[,-1], min, MARGIN = 1)
  umap.df$MaxCount <- apply(count.set[,-1], max, MARGIN = 1)
  umap.df$MedianCount <- apply(count.set[,-1], median, MARGIN = 1)
  umap.df$SD <- apply(count.set[,-1], sd, MARGIN = 1)
  
  results <- list(umap = umap.res,
                  umap.df = umap.df)

  return(results)
}

pred.umap <- 
  function(umap.res, taxon.counts){
    # Project the focal taxon-set onto the full umap space
    umap.pred <- predict(object = umap.res, data = as.matrix(taxon.counts[,-1]))
    umap.pred <- data.frame(umap.pred)
    colnames(umap.pred) <- c('UMAP_Axis1', 'UMAP_Axis2')
    
    # Calculate variances of rows, and use to color when plotting
    umap.pred$CountVariance <- apply(taxon.counts[,-1], var, MARGIN = 1)
    umap.pred$Log10CountVariance <- log10(umap.pred$CountVariance+1)
    umap.pred$MeanCount <- apply(taxon.counts[,-1], mean, MARGIN = 1)
    umap.pred$MinCount <- apply(taxon.counts[,-1], min, MARGIN = 1)
    umap.pred$MaxCount <- apply(taxon.counts[,-1], max, MARGIN = 1)
    umap.pred$MedianCount <- apply(taxon.counts[,-1], median, MARGIN = 1)
    umap.pred$SD <- apply(taxon.counts[,-1], sd, MARGIN = 1)

    return(list(prediction = umap.pred))
  }

# pull out and transform the count data for the full set, and then for each taxon-set.
all.counts <- list()
all.counts[[1]] <- format.counts('all', count.data = up.og.counts.mat)
for(i in 1:length(tax.sets)){
  taxon <- names(tax.sets)[i]
  all.counts[[i+1]] <- format.counts(taxon, count.data = all.counts[[1]], all.tax = F)
}
names(all.counts) <- c('all', names(tax.sets))

# Set some configurations in preparation of running UMAP
umap.config <- umap.defaults
umap.config$n_neighbors <- 25
umap.config$knn_repeats <- 10
umap.config$metric <- 'euclidean'
umap.config$n_epochs <- 333
umap.config$min_dist <- 0.95
umap.config$spread <- 1
umap.config$verbose <- TRUE

# Run UMAP on all taxa in UniProt
all.umap.res <- run.umap(all.counts$all, umap.config)

# Plot it
all.umap.res$umap.plot <-
  ggplot(data = all.umap.res$umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2,
                             color = CountVariance)) +
  geom_point(size = 0.75, alpha = 0.5) +
  scale_color_met_c('Hiroshige', direction = -1) +
  theme_bw(base_size = 14) +
  xlab('UMAP Axis 1') +
  ylab('UMAP Axis 2') +
  labs(color = "Count Variance") +
  theme(legend.position = 'bottom')

# Now project each taxon-set into this space. 
tax.umap.projs <- list()
for(i in 1:length(sup.grps)){
  taxon <- names(tax.sets)[i]
  taxa <- tax.sets[[taxon]]$CombinedID
  
  print(paste0('Working on projecting ', taxon, '. Hang in there.'))
  tax.umap.projs[[i]] <- pred.umap(umap.res = all.umap.res$umap, taxon.counts = all.counts[[taxon]])
}
names(tax.umap.projs) <- names(tax.sets)


# Then build into a big plot. 
all.umap.proj.res <- all.umap.res$umap.df
all.umap.proj.res$TaxonSet <- 'All'
for(i in 1:length(sup.grps)){
  taxon <- sup.grps[i]
  taxon.project <- tax.umap.projs[[taxon]]$prediction
  taxon.project$TaxonSet <- taxon
  all.umap.proj.res <- rbind(all.umap.proj.res, taxon.project)
}
all.umap.proj.res$TaxonSet <- 
  factor(all.umap.proj.res$TaxonSet, 
         levels = unique(all.umap.proj.res$TaxonSet)[order(summary(factor(all.umap.proj.res$TaxonSet)), decreasing = T)], 
         ordered = T)

all.umap.proj.res$PlotOrder <- 0
to.jitter <- which(all.umap.proj.res$TaxonSet != "All")
all.umap.proj.res[to.jitter,"PlotOrder"] <- sample(1:11, size = length(to.jitter), replace = T)
all.umap.proj.res$PlotOrder <- factor(all.umap.proj.res$PlotOrder, levels = 0:11, ordered = T)
all.umap.proj.res <- all.umap.proj.res[order(all.umap.proj.res$PlotOrder),]

# Now build a single, unified plot
ggplot(data = all.umap.proj.res[-which(all.umap.proj.res$TaxonSet == 'All'),], 
       aes(x = UMAP_Axis1, y = UMAP_Axis2,
           color = TaxonSet)) +
  geom_point(alpha = 0.5, size = 0.75) +
  scale_color_manual(values = c('black', met.brewer('Lakota', n = 11, direction = 1))) +
  theme_bw(base_size = 14) +
  xlab('UMAP Axis 1') +
  ylab('UMAP Axis 2') +
    labs(color = "Taxon Set") +
  theme(legend.position = 'bottom')












umap.config$n_neighbors <- 30
umap.config$min_dist <- 0.65
opis.umap.res <- format.counts('Opisthokonta', up.og.counts.mat, umap.config)
opis.umap.res$umap.plot

umap.config$n_neighbors <- 25
umap.config$min_dist <- 0.85
chlor.umap.res <- format.counts('Chloroplastida', up.og.counts.mat, umap.config)
chlor.umap.res$umap.plot

# umap.config$n_neighbors <- 30
# umap.config$min_dist <- 0.5
umap.config$n_neighbors <- 20
umap.config$min_dist <- 0.95
alve.umap.res <- format.counts('Alveolata', up.og.counts.mat, umap.config)
alve.umap.res$umap.plot

# umap.config$min_dist <- 0.5
# all.umap.res <- format.counts('all', up.og.counts.mat, umap.config)



# Plot them all together. 
umap.plts <- 
  plot_grid(all.umap.res$umap.plot + ggtitle('All UniProt TCS'), 
            opis.umap.res$umap.plot + ggtitle('Opisthokonta'), 
            chlor.umap.res$umap.plot + ggtitle('Chloroplastida'), 
            alve.umap.res$umap.plot + ggtitle('Alveolata'),
            nrow = 2, ncol = 2, align = 'hv')

setwd(gitdir)
ggsave(umap.plts, filename = 'TCS-OG-Count-UMAP-Plots.pdf',
       height = 12, width = 16)
ggsave(umap.plts, filename = 'TCS-OG-Count-UMAP-Plots.png',
       height = 12, width = 16)



# Plot them all together. 
umap.pred.plts <- 
  plot_grid(all.umap.res$umap.plot + ggtitle('All UniProt TCS'), 
            opis.pred.umap$plot + ggtitle('Opisthokonta Projection'), 
            chlor.pred.umap$plot + ggtitle('Chloroplastida Projection'), 
            alve.pred.umap$plot + ggtitle('Alveolata Projection'),
            nrow = 2, ncol = 2, align = 'hv')

setwd(gitdir)
ggsave(umap.pred.plts, filename = 'TCS-OG-Count-UMAP-Projection-Plots.pdf',
       height = 12, width = 16)
ggsave(umap.pred.plts, filename = 'TCS-OG-Count-UMAP-Projection-Plots.png',
       height = 12, width = 16)

comb.preds <- rbind(all.umap.res$umap.df, opis.pred.umap$prediction,
                    chlor.pred.umap$prediction,alve.pred.umap$prediction)
comb.preds$set <- c(rep('All', nrow(all.umap.res$umap.df)),
                    rep('Oposthokonta', nrow(opis.pred.umap$prediction)),
                    rep('Chloroplastida', nrow(chlor.pred.umap$prediction)),
                    rep('Alveolata', nrow(alve.pred.umap$prediction)))
comb.preds$set <- 
  factor(comb.preds$set, levels = c('All', 'Oposthokonta', 
                                    'Chloroplastida', 'Alveolata'), 
         ordered = T)

comb.pred.plt <- 
  ggplot(data = comb.preds, aes(x = UMAP_Axis1, y = UMAP_Axis2, 
                               color = set)) + 
  geom_point(size = 0.75, alpha = 0.5) + 
  scale_color_met_d('Egypt', direction = -1) + 
  theme_bw(base_size = 14) + 
  guides(colour = guide_legend(override.aes = 
                                 list(size=3))) +
  xlab('UMAP Axis 1') +
  ylab('UMAP Axis 2') + 
  labs(color = "Tax-Set") + 
  theme(legend.position = 'bottom')

# Perform some clustering of the count data, independent of the umap embedding
# To do we we first need to obtain a similarity matrix to assemble into a graph. 
library(proxy)
chlor.counts.simil <- simil(chlor.umap.res$Log10_CountData[,-1], method = 'correlation')

# Convert this to an adjacency matrix
chlor.counts.adjac <-
  adjacency.fromSimilarity(chlor.counts.simil,
                           type = "signed",
                           power = if (type=="distance") 1 else 6)

# Now cluster, using the similarity
chlor.og.mcl <- mcl(x = chlor.counts.adjac, addLoops=TRUE, ESM = TRUE)

########################################################################
# Let's similarly re-work the full count matrix for all taxa to remove # 
# any orthologs that aren't present in at least two specie, and to     #
# log-transform those data.                                            #
########################################################################
# And remove any orthogroup that is only observed in a single species
single.spp <- c()
for(og in 1:nrow(up.og.counts.mat)){
  single.spp[og] <- as.numeric(length(which(up.og.counts.mat[og,-1] > 0)) == 1)
}
up.og.counts.mat <- up.og.counts.mat[-which(single.spp == 1),]

# Log10 transform to normalize *slightly*
up.og.log10.counts.mat <- up.og.counts.mat
opis.log10.counts[,-1] <- log10(opis.counts[,-1]+1)

opis.log10.counts <- opis.counts





up.og.log10.counts.mat <- up.og.counts.mat
up.og.log10.counts.mat[,-1] <- log10(up.og.log10.counts.mat[,-1]+1)
up.og.counts.umap <- 

#save(up.og.counts.umap, file = 'UniProt-TCS-OG-Count-UMAP-learn.Rds')
save(up.og.counts.umap, file = 'UniProt-TCS-OG-Log10-Count-UMAP-Euclid-25nn-0.05mindis-5knn-300epochs.Rds')
#load(file = 'UniProt-TCS-OG-Count-UMAP-learn.Rds')
umap.df <- as.data.frame(up.og.counts.umap$layout)
colnames(umap.df) <- c('UMAP_Axis1', 'UMAP_Axis2')

# Calculate variances of rows, and use to color when plotting
umap.df$CountVariance <- apply(up.og.counts.mat[,-1], var, MARGIN = 1)

ggplot(data = umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2, color = log(CountVariance+1))) + 
  geom_point(size = 0.5, alpha = 0.1) + scale_color_viridis_c() + 
  theme_bw(base_size = 14)

# Cool. Can we cluster them using Gaussian mixture models?
# Find the optimal number of clusters
opt_gmm <- 
  Optimal_Clusters_GMM(umap.df[,1:2], max_clusters = c(20:60), criterion = "BIC", 
                       dist_mode = "eucl_dist", seed_mode = "random_subset",
                       km_iter = 50, em_iter = 25, var_floor = 1e-10, 
                       plot_data = T)
plot(opt_gmm)
# And fit the full model for that optimal number
# We're actually going to just go with 25, since it's nearly as good as 50, 
# but obviously a bit less extreme
umap_gmm <- 
  GMM(umap.df[,1:2], gaussian_comps = 58, 
      dist_mode = 'maha_dist', seed_mode = "random_subset",
      km_iter = 100, em_iter = 100, var_floor = 1e-10)

# And predict so that we can plot
umap_gmm_preds <- 
  predict_GMM(umap.df[,1:2], umap_gmm$centroids, 
              umap_gmm$covariance_matrices, 
              umap_gmm$weights)

# And add the cluster predictions to the umap dataframe
umap.df$GMM_Cluster <- as.factor(umap_gmm_preds$cluster_labels)

# And plot again, coloring by cluster
og.count.umap.clusts <- 
  ggplot(data = umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2, 
                             color = GMM_Cluster)) + 
  geom_point(size = 0.5, alpha = 0.1) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = 'none')

# And zoom in on that central region
og.count.umap.clusts.core <- 
  ggplot(data = umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2, 
                             color = GMM_Cluster)) + 
  geom_point(size = 0.1, alpha = 0.1) + 
  coord_cartesian(xlim = c(-0.5, 6), ylim = c(14, 21)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = 'none')

########
# Same but for k-means clustering
opt_km <- 
  Optimal_Clusters_KMeans(umap.df[,1:2], max_clusters = 50, criterion = 'BIC', 
                          max_iters = 10, plot_clusters = T, num_init = 10)

plot(opt_km)
# And fit the full model for that optimal number
# We're actually going to just go with 25, since it's nearly as good as 50, 
# but obviously a bit less extreme
umap_km <- 
  GMM(umap.df[,1:2], gaussian_comps = 50, 
      dist_mode = 'maha_dist', seed_mode = "random_subset",
      km_iter = 100, em_iter = 100, var_floor = 1e-10)

# And predict so that we can plot
umap_km_preds <- 
  predict_GMM(umap.df[,1:2], umap_gmm$centroids, 
              umap_gmm$covariance_matrices, 
              umap_gmm$weights)

# And add the cluster predictions to the umap dataframe
umap.df$GMM_Cluster <- as.factor(umap_gmm_preds$cluster_labels)

# And plot again, coloring by cluster
og.count.umap.kclusts <- 
  ggplot(data = umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2, 
                             color = GMM_Cluster)) + 
  geom_point(size = 0.5, alpha = 0.1) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = 'none')

# And zoom in on that central region
og.count.umap.kclusts.core <- 
  ggplot(data = umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2, 
                             color = GMM_Cluster)) + 
  geom_point(size = 0.1, alpha = 0.1) + 
  coord_cartesian(xlim = c(-1, 5), ylim = c(3, 9)) + 
  theme_bw(base_size = 14) + 
  theme(legend.position = 'none')


# Clean up
tcs <- colnames(counts.mat)[-c(1, ncol(counts.mat))]
spp.id <- gsub("_.*", "", tcs)
spp <- sub(".*?_", "", tcs)
metadat <- metadat[which(metadat$EukProt_ID %in% spp.id),]
metadat$ID <- tcs

# Read in a table that depicts, for each count of orthogroups, the number of 
# species in each (i.e. # of OGs ~ # of Species in OG)
spp.per.og <- 
  read.table('./NumSpeciesPerOrthogroup-Hist.tsv', 
             sep = "\t", header = T)

# Plot an overall summary - the distribution of the number of orthologs for each
# count of species. 
spp.per.og.plt <- 
  ggplot(data = spp.per.og, 
         aes(x = Number.of.species.in.orthogroup, 
             y = Number.of.orthogroups, 
             color = Number.of.species.in.orthogroup)) + 
  geom_point(alpha = 0.8, size = 2) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) +
  xlab('# Species in OG') +
  ylab('# of Orthogroups') + 
  theme(legend.position = 'none')
spp.per.og.plt <- 
  ggMarginal(spp.per.og.plt, type = 'histogram', 
             margins = 'y', fill ='white')

# We actually want to do this for each taxonomic group though. 
# So, get the names of the "Supergroups"
tax.names <- unique(metadat$Supergroup_UniEuk)

# And then pull out the species lists for each
num.spp.in.og <- list()
prop.spp.in.og <- list()

# The below loop goes through, for each taxonomic group, and pulls out
# the number of species in each orthogroup. 
# So, the result for each taxonomic groups is a vector the same length as
# the total number of inferred orthogroups. 
# We then can take these vectors to summarize, for each taxonomic group, the 
# number of orthogroups with each count of species. 

# Do this in parallel since it's.... slow
# Initialize the parallel cluster
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK"
)
doParallel::registerDoParallel(cl = my.cluster)

# Great, now do the work
tax.counts <- list()
for(i in 1:length(tax.names)){
  print(paste0("Working on Species ", tax.names[i]))
  # Get the species names in that taxonomic group
  spp.set <- metadat$ID[which(metadat$Supergroup_UniEuk == tax.names[i])]
  
  # Pull the counts out
  tax.counts[[i]] <- counts.mat[which(colnames(counts.mat) %in% spp.set)]
  
  # and then go through the set of orthogroups to figure out 
  # the number of species in each orthogroup. 
  num.spp.in.og[[i]] <- 
    rep(NA, nrow(tax.counts[[i]]))
  
  breaks <- 
    hist(1:nrow(counts.mat), 
         breaks = 100, plot = F)$mids # Keep track of progress
  
  num.spp.in.og[[i]] <- 
    foreach(x = 1:nrow(counts.mat), 
          .combine = 'c') %dopar% {
            length(which(tax.counts[[i]][x,] > 0))
          }
  prop.spp.in.og[[i]] <- num.spp.in.og[[i]] / length(spp.set)
}

# and close the parallel cluster
parallel::stopCluster(cl = my.cluster)

names(num.spp.in.og) <- tax.names
names(prop.spp.in.og) <- tax.names

# Reformat to be more useful. 
all.num.spp.in.og <- 
  as.data.frame(cbind(counts.mat$Orthogroup, rlist::list.cbind(num.spp.in.og)))
all.prop.spp.in.og <- 
  as.data.frame(cbind(counts.mat$Orthogroup, rlist::list.cbind(prop.spp.in.og)))
all.num.spp.in.og[,-1] <- apply(all.num.spp.in.og[,-1], 2, as.numeric)
all.prop.spp.in.og[,-1] <- apply(all.prop.spp.in.og[,-1], 2, as.numeric)

# Initialize a dataframe that will, for each taxonomic group, provide the 
# "# of Species per Orthogroup ~ # of Orthogroups" frequency distribution
spp.per.og.allsets <- 
  as.data.frame(matrix(nrow = max(all.num.spp.in.og[-1])+1, 
                       ncol = length(tax.names)+1,
                       dimnames = list(NULL, c('NumSppPerOG', tax.names))))
spp.per.og.allsets$NumSppPerOG <- 0:max(all.num.spp.in.og[-1])

# And fill in
for(taxon in tax.names){
  freqs <- summary(as.factor(all.num.spp.in.og[taxon][,1]))
  spp.per.og.allsets[1:length(freqs),taxon] <- freqs
}

spp.per.og.allsets <- 
  reshape::melt(spp.per.og.allsets, 
                id.vars = "NumSppPerOG",
                variable_name = 'Supergroup')

# Plot for all supergroups the Orthogroup species density frequencies 
all.spp.og.plt <- 
  ggplot(data = spp.per.og.allsets, 
         aes(x = NumSppPerOG, y = value, 
             fill = Supergroup)) + 
  geom_point(alpha = 0.7, size = 2.5, pch = 21) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) +
  xlab('# Species in OG') +
  ylab('# of Orthogroups') + 
  theme(legend.position = 'none')
all.spp.og.plt

# There are a number of Supergroups with only one sample - this isn't necessarily super useful 
tax.sizes <- apply(all.num.spp.in.og[-1], 2, max)
core.set <- tax.names[which(tax.sizes > 5)]
spp.per.og.coresets <- 
  spp.per.og.allsets[which(spp.per.og.allsets$Supergroup %in% core.set),]

core.spp.og.plt <- 
  ggplot(data = spp.per.og.coresets, 
         aes(x = NumSppPerOG, y = value, 
             colour = Supergroup)) + 
  geom_point(alpha = 0.8, size = 2.5) +
  scale_color_met_d('Signac') +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) +
  xlab('# Species in OG') +
  ylab('# of Orthogroups') + 
  theme(legend.position = 'top')

core.spp.og.plt <- 
  ggMarginal(core.spp.og.plt, type = 'density',
             margins = 'y', groupColour = T)
core.spp.og.plt

# Now, let's also summarize the number of species within each supergroup
tax.sizes <- data.frame(tax.sizes)
tax.sizes$Supergroup <- rownames(tax.sizes)
tax.sizes$NumOGs <- NA
tax.sizes$MedianCount <- NA
for(i in 1:nrow(tax.sizes)){
  tax <- tax.sizes$Supergroup[i]
  tax.sizes$NumOGs[i] <- length(which(all.num.spp.in.og[tax] > 0 ))
}

tax.counts <- list()
tax.counts.vec <- list()
for(i in 1:length(tax.names)){
  print(i)
  # Get the species names in that taxonomic group
  spp.set <- metadat$ID[which(metadat$Supergroup_UniEuk == tax.names[i])]
  
  # Pull the counts out - we will use this to plot histograms/gensity plots of 
  # gene counts per Supergroup. 
  tax.counts[[i]] <- counts.mat[,which(colnames(counts.mat) %in% spp.set)]
  tax.counts.vec[[i]] <- Reduce(c, c(tax.counts[[i]]))
}

ggplot(data = tax.sizes, aes(y = NumOGs, x = tax.sizes)) + 
  geom_point(alpha = 0.8, size = 2.5) + 
  theme_bw(base_size = 14) +
  theme(axis.text.x=element_text(angle=45))

ggplot(data = tax.sizes, aes(y = tax.sizes, Supergroup)) + 
  geom_bar(stat = 'identity') 

# Now plot the total set along with with core taxon set (those taxonomic groups 
# with at least 5 species included in the comparative set).
spp.per.og.plts <-
  plot_grid(spp.per.og.plt, core.spp.og.plt, 
            align = 'v', nrow = 2, label_size = 16,
            labels = c('a.', 'b.'))

ggsave(spp.per.og.plts, height = 8, width = 10,
       filename = 'NumSpeciesInOG-NumOrthogroups-FreqDists.pdf')
ggsave(spp.per.og.plts, height = 8, width = 10, dpi = 600,
       filename = 'NumSpeciesInOG-NumOrthogroups-FreqDists.png')
# Generate a dissimilarity matrix using the count data
count.dists <- dist(t(counts.mat[-c(1, ncol(counts.mat))]))
count.dists.tax1 <- count.dists
count.dists.tax2 <- count.dists
attributes(count.dists.tax1)$Labels <- metadat$Taxogroup1_UniEuk
attributes(count.dists.tax2)$Labels <- metadat$Taxogroup2_UniEuk

dist.clusts <- hclust(as.dist(count.dists))
dist.clusts.tax1 <- hclust(as.dist(count.dists.tax1))
dist.clusts.tax2 <- hclust(as.dist(count.dists.tax2))
plot(dist.clusts, cex = 0.5)
plot(dist.clusts.tax1, cex = 0.5)
plot(dist.clusts.tax2, cex = 0.5)
# 
# # Apply t-SNE to these count data - the approach should be less sensitive to outliers than PCA. 
# dt <- t(counts.mat[-c(1, ncol(counts.mat))])
# colnames(dt) <- counts.mat$Orthogroup
# og.tsne.50 <- tsne(t(counts.mat[-c(1, ncol(counts.mat))]), perplex = 50, check_duplicates=FALSE)
# og.tsne.100 <- tsne(t(counts.mat[-c(1, ncol(counts.mat))]), perplex = 100)
# og.tsne.500 <- tsne(t(counts.mat[-c(1, ncol(counts.mat))]), perplex = 500)


# Conduct a PCA on the full count matrix
dt <- t(counts.mat[-c(1, ncol(counts.mat))])
count.pca <- prcomp(dt, rank. = 20, center = T, scale = T)
pcs <- count.pca$x
pcs <- cbind(pcs, metadat)

ggplot(data = pcs, aes(x = PC9, y = PC10, fill = Supergroup_UniEuk)) + 
  geom_point(pch = 21, alpha = 0.5) +
  theme_bw(base_size = 16)

# Calculate the proportion of variance explained by each PC. 
pve <- count.pca$sdev^2 / sum(count.pca$sdev^2)

loadings <- count.pca$rotation
dt <- scale(x = dt) %>% data.table()

pve <- sapply(1:50, function(i){
  proto_x <- sweep(dt, MARGIN = 2, loadings[,i], "*")
  pc_x = apply(proto_x, 1, sum)
  sum(pc_x^2)
})

barplot(pve)


pve <- pve/sum(dt^2)
pve
count.pca.var <- count.pca$sdev^2
pve <- count.pca.var / sum(count.pca.var)
loadings <- count.pca$rotation
sumvar <- sum(apply(t(counts.mat[-c(1, ncol(counts.mat)-1)])^2, 2, sum))
pc.pve <- apply((as.matrix(t(counts.mat[-c(1, ncol(counts.mat)-1)])) %*% loadings)^2, 2, sum) / sumvar


# Lets make a single big dataframe of ortholog gene counts, 
# with species and taxonomic assignments as ID variables
counts <- melt(counts, id.vars = "Orthogroup", variable_name = 'ID')


counts <- merge(counts, metadat, by = 'ID')

ggplot(data = counts, aes(x = Total, fill = Taxogroup1_UniEuk)) + 
  geom_histogram(alpha = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  annotation_logticks() + 
  theme_bw(base_size = 16)




