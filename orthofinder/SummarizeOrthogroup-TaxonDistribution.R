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
gitdir <- '~/ArcadiaScience/github/protein-structural-evolution/orthofinder/'

setwd(gitdir)
source('ArcadiaColorBrewer.R')

# A General plotting function to generate many discrete colors
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

# Read in metadata - we want to know the number of species per supergroup
metadat <- read_tsv('../../../Protein-SeqStruct-Evolution/EukProt-V3-Data/TCS_EukProt_included_data_sets.v03.2021_11_22.txt')
supgrp.spp.counts <- summary(as.factor(metadat$Supergroup_UniEuk))

# Read in the summaries
spp.per.og <- read_tsv(paste0(gitdir, 'og_taxonomy_genecount_summaries/NumSpecies_PerSupergroup_InOGs.tsv'))
count.per.og <- read_tsv(paste0(gitdir, 'og_taxonomy_genecount_summaries/NumGeneCopies_PerSupergroup-InOGs.tsv'))
prop.count.per.og <- read_tsv(paste0(gitdir, 'og_taxonomy_genecount_summaries/PropTotalGeneCopies_PerSupergroup-InOGs.tsv'))

# How many orthogroups does each Supergroup have proteins assigned to?
ogs.per.grp <- data.frame(colSums(spp.per.og[,-c(1:3)] > 0))
ogs.per.grp$Supergroup <- rownames(ogs.per.grp)
colnames(ogs.per.grp)[1] <- "NumOGs"

spp.per.og$NumGrps <- rowSums(spp.per.og[,-c(1:3)] > 0)

# The number of species per orthogroup
spp.per.og.plt <- 
  ggplot(data = spp.per.og, aes(x = NumSppInOG, y = TotalNumGenes, color = NumGrps)) + 
  geom_point(size = 3, stroke = 0.1, alpha = 0.6) + 
  scale_y_log10() + 
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) + 
  xlab("# of Species in OG") +
  ylab("Total # of Gene Copies") + 
  scale_color_viridis_c(option = 'B') + 
  labs(color = "# of Supergroups") +
  theme(legend.position = 'top')


# Summarize how many orthogroups include each number of species
spp.og.hist.dat <- data.frame(NumSppInOG = unique(factor(spp.per.og$NumSppInOG)),
                              NumOGs = NA, NumGrps = NA)
spp.og.hist.dat <- spp.og.hist.dat[order(spp.og.hist.dat$NumSppInOG),]
for(i in 1:196){
  spp.og.hist.dat$NumOGs[i] <- length(which(spp.per.og$NumSppInOG == i))
  spp.og.hist.dat$NumGrps[i] <- mean(spp.per.og$NumGrps[which(spp.per.og$NumSppInOG == i)])
  
}
spp.og.hist.dat$NumSppInOG <- as.numeric(spp.og.hist.dat$NumSppInOG)
spp.per.og.hist <- 
  ggplot(data = spp.og.hist.dat, 
         aes(x = NumSppInOG, 
             y = NumOGs, 
             fill = NumGrps)) + 
  geom_point(pch = 21, alpha = 0.8, size = 3) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  scale_fill_viridis_c(option = 'B') +
  xlab("# of Species in OG") +
  ylab("# of Orthogroups") +
  theme_bw(base_size = 14) +
  labs(fill = "# of Supergroups in OG") +
  theme(legend.position = 'top')

# Summarize how many orthogroups include each number of supergroups
grp.og.hist.dat <- data.frame(NumGrpsInOG = unique(factor(spp.per.og$NumGrps)),
                              NumOGs = NA, NumSppInOG = NA)
grp.og.hist.dat <- grp.og.hist.dat[order(grp.og.hist.dat$NumGrpsInOG),]
for(i in 1:31){
  grp.og.hist.dat$NumOGs[i] <- length(which(spp.per.og$NumGrps == i))
  grp.og.hist.dat$NumSppInOG[i] <- mean(spp.per.og$NumSppInOG[which(spp.per.og$NumGrps == i)])
  
}
grp.og.hist.dat$NumGrpsInOG <- as.numeric(grp.og.hist.dat$NumGrpsInOG)
grp.per.og.hist <- 
  ggplot(data = grp.og.hist.dat, 
         aes(x = NumGrpsInOG, 
             y = NumOGs, 
             fill = NumSppInOG)) + 
  geom_point(pch = 21, alpha = 0.8, size = 3) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  scale_fill_viridis_c(option = 'B') +
  xlab("# of Supergroups in OG") +
  ylab("# of Orthogroups") +
  theme_bw(base_size = 14) +
  labs(fill = "# of Species in OG") +
  theme(legend.position = 'top')

# The relationship between number of gene copies, and the number of supergroups in that orthogroup 
grps.per.og.plt <- 
  ggplot(data = spp.per.og, aes(x = NumGrps, y = TotalNumGenes, color = NumSppInOG)) + 
  geom_point(size = 3, stroke = 0.1, alpha = 0.6) + 
  scale_y_log10() + 
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) + 
  xlab("Number of Supergroups in OG") +
  ylab("Total Number of Gene Copies") + 
  scale_color_viridis_c(option = 'B') + 
  labs(color = "# of Species") +
  theme(legend.position = 'top')

# The number of orthogroups per supergroup
ogs.per.grp <- ogs.per.grp[order(ogs.per.grp$NumOGs),]
ogs.per.grp$Supergroup <- factor(ogs.per.grp$Supergroup, levels = ogs.per.grp$Supergroup)
ogs.per.grp.plt <- 
  ggplot(data = ogs.per.grp, aes(x = Supergroup, y = NumOGs, fill = Supergroup, group = Supergroup)) + 
  geom_point(pch = 21, size = 5) + 
  scale_y_log10() + 
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) + 
  xlab("Supergroup") +
  ylab("Number of OGs") + 
  scale_fill_manual(values = getPalette(31)) +
  labs(fill = "Supergroup") +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())

# Now what if we account for the number of species in each orthogroup?
ogs.per.grp$SppInSupergroup <- 
  supgrp.spp.counts[match(ogs.per.grp$Supergroup, 
                          names(supgrp.spp.counts))]

num.spp.in.grp.v.num.ogs.plt <- 
  ggplot(data = ogs.per.grp, aes(x = SppInSupergroup, y = NumOGs, fill = Supergroup)) + 
  geom_point(pch = 21, size = 5, alpha = 0.8) + 
  scale_y_log10() + 
  annotation_logticks(sides = 'l') + 
  theme_bw(base_size = 14) + 
  xlab("Number of Species in Supergroup") +
  ylab("Number of OGs") + 
  scale_fill_manual(values = getPalette(31)) +
  labs(fill = "Supergroup") +
  theme(legend.position = 'top')

# Why don't we fit a linear model that accounts for multiple things
# predict number of orthogroups for each supergroup, given...
# The number of species in that supergroup 
# The total number of proteins in species within that supergroup
# The age of the supergroup
ogs.per.grp.lm <- lm(NumOGs ~ SppInSupergroup, data = ogs.per.grp)
jtools::plot_summs(ogs.per.grp.lm, plot.distributions = TRUE)


og.count.summ.plts <- 
  plot_grid(spp.per.og.plt, grps.per.og.plt, 
            spp.per.og.hist, grp.per.og.hist,
            labels = c('a.', 'b.', 'c.', 'd.'), 
            nrow = 2, align = 'hv')

ggsave(og.count.summ.plts, file = 'Orthogroup-TaxonomicDistribution-Summaries.pdf',
       height = 12, width = 12)
ggsave(og.count.summ.plts, file = 'Orthogroup-TaxonomicDistribution-Summaries.png',
       height = 12, width = 12, dpi = 600)

ggsave(ogs.per.grp.plt, file = 'NumOrthogroups-PerSupergroup.pdf', 
       height = 4, width = 8)
ggsave(ogs.per.grp.plt, file = 'NumOrthogroups-PerSupergroup.png', 
       height = 4, width = 8, dpi = 600)


# Cumulative distribution, all together, and for each supergroup, of the percent of proteins 
# included in orthogroups that include all species, down to just 1. 
cum.dist.og.size <-
  data.frame(NumSppInOG = 196:1, PercentOfAllProteins = NA)
cum.dist.og.size <- cbind(cum.dist.og.size, prop.count.per.og[1:196,c(4:ncol(prop.count.per.og))])
cum.dist.og.size[,3:ncol(cum.dist.og.size)] <- NA

total.num.prots <- colSums(count.per.og[,c(2,4:ncol(count.per.og))])
for(i in 1:196){
  og.set <- count.per.og[which(count.per.og$NumSppInOG >= cum.dist.og.size[i,1]),]
  
  cum.dist.og.size[i,-1] <- 
    colSums(og.set[,-c(1,3)]) / total.num.prots
}

cum.dist.og.size <- 
  melt(cum.dist.og.size, id.vars = 'NumSppInOG')
cum.dist.og.size.all <- cum.dist.og.size[which(cum.dist.og.size$variable == 'PercentOfAllProteins'),]
cum.dist.og.size <- cum.dist.og.size[-which(cum.dist.og.size$variable == 'PercentOfAllProteins'),]


og.spp.inclusion.plt <- 
  ggplot(data = cum.dist.og.size, aes(x = NumSppInOG, y = value, color = variable)) + 
  geom_line(size = 1.5) + 
  geom_line(data = cum.dist.og.size.all, aes(x = NumSppInOG, y = value),
            size = 3, color = 'black') + 
  geom_vline(xintercept = 5, color = 'blue', lty = 2, size = 1, alpha = 0.5) +
  geom_vline(xintercept = 10, color = 'red', lty = 2, size = 1, alpha = 0.5) +
  ylab('Cumul. Prop. of Orthologs') +
  xlab('Number of Species in Orthogroup') +
  scale_color_manual(values = getPalette(31)) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'top') +
  labs(color = "Supergroup")


ggsave(og.spp.inclusion.plt, file = 'CumulDist-PropOrthologsIn-N-SppOGs.pdf', 
       height = 8, width = 10)
ggsave(og.spp.inclusion.plt, file = 'CumulDist-PropOrthologsIn-N-SppOGs.png', 
       height = 8, width = 10, dpi = 600)

# And another plot of the # of OGs included for each supergroup when using cutoffs of 4 or 10 spp per og. 
og.cutoff.effect.per.grp <- 
  data.frame(Supergroup = unique(cum.dist.og.size$variable),
             UniqueSpp_Cut4 = NA, 
             UniqueSpp_Cut10 = NA)

grps <- unique(cum.dist.og.size$variable)
for(grp in 1:length(grps)){
  grp.props <- cum.dist.og.size[which(cum.dist.og.size$variable == grps[grp]),]
  og.cutoff.effect.per.grp$UniqueSpp_Cut4[grp] <- 
    grp.props$value[which(grp.props$NumSppInOG == 4)]
  og.cutoff.effect.per.grp$UniqueSpp_Cut10[grp] <- 
    grp.props$value[which(grp.props$NumSppInOG == 10)]
}

# Reorder to be increasing
og.cutoff.effect.per.grp$Rank <- rank(og.cutoff.effect.per.grp$UniqueSpp_Cut4)
og.cutoff.effect.per.grp <- melt(og.cutoff.effect.per.grp, id.vars = c('Supergroup', 'Rank'))
og.cutoff.effect.per.grp <- 
  og.cutoff.effect.per.grp[order(og.cutoff.effect.per.grp$Rank),]
og.cutoff.effect.per.grp$Supergroup <- 
  factor(og.cutoff.effect.per.grp$Supergroup, 
         levels = unique(og.cutoff.effect.per.grp$Supergroup[order(og.cutoff.effect.per.grp$Rank)]))

og.cutoff.effect.per.grp.plt <- 
  ggplot(data = og.cutoff.effect.per.grp,
         aes(x = Supergroup, y = value, fill = variable)) + 
  geom_point(pch = 21, size = 3) + 
  scale_fill_manual(values = unname(arcadia.pal(n = 3, name = 'Accent'))) +
  theme_bw(base_size = 14) + 
  xlab("Supergroup") +
  ylab("Proportion of OGs\nRetained at Cutoff") + 
  labs(fill = "Cutoff") +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())

ggsave(og.cutoff.effect.per.grp.plt,
       file = 'PropOGs_Retained_PerCutoff_PerSupergroup.pdf',
       height = 6, width = 6)
ggsave(og.cutoff.effect.per.grp.plt,
       file = 'PropOGs_Retained_PerCutoff_PerSupergroup.png',
       height = 6, width = 6, dpi = 600)

####################################################################
# Run UMAP

run.umap <- function(count.set, umap.config, is.prop = F){
  if(is.prop == F){
    count.set[,-1] <- center_scale(count.set[,-1])
  }
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

umap.config <- umap.defaults
umap.config$n_neighbors <- 25
#umap.config$knn_repeats <- 10
#umap.config$metric <- 'euclidean'
umap.config$n_epochs <- 1000
umap.config$min_dist <- 0.25
umap.config$spread <- 1.75
#umap.config$verbose <- TRUE

# Umap using the cutoff of ten species minimum per orthogroup
core.counts <- count.per.og[which(count.per.og$NumSppInOG >= 10),]
per.tax.count.umap <- 
  run.umap(core.counts[,-c(2:3)], umap.config)
per.tax.count.umap$umap.plot <-
  ggplot(data = per.tax.count.umap$umap.df, aes(x = UMAP_Axis1, y = UMAP_Axis2,
                                          color = core.counts$NumSppInOG)) +
  geom_point(size = 0.75, alpha = 0.5) +
  scale_color_met_c('Hiroshige', direction = -1) +
  theme_bw(base_size = 14) +
  xlab('UMAP Axis 1') +
  ylab('UMAP Axis 2') +
  labs(color = "Number of Species\nin Orthogroup") +
  theme(legend.position = 'bottom')
per.tax.count.umap$umap.plot

ggsave(per.tax.count.umap$umap.plot,
       file = "PerSupergroup_GeneCount_UMAP.pdf",
       height = 8, width = 8)
ggsave(per.tax.count.umap$umap.plot,
       file = "PerSupergroup_GeneCount_UMAP.png",
       height = 8, width = 8)


############################################################################
# What do things look like when we consider orthogroups that have on average <= 3 gene copies, 
# per species, and include at least 10 species?
core.ogs <- 
  spp.per.og$Orthogroup[which(c(spp.per.og$TotalNumGenes / spp.per.og$NumSppInOG) >= 1 & 
                                c(spp.per.og$TotalNumGenes / spp.per.og$NumSppInOG) <= 3 & 
                                spp.per.og$NumSppInOG >= 10)] 

# Read in the summaries
spp.per.core.og <- spp.per.og[which(spp.per.og$Orthogroup %in% core.ogs),]
count.per.core.og <- count.per.og[which(count.per.og$Orthogroup %in% core.ogs),]
prop.count.per.core.og <- prop.count.per.og[which(prop.count.per.og$Orthogroup %in% core.ogs),]

# The number of gene copies per orthogroup size (number of species included)
spp.per.core.og.plt <- 
  ggplot(data = spp.per.core.og, aes(x = NumSppInOG, y = TotalNumGenes, color = NumGrps)) + 
  geom_point(size = 3, stroke = 0.1, alpha = 0.75) + 
  theme_bw(base_size = 14) + 
  xlab("# of Species in OG") +
  ylab("Total # of Gene Copies") + 
  scale_color_viridis_c(option = 'B') + 
  labs(color = "# of Supergroups") +
  theme(legend.position = 'top')


# Summarize how many orthogroups include each number of species
spp.core.og.hist.dat <- data.frame(NumSppInOG = unique(factor(spp.per.core.og$NumSppInOG)),
                              NumOGs = NA, NumGrps = NA)
spp.core.og.hist.dat <- spp.core.og.hist.dat[order(spp.core.og.hist.dat$NumSppInOG),]
for(i in 1:nrow(spp.core.og.hist.dat)){
  spp.core.og.hist.dat$NumOGs[i] <- length(which(spp.per.core.og$NumSppInOG == i+9))
  spp.core.og.hist.dat$NumGrps[i] <- mean(spp.per.core.og$NumGrps[which(spp.per.core.og$NumSppInOG == i+9)])
}
spp.core.og.hist.dat$NumSppInOG <- as.numeric(as.character(spp.core.og.hist.dat$NumSppInOG))
spp.per.core.og.hist <- 
  ggplot(data = spp.core.og.hist.dat, 
         aes(x = NumSppInOG, 
             y = NumOGs, 
             fill = NumGrps)) + 
  geom_point(pch = 21, alpha = 0.8, size = 3) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  scale_fill_viridis_c(option = 'B') +
  xlab("# of Species in OG") +
  ylab("# of Orthogroups") +
  theme_bw(base_size = 14) +
  labs(fill = "# of Supergroups in OG") +
  theme(legend.position = 'top')

