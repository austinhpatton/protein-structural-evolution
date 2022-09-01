# Summaries of orthogroups and comparative genomic statistics for The Comparative SeT
# as inferred using OrthoFinder
library(tidyverse)
library(reshape2)
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
                              NumOGs = NA)
spp.og.hist.dat <- spp.og.hist.dat[order(spp.og.hist.dat$NumSppInOG),]
for(i in 1:196){
  spp.og.hist.dat$NumOGs[i] <- length(which(spp.per.og$NumSppInOG == i))
}
spp.og.hist.dat$NumSppInOG <- as.numeric(spp.og.hist.dat$NumSppInOG)
spp.per.og.hist <- 
  ggplot(data = spp.og.hist.dat, 
         aes(x = NumSppInOG, 
             y = NumOGs, 
             fill = NumSppInOG)) + 
  geom_point(pch = 21, alpha = 0.8, size = 3) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  scale_fill_viridis_c(option = 'B') +
  xlab("# of Species in OG") +
  ylab("Total # of Gene Copies") +
  theme_bw(base_size = 14) +
  labs(fill = "# of Species") +
  theme(legend.position = 'top')

# Summarize how many orthogroups include each number of supergroups
grp.og.hist.dat <- data.frame(NumGrpsInOG = unique(factor(spp.per.og$NumGrps)),
                              NumOGs = NA)
grp.og.hist.dat <- grp.og.hist.dat[order(grp.og.hist.dat$NumGrpsInOG),]
for(i in 1:31){
  grp.og.hist.dat$NumOGs[i] <- length(which(spp.per.og$NumGrps == i))
}
grp.og.hist.dat$NumGrpsInOG <- as.numeric(grp.og.hist.dat$NumGrpsInOG)
grp.per.og.hist <- 
  ggplot(data = grp.og.hist.dat, 
         aes(x = NumGrpsInOG, 
             y = NumOGs, 
             fill = NumGrpsInOG)) + 
  geom_point(pch = 21, alpha = 0.8, size = 3) +
  scale_y_log10() +
  annotation_logticks(sides = 'l') + 
  scale_fill_viridis_c(option = 'B') +
  xlab("# of Supergroups in OG") +
  ylab("# of Orthogroups") +
  theme_bw(base_size = 14) +
  labs(fill = "# of Supergroups") +
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
  scale_color_viridis_d(option = 'B') +
  labs(fill = "Supergroup") +
  theme(legend.position = 'none', 
        axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank())

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

getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))

og.spp.inclusion.plt <- 
  ggplot(data = cum.dist.og.size, aes(x = NumSppInOG, y = value, color = variable)) + 
  geom_line(size = 1.5) + 
  geom_line(data = cum.dist.og.size.all, aes(x = NumSppInOG, y = value),
            size = 3, color = 'black') + 
  scale_x_reverse() + 
  xlab('Cumul. Prop. of Orthologs') +
  ylab('Number of Species in Orthogroup') +
  scale_color_manual(values = getPalette(31)) +
  theme_bw(base_size = 14) +
  theme(legend.position = 'top') +
  labs(color = "Supergroup")

ggsave(og.spp.inclusion.plt, file = 'CumulDist-PropOrthologsIn-N-SppOGs.pdf', 
       height = 8, width = 10)
ggsave(og.spp.inclusion.plt, file = 'CumulDist-PropOrthologsIn-N-SppOGs.png', 
       height = 8, width = 10, dpi = 600)
