#!/bin/bash
conda activate orthofinder

cd ~/environment/analyses/speciesrax 

# Store some paths as variables to make things more legible
generax=/home/ubuntu/environment/software/GeneRax/build/bin/generax
msaDir=/home/ubuntu/environment/analyses/pargenes/Extreme-Core-OGs/extr-core-og-trimmed-msas
treeDir=/home/ubuntu/environment/analyses/speciesrax/ExtremeCoreOG_UnreconciledTrees
mapDir=/home/ubuntu/environment/analyses/speciesrax/ExtremeCoreOG_SppGeneMaps
families=/home/ubuntu/environment/analyses/speciesrax/ExtremeCoreOG.families
outDir=/home/ubuntu/environment/analyses/speciesrax/EukProt_TCS_ExtremeCoreOGs

# If we haven't created the directory of unreconciled trees, do so now
if [ ! -d $treeDir ]
then
    # If we haven't done so yet, create the directory
    mkdir $treeDir

    # Read in the list of trees
    origTrees=$(ls /home/ubuntu/environment/analyses/pargenes/EukProt_TCS_ExtremeCoreOGs/mlsearch_run/results/*/*bestTree)
    
    # And start copying the MSAs over here. 
    for og in $origTrees
    do
        cp $og $treeDir
    done
fi

# If we haven't created the directory of gene-species maps, do so now
if [ ! -d $mapDir ]
then
    # If we haven't done so yet, create the directory
    mkdir $mapDir

    # Read in the list of trees
    trees=$(ls $treeDir)
    
    # And start copying the MSAs over here. 
    for tree in $trees
    do
        # Get the OG name
        og=$(echo $tree | sed "s/-clipkit_fa.raxml.bestTree//g")
        
        # Now pull out the sequences, and split into a TreeRecs format mapping
        # file, where each protein in the tree is a new line, listing species 
        # and then the protein
        grep ">" $msaDir/${og}* | sed "s/>//g" > prot
        sed "s/^[^_]*_//" prot | sed "s/_P.*//g" > spp
        paste prot spp > $mapDir/${og}_SppProtMap.link
        rm prot && rm spp
    done
fi


# Lastly, create the families file if we haven't yet 
if [ ! -f $families ]
then
    # Initialize 
    echo "[FAMILIES]" > $families
    
    # Now loop through, and for each orthogroup write out the 
    # location of the tree file, species-sequence mapping, MSA, 
    # and substitution model
    for og in $(ls $treeDir | sed "s/-clipkit_fa.raxml.bestTree//g")
    do
        subst=/home/ubuntu/environment/analyses/pargenes/EukProt_TCS_ExtremeCoreOGs/mlsearch_run/results/$og-clipkit_fa/$og-clipkit_fa.raxml.bestModel
        echo "- $og" >> $families
        echo "starting_gene_tree = $treeDir/$og-clipkit_fa.raxml.bestTree" >> $families
        echo "mapping = $mapDir/${og}_SppProtMap.link" >> $families
        echo "alignment = $msaDir/$og-clipkit.fa" >> $families
        echo "subst_model = $subst" >> $families
    done
fi

# Infer the species tree
mpiexec -np 4 $generax --families $families --strategy SPR \
--si-strategy HYBRID --species-tree MiniNJ --rec-model UndatedDTL \
--per-family-rates --prune-species-tree --si-estimate-bl \
--si-spr-radius 5 --max-spr-radius 5 --si-quartet-support \
--prefix SpeciesRax