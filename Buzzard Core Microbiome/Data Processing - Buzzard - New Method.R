
#Load Libraries 
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library(ggVennDiagram)

#Create phyloseq object
buzzard_mpt <- readRDS("16S_phyloseq.rds")

#Convert to relative abundance 
mpt_RA_buzzard <- transform_sample_counts(buzzard_mpt, fun=function(x) x/sum(x))


# Set a prevalence threshold and abundance threshold. 
# Justification:Detection at 0.001 to filter out rare taxa that are irrelavent to our "core microbiome". 
# Prevalence at 0.8 to include taxa that appear in 80% of samples, thus appearing in most samples making it fair to call the "core microbiome"
immature_ASVs <- core_members(mpt_RA_buzzard, detection=0, prevalence = 0.7)


# What are these ASVs? you can code it in two different ways to see the same things

mostabundantimmature<-tax_table(prune_taxa(immature_ASVs,buzzard_mpt))

# can plot those ASVs' relative abundance

prune_taxa(immature_ASVs,buzzard_mpt) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`age_days`, scales ="free")

#Save table
write.csv(mostabundantimmature, "top_10_percent_taxa_immature.csv", row.names = FALSE)

