
#Load Libraries 
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library(ggVennDiagram)
#Load data
# Change file paths as necessary
metafp <- "WarblerMetadata.csv"
meta <- read_csv(metafp)

otufp <- "WarblerTable.csv"
otu <- read_csv(otufp)

taxfp <- "WarblerTaxonomy.csv"
tax <- read_csv(taxfp)

phylotreefp <- "warblertree.nwk"
phylotree <- read.tree(phylotreefp)

#Standardize Files
colnames(otu)[1] <- "Feature.ID"
colnames(tax)[1] <- "Feature.ID"
colnames(meta)[1] <- "SampleID"

#Format OTU table
# save everything except first column (OTU ID) into a matrix
otu_mat <- as.matrix(otu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(otu_mat) <- otu$Feature.ID
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) 
class(OTU)

#Format sample metadata
# Save everything except sampleid as new data frame
samp_df <- as.data.frame(meta[,-1])
# Make sampleids the rownames
rownames(samp_df)<- meta$SampleID
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(samp_df)
class(SAMP)

#Format taxonomy
# Convert taxon strings to a table with separate taxa rank columns
# Might have to change Confidence
tax_mat <- as.matrix(tax)
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$Feature.ID
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#Create phyloseq object
warbler_mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

#Convert to relative abundance 
mpt_RA_warbler <- transform_sample_counts(warbler_mpt, fun=function(x) x/sum(x))
# Subset dataset into treatment and control groups
# Change based on age, look at metadata
fletchling_stat <- subset_samples(mpt_RA_warbler, Ageclass=="FL")
adult_stat <- subset_samples(mpt_RA_warbler, Ageclass=="A")

# Set a prevalence threshold and abundance threshold. 
# Justification:Detection at 0.001 to filter out rare taxa that are irrelavent to our "core microbiome". 
# Prevalence at 0.8 to include taxa that appear in 80% of samples, thus appearing in most samples making it fair to call the "core microbiome"
fletchling_ASVs <- core_members(fletchling_stat, detection=0, prevalence = 0.7)
adult_ASVs <- core_members(adult_stat, detection=0, prevalence = 0.7)

# Make a Venn-diagram
age_venn <- ggVennDiagram(x=list(Immature = fletchling_ASVs, Mature = adult_ASVs))

# Make sure that you have a line that saves the Venn diagram as a png and this file is present within your project folder
ggsave("Warbler_Venn_Age.png", age_venn)


##Curiosity##

# What are these ASVs? you can code it in two different ways to see the same things
mostabundantmature<- prune_taxa(adult_ASVs,warbler_mpt) %>%
  tax_table()

mostabundantimmature<-tax_table(prune_taxa(fletchling_ASVs,warbler_mpt))

# can plot those ASVs' relative abundance
prune_taxa(adult_ASVs,warbler_mpt) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Ageclass`, scales ="free")

prune_taxa(fletchling_ASVs,warbler_mpt) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Ageclass`, scales ="free")

#Save table
write.csv(mostabundantmature, "top_10_percent_taxa_mature.csv", row.names = FALSE)
write.csv(mostabundantimmature, "top_10_percent_taxa_immature.csv", row.names = FALSE)

