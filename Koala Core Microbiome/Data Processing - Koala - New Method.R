
#Load Libraries 
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(microbiome)
library(ggVennDiagram)
#Load data
# Change file paths as necessary
metafp <- "Koala_Metadata.tsv"
meta <- read_delim(metafp, delim="\t")

otufp <- "feature-table.txt"
otu <- read_delim(file = otufp, delim="\t", skip=1)

taxfp <- "taxonomy.tsv"
tax <- read_delim(taxfp, delim="\t")

phylotreefp <- "tree.nwk"
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
tax_mat <- tax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]
# Make sampleids the rownames
rownames(tax_mat) <- tax$Feature.ID
# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)

#Create phyloseq object
koala_mpt <- phyloseq(OTU, SAMP, TAX, phylotree)

#Convert to relative abundance 
mpt_RA_koala <- transform_sample_counts(koala_mpt, fun=function(x) x/sum(x))
# Subset dataset into treatment and control groups
# Change based on age, look at metadata
mature_stat <- subset_samples(mpt_RA_koala, Age=="Mature")
immature_stat <- subset_samples(mpt_RA_koala, Age=="Immature")

# Set a prevalence threshold and abundance threshold. 
# Justification:Detection at 0.001 to filter out rare taxa that are irrelavent to our "core microbiome". 
# Prevalence at 0.8 to include taxa that appear in 80% of samples, thus appearing in most samples making it fair to call the "core microbiome"
mature_ASVs <- core_members(mature_stat, detection=0, prevalence = 0.7)
immature_ASVs <- core_members(immature_stat, detection=0, prevalence = 0.7)

# Make a Venn-diagram
age_venn <- ggVennDiagram(x=list(Immature = immature_ASVs, Mature = mature_ASVs))

# Make sure that you have a line that saves the Venn diagram as a png and this file is present within your project folder
ggsave("Venn_Age.png", age_venn)


##Curiosity##

# What are these ASVs? you can code it in two different ways to see the same things
mostabundantmature<- prune_taxa(mature_ASVs,koala_mpt) %>%
  tax_table()

mostabundantimmature<-tax_table(prune_taxa(immature_ASVs,koala_mpt))

# can plot those ASVs' relative abundance
prune_taxa(mature_ASVs,koala_mpt) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Age`, scales ="free")

prune_taxa(immature_ASVs,koala_mpt) %>% 
  plot_bar(fill="Genus") + 
  facet_wrap(.~`Age`, scales ="free")

#Save table
write.csv(mostabundantmature, "top_10_percent_taxa_mature.csv", row.names = FALSE)
write.csv(mostabundantimmature, "top_10_percent_taxa_immature.csv", row.names = FALSE)

