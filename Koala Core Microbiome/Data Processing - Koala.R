#Load Library
library(tidyverse)
library(biomformat)

# Load your files
taxonomy <- read_delim(file = "taxonomy.tsv", delim="\t")
feature_table <- read_tsv("feature-table.txt", skip = 1)
metadata <- read_tsv("Koala_Metadata.tsv")

#Standardize Files
colnames(feature_table)[1] <- "Feature.ID"
colnames(taxonomy)[1] <- "Feature.ID"
colnames(metadata)[1] <- "SampleID"
#Separate and filter
taxonomy <- taxonomy %>%
  separate(Taxon,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";",
           fill = "right")

feature_long <- feature_table %>%
  pivot_longer(
    cols = -Feature.ID,
    names_to = "SampleID",
    values_to = "Count"
  )

merged <- feature_long %>%
  left_join(metadata, by = "SampleID") %>%
  left_join(taxonomy, by = "Feature.ID")

#Might have to change Age variable
taxa_abundance <- merged %>%
  group_by(Age, Genus) %>%
  summarise(TotalAbundance = sum(Count), .groups = "drop")

taxa_abundance <- taxa_abundance %>%
  group_by(Age) %>%
  mutate(RelAbundance = TotalAbundance / sum(TotalAbundance))

top_taxa <- taxa_abundance %>%
  group_by(Age) %>%
  filter(RelAbundance >= quantile(RelAbundance, 0.90))

top_taxa %>%
  arrange(Age, desc(RelAbundance))

top_taxa_table <- top_taxa %>%
  select(Age, Genus, RelAbundance) %>%
  pivot_wider(
    names_from = Age,
    values_from = RelAbundance
  )
write.csv(top_taxa_table, "top_10_percent_taxa.csv", row.names = FALSE)
