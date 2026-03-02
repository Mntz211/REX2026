#Load Library
library(tidyverse)
library(biomformat)
# Load your files
taxonomy <- read_delim(file = "taxonomy.tsv", delim="\t")
feature_table <- read_tsv("otutable.tsv", skip = 1)
metadata <- read_tsv("Gibbon_metadata.tsv")

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
  group_by(AGE, Genus) %>%
  summarise(TotalAbundance = sum(Count), .groups = "drop")

taxa_abundance <- taxa_abundance %>%
  group_by(AGE) %>%
  mutate(RelAbundance = TotalAbundance / sum(TotalAbundance))

top_taxa <- taxa_abundance %>%
  group_by(AGE) %>%
  filter(RelAbundance >= quantile(RelAbundance, 0.90))

top_taxa %>%
  arrange(AGE, desc(RelAbundance))

top_taxa_table <- top_taxa %>%
  select(AGE, Genus, RelAbundance) %>%
  pivot_wider(
    names_from = AGE,
    values_from = RelAbundance
  )
write.csv(top_taxa_table, "adultgibbons_top_10_percent_taxa.csv", row.names = FALSE)
