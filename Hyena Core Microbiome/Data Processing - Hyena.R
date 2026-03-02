#Load Library
library(tidyverse)
library(biomformat)

# Load your files
taxonomy <- read.csv(file = "hyenataxonomy.csv")
feature_table <- load("hyenatable.Rdata")
metadata <- read_csv("hyenametadata.csv")

#Standardize Files
asvf <- cbind("Feature.ID" = rownames(asvf), asvf)
rownames(asvf) <- NULL
colnames(taxonomy)[1] <- "Feature.ID"
colnames(metadata)[1] <- "SampleID"
#Separate and filter
taxonomy <- taxonomy %>%
  separate(Taxon,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";",
           fill = "right")

feature_long <- asvf %>%
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
  group_by(age_yrs, Genus) %>%
  summarise(TotalAbundance = sum(Count), .groups = "drop")

taxa_abundance <- taxa_abundance %>%
  group_by(age_yrs) %>%
  mutate(RelAbundance = TotalAbundance / sum(TotalAbundance))

top_taxa <- taxa_abundance %>%
  group_by(age_yrs) %>%
  filter(RelAbundance >= quantile(RelAbundance, 0.90))

top_taxa %>%
  arrange(age_yrs, desc(RelAbundance))

top_taxa_table <- top_taxa %>%
  select(age_yrs, Genus, RelAbundance) %>%
  pivot_wider(
    names_from = age_yrs,
    values_from = RelAbundance
  )
write.csv(top_taxa_table, "hyena_top_10_percent_taxa.csv", row.names = FALSE)
