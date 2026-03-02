#Load Libraries
library(tidyverse)
library(janitor)
library(broom)
library(patchwork)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)


#Functions

# Extract species name from filename
get_species_from_filename <- function(path) {
  nm <- basename(path)
  sp <- str_split(nm, "_top_", simplify = TRUE)[1, 1]
  sp <- str_replace_all(sp, "-", "_")
  str_to_title(sp)
}

# Detect maturity from filename if needed
get_maturity_hint_from_filename <- function(path) {
  nm <- tolower(basename(path))
  if (str_detect(nm, "immature")) return("Immature")
  if (str_detect(nm, "mature"))   return("Mature")
  return(NA_character_)
}

# Find genus column robustly
pick_genus_col <- function(df) {
  cn <- names(df)
  if ("genus" %in% cn) return("genus")
  genus_like <- cn[str_detect(cn, "genus")]
  if (length(genus_like) >= 1) return(genus_like[1])
  return(NA_character_)
}
# Read one CSV and reshape into long format
read_one_file <- function(path) {
  
  raw <- readr::read_csv(path, show_col_types = FALSE) %>%
    clean_names()
  
  species <- get_species_from_filename(path)
  maturity_hint <- get_maturity_hint_from_filename(path)
  genus_col <- pick_genus_col(raw)
  
  if (is.na(genus_col)) {
    stop("No genus column found in: ", basename(path))
  }
  
  cn <- names(raw)
  has_imm <- "immature" %in% cn
  has_mat <- "mature" %in% cn
  
  # Case 1: Has explicit Immature/Mature columns
  if (has_imm || has_mat) {
    long <- raw %>%
      transmute(
        genus = .data[[genus_col]],
        immature = if (has_imm) as.numeric(.data[["immature"]]) else NA_real_,
        mature   = if (has_mat) as.numeric(.data[["mature"]]) else NA_real_
      ) %>%
      pivot_longer(c(immature, mature),
                   names_to = "maturity",
                   values_to = "abundance") %>%
      mutate(
        maturity = str_to_title(maturity),
        species = species
      )
    return(long)
  }
  
  # Case 2: One numeric abundance column only
  numeric_cols <- raw %>%
    select(where(is.numeric)) %>%
    names()
  
  if (length(numeric_cols) == 1) {
    long <- raw %>%
      transmute(
        genus = .data[[genus_col]],
        abundance = as.numeric(.data[[numeric_cols[1]]]),
        maturity = maturity_hint %||% "Unknown",
        species = species
      )
    return(long)
  }
  
  # Case 3: Taxonomy only (no abundance)
  raw %>%
    transmute(
      genus = .data[[genus_col]],
      abundance = NA_real_,
      maturity = maturity_hint %||% "Unknown",
      species = species
    )
}
files <- list.files(pattern = "\\.csv$", full.names = TRUE)

combined <- purrr::map_dfr(files, read_one_file) %>%
  mutate(
    genus = as.character(genus),
    maturity = factor(maturity, levels = c("Immature","Mature","Unknown"))
  ) %>%
  filter(!is.na(genus), genus != "")

write_csv(combined, "combined_genus_maturity.csv")

# Statistics

abund_df <- combined %>%
  filter(maturity %in% c("Immature","Mature")) %>%
  filter(!is.na(abundance)) %>%
  mutate(
    log_abund = log10(abundance + 1e-6)
  )
## Paired Test
paired_df <- abund_df %>%
  select(species, genus, maturity, log_abund) %>%
  distinct() %>%
  pivot_wider(names_from = maturity,
              values_from = log_abund) %>%
  filter(!is.na(Immature) & !is.na(Mature)) %>%
  mutate(diff = Mature - Immature)

genus_stats <- paired_df %>%
  group_by(genus) %>%
  summarise(
    n_species = n(),
    mean_diff = mean(diff),
    median_diff = median(diff),
    p_value = if (n_species >= 3) t.test(diff)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      is.na(p_adj) ~ "Insufficient pairs",
      p_adj < 0.05 & mean_diff > 0 ~ "More common in Mature",
      p_adj < 0.05 & mean_diff < 0 ~ "More common in Immature",
      TRUE ~ "No clear difference"
    )
  ) %>%
  arrange(p_adj, desc(abs(mean_diff)))

write_csv(genus_stats, "genus_maturity_stats.csv")


## Different Model Accounting For Lack of Pairs

model_df <- abund_df %>%
  filter(maturity %in% c("Immature","Mature")) %>%
  drop_na(log_abund)

model_df$maturity <- factor(model_df$maturity, levels = c("Immature","Mature"))
model_df$species  <- factor(model_df$species)
model_df$genus    <- factor(model_df$genus)

model <- lmer(
  log_abund ~ maturity + (1 | species) + (1 | genus),
  data = model_df
)

summary(model)

p_by_species <- ggplot(model_df, aes(x = maturity, y = log_abund)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  facet_wrap(~ species) +
  labs(
    title = "Raw log-abundance by maturity within each species",
    subtitle = "Shows whether maturity shifts are consistent across species",
    x = NULL,
    y = "log10(abundance + pseudocount)"
  ) +
  theme_minimal()

p_by_species

## Wilcox Test
wilcox.test(log_abund ~ maturity, data = model_df)
wtest <- wilcox.test(log_abund ~ maturity, data = model_df)

p_w <- signif(wtest$p.value, 3)

p_wilcox <- ggplot(model_df, aes(x = maturity, y = log_abund)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, alpha = 0.25) +
  labs(
    title = "Pooled log-abundance by maturity (Wilcoxon test)",
    subtitle = paste0("Wilcoxon rank-sum p = ", p_w,
                      " (pooled; ignores species/genus)"),
    x = NULL,
    y = "log10(abundance + pseudocount)"
  ) +
  theme_minimal()

p_wilcox

## General BoxPlot

ggplot(model_df, aes(x = maturity, y = log_abund)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.4) +
  labs(
    title = "Abundance Differences by Maturity (Mixed Effects Model)",
    y = "log10(Abundance)"
  ) +
  theme_minimal()
