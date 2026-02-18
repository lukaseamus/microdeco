# 1. Prepare data ####
# 1.1 Load ASVs ####
require(tidyverse)
require(magrittr)
require(here)
ASV <- here("DNA", "RDS", "ASV.rds") %>%
  read_rds() %T>%
  print()

# 1.2 Remove contamination ####
# 1.2.1 decontam ####
# decontam needs a logical variable identifying blanks
require(phyloseq)
sample_data(ASV)$Blank <- sample_data(ASV)$Sample == "Blank"
sample_data(ASV)

require(decontam)
contaminant <- isContaminant(
  ASV, method = "prevalence",
  neg = "Blank"
)

contaminant %>%
  count(contaminant)
# No identified contaminants

# 1.2.2 microDecon ####
require(microDecon)
# Reformat
ASV_df <- otu_table(ASV) %>%
  as.data.frame() %>%
  rownames_to_column("File") %>%
  as_tibble() %>%
  mutate(Sample = File %>% 
           str_split_i("_", 1)) %>%
  select(-File) %>%
  pivot_longer(cols = -Sample,
               names_to = "ASV",
               values_to = "Abundance") %>%
  pivot_wider(names_from = Sample,
              values_from = Abundance) %>%
  as.data.frame() %T>%
  print()

# Check numbers
sample_data(ASV) %>%
  as_tibble() %>%
  count(Blank)

ASV_decon <- ASV_df %>%
  decon(numb.blanks = 1, # There is 1 blank
        numb.ind = 64, # There are 64 samples
        taxa = FALSE) %T>%
  print()

ASV_decon$OTUs.removed
# Unlike decontam, microDecon identified the only two balnks ASVs
# as contaminants. microDecon suggests to remove all reads of the 
# two blank ASVs from samples. This is unsurprising when looking
# at reads of the two ASVs in question across samples:
ASV_df %>% filter(Blank > 1) # They're most abundant in Blank

# 1.2.3 Manual removal ####
# Identify blank ASVs
ASV %>%
  subset_samples(Blank) %>%
  prune_taxa(taxa = taxa_sums(.) > 0) %>%
  tax_table()
# These are pretty common, human-associated bacteria.
# Remove them completely from all samples.
ASV %<>%
  prune_taxa(taxa = !taxa_names(ASV) %in% c("ASV176", "ASV756")) %T>%
  print()

# Remove blank and blank indicator
ASV %<>%
  subset_samples(!Blank) %T>%
  print()

sample_data(ASV)$Blank <- NULL
ASV

# 1.3 Remove chloroplasts ####
ASV %>%
  subset_taxa(Order == "Chloroplast") %>%
  tax_table()
# 175 chloroplast ASVs

# Remove
ASV %<>% # Importantly, don't drop NAs
  subset_taxa(is.na(Order) | Order != "Chloroplast") %T>%
  print()

# 1.4 Assign treatments ####
# Add treatment variable
sample_data <- sample_data(ASV) %>%
  data.frame() %>% # Note data.frame() instead of as.data.frame()
  rownames_to_column("File") %>%
  as_tibble() %>%
  mutate(
    Treatment = case_when(
      PAR == 0 ~ "Dark 15°C",
      Temperature == 15 ~ "Light 15°C",
      Temperature == 20 ~ "Light 20°C",
      Temperature == 25 ~ "Light 25°C"
    )
  ) %T>%
  print()

# Assign baseline samples random treatments
set.seed(1)
random <- sample_data %>%
  filter(is.na(Treatment)) %>%
  distinct(Sample) %>%
  mutate(
    Treatment = sample_data %$% Treatment %>% na.omit() %>% unique() %>%
      sample(., 4) %>% c(sample(., 1)),
    Temperature = Treatment %>% str_extract("\\d+") %>% as.numeric(),
    PAR = if_else(Treatment %>% str_detect("Dark"), 0, 8)
  ) %>%
  mutate(
    Tank = case_when(
      Treatment == "Light 15°C" ~ sample(1:4, n()),
      Treatment == "Dark 15°C" ~ sample(5:8, n()),
      Treatment == "Light 20°C" ~ sample(9:12, n()),
      Treatment == "Light 25°C" ~ sample(13:16, n())
    ),
    .by = Treatment
  ) %T>%
  print()

# Add randomly generated treatments to sample_data
sample_data %<>%
  left_join(random, by = "Sample") %>%
  mutate(
    Treatment = coalesce(Treatment.x, Treatment.y),
    Tank = coalesce(Tank.x, Tank.y),
    Temperature = coalesce(Temperature.x, Temperature.y),
    PAR = coalesce(PAR.x, PAR.y)
  ) %>%
  select(-c(ends_with(".x"), ends_with(".y"))) %T>%
  print()

# Update sample_data(ASV)
sample_data(ASV)
sample_data(ASV) <- sample_data %>% column_to_rownames("File")
sample_data(ASV)

# Save decontaminated phyloseq object
ASV %>%
  write_rds(here("DNA", "RDS", "ASV_clean.rds"))

# 1.5 Tidy data ####
# 1.5.1 ASV level ####
# Phyloseq is helpful, but in some cases I'll want all info
# in place. otu_table() is already a species-abundance matrix:
taxa_are_rows(ASV)
otu_table(ASV)

# But I also want a tidy dataframe
ASV_tidy <- otu_table(ASV) %>%
  as.data.frame() %>%
  rownames_to_column("File") %>%
  as_tibble() %>%
  full_join(
    sample_data(ASV) %>%
      data.frame() %>%
      rownames_to_column("File") %>%
      as_tibble()
  ) %>%
  pivot_longer(cols = starts_with("ASV"),
               names_to = "ASV",
               values_to = "Reads") %>%
  full_join(
    tax_table(ASV) %>% 
      as.data.frame() %>%
      rownames_to_column("ASV") %>%
      as_tibble()
  ) %T>%
  print()

# Calculate relative abundance and presence-absence
ASV_tidy %<>%
  mutate(
    Abundance = Reads / sum(Reads),
    Presence = Reads > 0,
    .by = Sample # grouping needed for relative abundance
  ) %T>%
  print()

# 1.5.2 Sample-level summary ####
require(vegan)
ASV_summary <- ASV_tidy %>%
  summarise(
    Total = sum(Reads),
    Richness = sum(Presence),
    D = diversity(Reads, index = "invsimpson"), # inverse Simpson index
    G = diversity(Reads, index = "simpson"), # Gini-Simpson index
    E = D / Richness, # Simpson evenness
    H = diversity(Reads, index = "shannon"), # Shannon index
    J = H / log(Richness), # Pielou evenness
    .by = c(File, Number, Sample, Date, Days, Tank, Temperature, PAR, Treatment)
  ) %T>%
  print(n = 64)

# 2. Explore data ####
# 2.1 phyloseq ####
ASV %>%
  plot_richness(
    x = "Days", color = "Treatment",
    measures = c("Observed", "Shannon", "InvSimpson", "Simpson")
  )
# Some form of decline in diversity with detrital age, but not much.

ASV %>%
  plot_bar(x = "Sample", y = "Abundance", fill = "Class") +
  scale_fill_manual(
    values = c(
      "Alphaproteobacteria" = "purple",
      "Gammaproteobacteria" = "goldenrod"
    ),
    na.value = "grey" # other taxa
  )
# Mostly Alpha- and Gammaproteobacteria as expected with
# some highly abundant ASVs (boxes).

ASV %>%
  ordinate(method = "NMDS", distance = "bray") %>%
  plot_ordination(physeq = ASV, color = "Treatment")
# Bit of change in community composition with detrital age

ASV %>%
  plot_heatmap(method = "NMDS", distance = "bray")
# Some shift in abundance of certain ASVs is evident but
# not very intelligible.

# Phyloseq is not customisable enough

# 2.2 Indices ####
# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(1, "cm"),
                 text = element_text(family = "Futura"))

# Total abundance
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, Total, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(Total = mean(Total), 
                          .by = c(Treatment, Days)),
              aes(Days, Total, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Richness
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, Richness, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(Richness = mean(Richness), 
                          .by = c(Treatment, Days)),
              aes(Days, Richness, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Spikes in total abundance and richness during the
# initial deocmposition phase.

# Simpson index
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, D, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(D = mean(D), 
                          .by = c(Treatment, Days)),
              aes(Days, D, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Gini-Simpson index
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, G, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(G = mean(G), 
                          .by = c(Treatment, Days)),
              aes(Days, G, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Shannon index
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, H, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(H = mean(H), 
                          .by = c(Treatment, Days)),
              aes(Days, H, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Simpson evenness
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, E, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(E = mean(E), 
                          .by = c(Treatment, Days)),
              aes(Days, E, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Pielou evenness
ASV_summary %>%
  ggplot() +
    geom_point(aes(Days, J, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>% 
                summarise(J = mean(J), 
                          .by = c(Treatment, Days)),
              aes(Days, J, colour = Treatment)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# General decline in diversity indices with time, with the
# exception of the Light 15°C treatment which I have shown
# to have remained physiologically healthy throughout.

# 2.3 Taxa ####
# 2.3.1 Phyla ####
ASV_tidy %>%
  drop_na(Phylum) %>%
  summarise(Reads = sum(Reads),
            Presence = sum(Presence),
            ASVs = n_distinct(ASV),
            .by = Phylum) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Pseudomonadota and Bacteroidota dominate

ASV_tidy %>%
  drop_na(Phylum) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Phylum, Sample, Days, Treatment)) %>%
  mutate(Abundance_mean = mean(Abundance), .by = Phylum) %>%
  # Aggregate phyla that contribute less than 1% to samples on average
  mutate(Phylum = if_else(Abundance_mean < 0.01, "Other", Phylum)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Phylum, Sample, Days, Treatment)) %>%
  ggplot() +
    geom_point(aes(Days, Abundance, colour = Phylum), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>%
                summarise(Abundance = mean(Abundance),
                          .by = c(Phylum, Treatment, Days)),
              aes(Days, Abundance, colour = Phylum)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme
# Pseudomonadota and Bacteroidota seem to be most dynamic 

# 2.3.2 Classes ####
ASV_tidy %>%
  drop_na(Class) %>%
  summarise(Reads = sum(Reads), 
            Presence = sum(Presence),
            ASVs = n_distinct(ASV),
            .by = Class) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Gammaproteobacteria, Alphaproteobacteria and Bacteroidia dominate

ASV_tidy %>%
  drop_na(Class) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Class, Sample, Days, Treatment)) %>%
  mutate(Abundance_mean = mean(Abundance), .by = Class) %>%
  mutate(Class = if_else(Abundance_mean < 0.01, "Other", Class)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Class, Sample, Days, Treatment)) %>%
  ggplot() +
    geom_point(aes(Days, Abundance, colour = Class), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>%
                summarise(Abundance = mean(Abundance),
                          .by = c(Class, Treatment, Days)),
              aes(Days, Abundance, colour = Class)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme
# Gammaproteobacteria, Alphaproteobacteria and Bacteroidia are most dynamic

# 2.3.3 Families ####
ASV_tidy %>%
  drop_na(Family) %>%
  summarise(Reads = sum(Reads), 
            Presence = sum(Presence),
            ASVs = n_distinct(ASV),
            .by = Family) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Paracoccaceae, Hyphomonadaceae, Cellvibrionaceae etc. dominate

ASV_tidy %>%
  drop_na(Family) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Family, Sample, Days, Treatment)) %>%
  mutate(Abundance_mean = mean(Abundance), .by = Family) %>%
  mutate(Family = if_else(Abundance_mean < 0.01, "Other", Family)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Family, Sample, Days, Treatment)) %>%
  ggplot() +
    geom_point(aes(Days, Abundance, colour = Family), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>%
                summarise(Abundance = mean(Abundance),
                          .by = c(Family, Treatment, Days)),
              aes(Days, Abundance, colour = Family)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme
# Quite a few above 1%

# 2.3.4 Genera ####
ASV_tidy %>%
  drop_na(Genus) %>%
  summarise(Reads = sum(Reads), 
            Presence = sum(Presence),
            ASVs = n_distinct(ASV),
            .by = Genus) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Colwellia, Pseudahrensia, Arenicella, Hellea, Vibrio and two 
# unknown genera (Incertae Sedis_939, Incertae Sedis_761) dominate

ASV_tidy %>%
  drop_na(Genus) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Genus, Sample, Days, Treatment)) %>%
  mutate(Abundance_mean = mean(Abundance), .by = Genus) %>%
  mutate(Genus = if_else(Abundance_mean < 0.01, "Other", Genus)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Genus, Sample, Days, Treatment)) %>%
  ggplot() +
    geom_point(aes(Days, Abundance, colour = Genus), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>%
                summarise(Abundance = mean(Abundance),
                          .by = c(Genus, Treatment, Days)),
              aes(Days, Abundance, colour = Genus)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme
# Some seem to spike in the early decomposition phase

# 2.3.5 Species ####
ASV_tidy %>%
  drop_na(Species) %>%
  summarise(Reads = sum(Reads), 
            Presence = sum(Presence),
            ASVs = n_distinct(ASV),
            .by = c(Genus, Species)) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Lots of classic heterotrophs. I'll make functional groups later.

ASV_tidy %>%
  drop_na(Species) %>%
  mutate(Binomial = str_c(Genus, Species, sep = " ")) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Binomial, Sample, Days, Treatment)) %>%
  mutate(Abundance_mean = mean(Abundance), .by = Binomial) %>%
  mutate(Binomial = if_else(Abundance_mean < 0.01, "Other", Binomial)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Binomial, Sample, Days, Treatment)) %>%
  ggplot() +
    geom_point(aes(Days, Abundance, colour = Binomial), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>%
                summarise(Abundance = mean(Abundance),
                          .by = c(Binomial, Treatment, Days)),
              aes(Days, Abundance, colour = Binomial)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme


###################

# 2.4 ASVs ####
# Abundant ASVs
ASV_tidy %>%
  filter()
  arrange(desc(Abundance)) %>%
  print(n = 100)



# 3. nMDS ####
nMDS <- metaMDS()




# 4. Diversity ####


# 5. Saprotrophs ####