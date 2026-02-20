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
ASV_sample <- ASV_tidy %>%
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

# 1.6 Functional classification ####
# 1.6.1 Check taxonomic resolution ####
ASV_tidy %>%
  summarise(
    Prop_Phylum = sum(!is.na(Phylum)) / n(),
    Prop_Class = sum(!is.na(Class)) / n(),
    Prop_Family = sum(!is.na(Family)) / n(),
    Prop_Genus = sum(!is.na(Genus)) / n(),
    Prop_Species = sum(!is.na(Species)) / n()
  )
# It drops off after genus (70% classified),
# so this is the most sensible level to
# functionally classify.

# 1.6.2 Load literature review ####
# List of bacterial genera that degrade macroalgae
sapro <- read_csv("Saprotrophs.csv") %T>%
  print()

sapro %>%
  summarise(
    Studies = n_distinct(Reference),
    Mentions = n(),
    .by = Genus
  ) %>%
  arrange(desc(Studies)) %>%
  print(n = 200)
# 151 distinct genera with by far the most commonly
# studied being Vibrio, Pseudoalteromonas, Pseudomonas,
# Alteromonas, Bacillus, Cellulophaga, Zobellia,
# Halomonas, Shewanella, Cobetia, Flavobacterium,
# Colwellia, Photobacterium and Algibacter

# I'll filter by this genus list:
sapro_genus <- sapro %>% pull(Genus) %>% 
  na.omit() %>% unique() %T>% print()

# But let's look at some of the species names too.
sapro %>%
  count(Species) %>%
  arrange(desc(n)) %>%
  print(n = 300)
# Species names such as algicola, fucicola, galactanivorans, alginolyticus,
# carrageenovora, agarivorans, alginovora, fucanivorans, agarilytica,
# agarlytica, fucoidanolyticus, ulvanivorans etc. indicate their function.

# Are any of such give-away species names represented in my data?
sapro_species <- ASV_tidy %>%
  filter(
    Species %>% str_detect(
      "vora|lytic|agar|alg|fuc|lami|ulva|carra"
    )
  ) %>%
  distinct(Genus, Species) %T>%
  print(n = 40)
# Some of these obvious algae-degraders are not in any of the genera
# on the list:
sapro_species %<>%
  mutate(Represented = Genus %in% sapro_genus) %>%
  filter(!Represented) %T>%
  print()
# These species are not represented in sapro_genus, because their genus
# isn't. I'll include all of these as well, but only at species level.

sapro_species %<>%
  mutate(Binomial = str_c(Genus, Species, sep = " ")) %>%
  pull(Binomial) %T>%
  print()

# 1.6.3 Identify saprotrophs ####
# Identify genera that are likely algae saprotrophs
ASV_tidy %<>%
  mutate(
    Saprotroph = Genus %in% sapro_genus | 
      str_c(Genus, Species, sep = " ") %in% sapro_species
  ) %T>%
  print()

# 1.7 Summarise data ####
ASV_summary
###############

# 2. Explore data ####
# 2.1 phyloseq ####
ASV %>%
  plot_richness(
    x = "Days", color = "Treatment",
    measures = c("Observed", "Chao1", "Shannon", "InvSimpson", "Simpson")
  )
# Some form of decline in diversity with detrital age, but not much.

ASV %>%
  plot_bar(x = "Sample", y = "Abundance", fill = "Class") +
  scale_fill_manual(
    values = c(
      "Alphaproteobacteria" = "purple",
      "Gammaproteobacteria" = "goldenrod"
    )
  )
# Mostly Alpha- and Gammaproteobacteria as expected with
# some highly abundant ASVs (boxes). NA are other taxa.

ASV %>%
  ordinate(method = "NMDS", distance = "bray") %>%
  plot_ordination(physeq = ASV, color = "Treatment")

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
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  summarise(Total = mean(Total), 
            .by = c(Treatment, Days)) %>%
  mutate(count = if_else(Days == 0, 4, 1)) %>% uncount(count) %>%
  mutate(
    Treatment = case_when(
      row_number() == 1 ~ "Light 15°C",
      row_number() == 2 ~ "Dark 15°C",
      row_number() == 3 ~ "Light 20°C",
      row_number() == 4 ~ "Light 25°C",
      TRUE ~ Treatment
    )
  ) %>%
  ggplot() +
    geom_line(aes(Days, Total, colour = Treatment)) +
    geom_point(data = ASV_summary %>% filter(Days != 0),
               aes(Days, Total, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Richness
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  summarise(Richness = mean(Richness), 
            .by = c(Treatment, Days)) %>%
  mutate(count = if_else(Days == 0, 4, 1)) %>% uncount(count) %>%
  mutate(
    Treatment = case_when(
      row_number() == 1 ~ "Light 15°C",
      row_number() == 2 ~ "Dark 15°C",
      row_number() == 3 ~ "Light 20°C",
      row_number() == 4 ~ "Light 25°C",
      TRUE ~ Treatment
    )
  ) %>%
  ggplot() +
    geom_line(aes(Days, Richness, colour = Treatment)) +
    geom_point(data = ASV_summary %>% filter(Days != 0),
               aes(Days, Richness, colour = Treatment), 
               shape = 16, alpha = 0.5) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# Spikes in total abundance and richness during the
# initial decomposition phase.

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
# unknown genera (Incertae Sedis_939, Incertae Sedis_761) dominate.

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
# Some seem to spike in the early decomposition phase, some later.

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
# Aquimarina latercula, Sulfitobacter donghicola, 
# Shewanella surugensis dominate described species.

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
# Shewanella surugensis is only really present at the beginning
# in the Light 15°C treatment.

# 2.3.6 ASVs ####
ASV_tidy %>%
  summarise(Reads = sum(Reads), 
            Presence = sum(Presence),
            .by = c(ASV, Family, Genus, Species)) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Lots of highly abundant ASVs that are present in most samples
# (up to 61 out of 64). Most from unknown species with the most
# abundant ASVs belonging to genera Pseudahrensia, Hellea, 
# Leucothrix, Arenicella, Litoreibacter, Tateyamaria, Aquimarina etc.

ASV_tidy %>%
  summarise(Abundance = sum(Abundance),
            .by = c(ASV, Sample, Days, Treatment)) %>%
  mutate(Abundance_mean = mean(Abundance), .by = ASV) %>%
  mutate(ASV = if_else(Abundance_mean < 0.01, "Other", ASV)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(ASV, Sample, Days, Treatment)) %>%
  ggplot() +
    geom_point(aes(Days, Abundance, colour = ASV), 
               shape = 16, alpha = 0.5) +
    geom_line(data = . %>%
                summarise(Abundance = mean(Abundance),
                          .by = c(ASV, Treatment, Days)),
              aes(Days, Abundance, colour = ASV)) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# 2.4 Saprotrophs ####
ASV_tidy %>%
  filter(Saprotroph) %>%
  summarise(Reads = sum(Reads),
            Presence = sum(Presence),
            .by = c(Genus, Species)) %>%
  mutate(Proportion = Reads / sum(Reads)) %>%
  arrange(desc(Reads)) %>%
  print(n = 100)
# Colwellia, Tateyamaria, Leucothrix, Vibrio, Litoreibacter, Sulfitobacter,
# Shewanella, Pseudophaeobacter are the most abundant. Surprisingly,
# Pseudoalteromonas, Pseudomonas and Alteromonas are further down the list.

Shew_list <- ASV_tidy %>%
  filter(Genus == "Shewanella", Species == "surugensis") %>%
  distinct(ASV, Genus, Species) %>% pull(ASV)

ASV_tidy %>%
  filter(Saprotroph & !ASV %in% Shew_list) %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  summarise(Abundance = sum(Abundance),
            .by = c(Sample, Days, Treatment)) %>%
  summarise(Abundance = mean(Abundance),
            .by = c(Treatment, Days)) %>%
  mutate(count = if_else(Days == 0, 4, 1)) %>% uncount(count) %>%
  mutate(
    Treatment = case_when(
      row_number() == 1 ~ "Light 15°C",
      row_number() == 2 ~ "Dark 15°C",
      row_number() == 3 ~ "Light 20°C",
      row_number() == 4 ~ "Light 25°C",
      TRUE ~ Treatment
    )
  ) %>%
  ggplot() +
    geom_line(aes(Days, Abundance)) +
    geom_point(data = ASV_tidy %>% filter(Saprotroph & Days != 0) %>%
                 summarise(Abundance = sum(Abundance),
                           .by = c(Sample, Days, Treatment)),
               aes(Days, Abundance), 
               shape = 16, alpha = 0.5) +
    facet_grid(~ Treatment, scales = "free", space = "free") +
    mytheme

# 3. NMDS ####
# 3.1 Calculate ####
# Raw counts (reads)
NMDS_reads <- otu_table(ASV) %>%
  metaMDS(distance = "bray", autotransform = FALSE)
NMDS_reads$stress # stress = 0.1713856

# Relative abundance (normalised reads per sample)
NMDS_abund <- otu_table(ASV) %>%
  # Calculate relative abundance with phyloseq function
  transform_sample_counts(function(x) x / sum(x)) %>%
  metaMDS(distance = "bray", autotransform = FALSE)
NMDS_abund$stress # stress = 0.1685437

# 3.2 Add axes to ASV_summary ####
ASV_summary %<>%
  full_join(
    NMDS_reads$points %>%
      as.data.frame() %>%
      rownames_to_column("File") %>%
      as_tibble() %>%
      rename(MDS1_reads = MDS1, MDS2_reads = MDS2)
  ) %>%
  full_join(
    NMDS_abund$points %>%
      as.data.frame() %>%
      rownames_to_column("File") %>%
      as_tibble() %>%
      rename(MDS1_abund = MDS1, MDS2_abund = MDS2)
  ) %T>%
  print()

# 3.3 Visualise ####
# 3.3.1 Treatments ####
# Reads
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  ggplot() +
    geom_point(aes(MDS1_reads, MDS2_reads, colour = Treatment)) +
    geom_polygon(data = . %>% group_by(Treatment) %>%
                   slice(chull(MDS1_reads, MDS2_reads)),
                 aes(MDS1_reads, MDS2_reads,
                     fill = Treatment, colour = Treatment),
                 alpha = 0.2, linewidth = 0.7) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())

# Relative abundance
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  ggplot() +
    geom_point(aes(MDS1_abund, MDS2_abund, colour = Treatment)) +
    geom_polygon(data = . %>% group_by(Treatment) %>%
                   slice(chull(MDS1_abund, MDS2_abund)),
                 aes(MDS1_abund, MDS2_abund,
                     fill = Treatment, colour = Treatment),
                 alpha = 0.2, linewidth = 0.7) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
# Hardly different but relative abundance has lower stress and
# accounts for minor differences in sequencing depth, so I'll 
# probably use that.

# 3.3.2 Treatments and time ####
ordisurf(NMDS_reads, sample_data(ASV)$Days)
ordisurf(NMDS_abund, sample_data(ASV)$Days)
# Detrital age is fairly linear along NMDS1 axis

# Reads
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  left_join( # Calculate centroids
    ASV_summary %>%
      mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
      summarise(MDS1_mean = mean(MDS1_reads),
                MDS2_mean = mean(MDS2_reads),
                .by = c(Treatment, Days))
  ) %>%
  ggplot() +
    geom_point(aes(MDS1_reads, MDS2_reads, colour = Treatment),
               shape = 16, alpha = 0.3) +
    # geom_point(aes(MDS1_mean, MDS2_mean, colour = Treatment),
    #            shape = 16, size = 2.5) +
    geom_segment(aes(x = MDS1_reads, y = MDS2_reads, 
                     xend = MDS1_mean, yend = MDS2_mean,
                     colour = Treatment), alpha = 0.3) +
    geom_path(data = . %>% distinct(Days, Treatment, MDS1_mean, MDS2_mean) %>% 
                mutate(count = if_else(Days == 0, 4, 1)) %>% uncount(count) %>%
                mutate(
                  Treatment = case_when(
                    row_number() == 1 ~ "Light 15°C",
                    row_number() == 2 ~ "Dark 15°C",
                    row_number() == 3 ~ "Light 20°C",
                    row_number() == 4 ~ "Light 25°C",
                    TRUE ~ Treatment
                  )
                ),
              aes(MDS1_mean, MDS2_mean, colour = Treatment),
              arrow = arrow(length = unit(0.3, "cm"), 
                            type = "closed", angle = 20)) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())

# Relative abundance
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  left_join(
    ASV_summary %>%
      mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
      summarise(MDS1_mean = mean(MDS1_abund),
                MDS2_mean = mean(MDS2_abund),
                .by = c(Treatment, Days))
  ) %>%
  ggplot() +
    geom_point(aes(MDS1_abund, MDS2_abund, colour = Treatment),
               shape = 16, alpha = 0.3) +
    # geom_point(aes(MDS1_mean, MDS2_mean, colour = Treatment),
    #            shape = 16, size = 2.5) +
    geom_segment(aes(x = MDS1_abund, y = MDS2_abund, 
                     xend = MDS1_mean, yend = MDS2_mean,
                     colour = Treatment), alpha = 0.3) +
    geom_path(data = . %>% distinct(Days, Treatment, MDS1_mean, MDS2_mean) %>% 
                mutate(count = if_else(Days == 0, 4, 1)) %>% uncount(count) %>%
                mutate(
                  Treatment = case_when(
                    row_number() == 1 ~ "Light 15°C",
                    row_number() == 2 ~ "Dark 15°C",
                    row_number() == 3 ~ "Light 20°C",
                    row_number() == 4 ~ "Light 25°C",
                    TRUE ~ Treatment
                  )
                ),
              aes(MDS1_mean, MDS2_mean, colour = Treatment),
              arrow = arrow(length = unit(0.3, "cm"), 
                            type = "closed", angle = 20)) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())

# 4. dbRDA ####
# 4.1 Calculate ####
# Reads
dbRDS_reads <- dbrda(
  otu_table(ASV) ~ Days + Temperature + PAR,
  data = sample_data(ASV) %>% data.frame(), 
  distance = "bray"
)
dbRDS_reads

# Relative abundance
dbRDS_abund <- dbrda(
  otu_table(ASV) %>% 
    transform_sample_counts(function(x) x / sum(x)) ~ 
    Days + Temperature + PAR,
  data = sample_data(ASV) %>% data.frame(), 
  distance = "bray"
)
dbRDS_abund
# Results are again nearly identical

# 4.2 Extract data ####
# Add choices = 1:3 for all three axes
scores_reads <- scores(dbRDS_reads, tidy = TRUE) %T>%
  print()
scores_abund <- scores(dbRDS_abund, tidy = TRUE) %T>%
  print()

# Save effects
effects <- bind_rows(
  Reads = scores_reads %>%
    filter(score == "biplot"),
  Abundance = scores_abund %>%
    filter(score == "biplot"),
  .id = "Response"
) %>%
  select(-score) %>%
  as_tibble() %>%
  mutate(
    label = label %>% replace_values(
      "PAR" ~ "Light",
      "Days" ~ "Age"
    )
  ) %T>%
  print()

# 4.3 Add axes to ASV_summary ####
ASV_summary %<>%
  full_join(
    scores_reads %>%
      filter(score == "sites") %>%
      select(-score) %>%
      as_tibble() %>%
      rename(File = label, 
             dbRDA1_reads = dbRDA1, 
             dbRDA2_reads = dbRDA2)
  ) %>%
  full_join(
    scores_abund %>%
      filter(score == "sites") %>%
      select(-score) %>%
      as_tibble() %>%
      rename(File = label, 
             dbRDA1_abund = dbRDA1, 
             dbRDA2_abund = dbRDA2)
  ) %T>%
  print()

# 4.4 Visualise ####
# Reads
require(geomtextpath)
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  ggplot() +
    geom_point(aes(dbRDA1_reads, dbRDA2_reads, colour = Treatment)) +
    geom_polygon(data = . %>% group_by(Treatment) %>%
                   slice(chull(dbRDA1_reads, dbRDA2_reads)),
                 aes(dbRDA1_reads, dbRDA2_reads,
                     fill = Treatment, colour = Treatment),
                 alpha = 0.2, linewidth = 0.7) +
    geom_textsegment(data = effects %>% filter(Response == "Reads"),
              aes(x = 1.5, y = 2, xend = dbRDA1+1.5, yend = dbRDA2+2, 
                  label = label), hjust = 1, family = "Futura", size = 3.5,
              arrow = arrow(length = unit(0.3, "cm"), 
                            type = "closed", angle = 20)) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())

# Relative abundance
ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  ggplot() + # dbRDA randomly flips axis signs, so I am reversing that
    geom_point(aes(dbRDA1_abund, dbRDA2_abund*-1, colour = Treatment)) +
    geom_polygon(data = . %>% group_by(Treatment) %>%
                   slice(chull(dbRDA1_abund, dbRDA2_abund)),
                 aes(dbRDA1_abund, dbRDA2_abund*-1,
                     fill = Treatment, colour = Treatment),
                 alpha = 0.2, linewidth = 0.7) +
    geom_textsegment(data = effects %>% filter(Response == "Abundance"),
              aes(x = 1.5, y = 2, xend = dbRDA1+1.5, yend = dbRDA2*-1+2,
                  label = label), hjust = 1, family = "Futura", size = 3.5,
              arrow = arrow(length = unit(0.3, "cm"),
                            type = "closed", angle = 20)) +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank())
# Again, hardly different.

# 5. SIMPER ####
# 5.1 Calculate ####
# Treatment is the only sensible categorical variable
SIMPER_reads <- simper(
  otu_table(ASV), # First five samples are baseline
  sample_data(ASV)$Treatment %>% replace(1:5, "Baseline")
) %T>%
  print()

summary(SIMPER_reads)

SIMPER_abund <- simper(
  otu_table(ASV) %>% transform_sample_counts(function(x) x / sum(x)),
  sample_data(ASV)$Treatment %>% replace(1:5, "Baseline")
) %T>%
  print()

summary(SIMPER_abund)
# Both are similar. Proceed with relative abundance version.

# 5.2 Filter ####
# cusum gives the cumulative proportion of the difference
# explained by that ASV, so I can filter by a threshold.

# Tidy up
SIMPER_tidy <- SIMPER_abund %>%
  map(as_tibble) %>%
  list_rbind(names_to = "Contrast") %>%
  mutate(Difference = abs(avb - ava)) %T>%
  print()

# Filter the fewest ASVs that explain 50% of each difference
SIMPER_tidy %<>%
  filter(cusum <= 0.5) %T>%
  print()

# 5.3 Identify ####
SIMPER_tidy %>%
  summarise(Difference_mean = mean(Difference),
            Difference_sd = sd(Difference),
            Contrasts = n_distinct(Contrast),
            .by = species) %>%
  rename(ASV = species) %>%
  left_join(
    tax_table(ASV) %>% 
      as.data.frame() %>%
      rownames_to_column("ASV") %>%
      as_tibble()
  ) %>%
  select(-c(Kingdom, Phylum)) %>%
  arrange(desc(Difference_mean)) %>%
  print(n = 50)
# These are also the most abundant ASVs, e.g. Pseudahrensia,
# Shewanella surugensis, Colwellia, Hellea, Aquimarina latercula etc.

# 6. Figure 3 ####
# 6.1 Figure 3a ####
# This is the multidimensional community structure figure.
# I am using relative abundance in all cases.
Fig_3a <- ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  ggplot() +
  geom_point(aes(MDS1_abund, MDS2_abund, colour = Treatment),
             shape = 16) +
  geom_polygon(data = . %>% group_by(Treatment) %>%
                 slice(chull(MDS1_abund, MDS2_abund)),
               aes(MDS1_abund, MDS2_abund,
                   fill = Treatment, colour = Treatment),
               alpha = 0.2, linejoin = "round") +
  scale_colour_manual(values = c("#C3B300", "#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c")) +
  scale_fill_manual(values = c("#C3B300", "#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c")) +
  coord_cartesian(clip = "off") +
  mytheme +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(linejoin = "mitre"))
Fig_3a

# 6.2 Figure 3b ####
Fig_3b <- ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  left_join(
    ASV_summary %>%
      mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
      summarise(MDS1_mean = mean(MDS1_abund),
                MDS2_mean = mean(MDS2_abund),
                .by = c(Treatment, Days))
  ) %>%
  ggplot() +
    geom_point(aes(MDS1_abund, MDS2_abund, colour = Treatment),
               shape = 16, alpha = 0.2) +
    geom_segment(aes(x = MDS1_abund, y = MDS2_abund, 
                     xend = MDS1_mean, yend = MDS2_mean,
                     colour = Treatment), alpha = 0.2,
                 lineend = "round") +
    geom_path(data = . %>% distinct(Days, Treatment, MDS1_mean, MDS2_mean) %>% 
                mutate(count = if_else(Days == 0, 4, 1)) %>% uncount(count) %>%
                mutate(
                  Treatment = case_when(
                    row_number() == 1 ~ "Light 15°C",
                    row_number() == 2 ~ "Dark 15°C",
                    row_number() == 3 ~ "Light 20°C",
                    row_number() == 4 ~ "Light 25°C",
                    TRUE ~ Treatment
                  )
                ),
              aes(MDS1_mean, MDS2_mean, colour = Treatment),
              arrow = arrow(length = unit(0.2, "cm"), 
                            type = "closed", angle = 20),
              lineend = "round", linejoin = "round") +
    scale_colour_manual(values = c("#C3B300", "#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c")) +
    scale_fill_manual(values = c("#C3B300", "#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c")) +
    coord_cartesian(clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(linejoin = "mitre"))

Fig_3b

# 6.3 Figure 3c ####
Fig_3c <- ASV_summary %>%
  mutate(Treatment = if_else(Days == 0, "Baseline", Treatment)) %>%
  ggplot() + # dbRDA randomly flips axis signs, so I am reversing that
    geom_point(aes(dbRDA1_abund, dbRDA2_abund*-1, colour = Treatment),
               shape = 16) +
    geom_polygon(data = . %>% group_by(Treatment) %>%
                   slice(chull(dbRDA1_abund, dbRDA2_abund)),
                 aes(dbRDA1_abund, dbRDA2_abund*-1,
                     fill = Treatment, colour = Treatment),
                 alpha = 0.2, linejoin = "round") +
    # geom_textsegment(data = effects %>% filter(Response == "Abundance"),
    #           aes(x = 1, y = 2, xend = dbRDA1+1, yend = dbRDA2*-1+2,
    #               label = label), hjust = 1, family = "Futura", size = 3.5, # ~10 pt
    #           arrow = arrow(length = unit(0.3, "cm"),
    #                         type = "closed", angle = 20)) +
    geom_segment(data = effects %>% filter(Response == "Abundance"),
                 aes(x = 0.9, y = 1.55, xend = dbRDA1+0.9, yend = dbRDA2*-1+1.55),
                 arrow = arrow(length = unit(0.2, "cm"),
                               type = "closed", angle = 20),
                 lineend = "round", linejoin = "round") +
    geom_text(data = effects %>% filter(Response == "Abundance"),
              aes(dbRDA1+0.9, dbRDA2*-1+1.55, label = label),
              hjust = c(0.4, 0.5, 0), vjust = c(1.3, -0.3, -0.3),
              family = "Futura", size.unit = "pt", size = 10) +
    scale_colour_manual(values = c("#C3B300", "#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c")) +
    scale_fill_manual(values = c("#C3B300", "#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c")) +
    coord_cartesian(clip = "off") +
    mytheme +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(linejoin = "mitre"))
Fig_3c

# 6.4 Combined panels ####
require(patchwork) # Fig_3b has different colour guides, which I remove
Fig_3 <- ( Fig_3a | Fig_3b + guides(colour = "none") | Fig_3c ) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = c("a", "b")) &
  theme(legend.position = "top", legend.justification = 0,
        plot.tag = element_text(family = "Futura", size = 12, face = "bold"),
        plot.tag.position = c(0.05, 0.83))

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", path = "Figures",
         device = cairo_pdf, width = 20, height = 8, units = "cm")

# 7. Figure 4 ####
# 7.1 Figure 4a ####
Fig_4a # Total abundance

# 7.2 Figure 4b ####
Fig_4b # Richness

# 7.3 Figure 4d ####
Fig_4c # Classes > 1% relative abundance

# 7.4 Figure 4d ####
Fig_4d # Relative abundance of algal hetero-/saprotrophs