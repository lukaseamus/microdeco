# 1. Check example sequences ####
require(tidyverse)
require(magrittr)
require(here)
require(Biostrings)

# Forward example
example_forward <- here("DNA", "Sequences", 
     "J173-16S_V3-V4_L2H8N_TCGTTGCTGC-GTCTCGTGAA_L001_R1.fastq.gz") %>%
  readDNAStringSet(format = "fastq") %T>%
  print()

# Reverse example
example_reverse <- here("DNA", "Sequences", 
     "J173-16S_V3-V4_L2H8N_TCGTTGCTGC-GTCTCGTGAA_L001_R2.fastq.gz") %>%
  readDNAStringSet(format = "fastq") %T>%
  print()

# 2. Load all sequences ####
require(dada2)
sequences <- here("DNA", "Sequences") %>%
  list.files(pattern = "\\.fastq.gz$", full.names = TRUE) %>%
  tibble(Path = .) %>%
  mutate(
    Number = Path %>% basename() %>% 
      str_split_i("-", 1) %>% 
      str_extract("\\d+") %>%
      as.numeric(),
    Direction = if_else(
      Path %>% basename() %>% str_detect("R1"),
      "Forward", "Reverse"
    ) %>% fct()
  ) %T>%
  print(n = 140)

# 3. Add sample names ####
samples <- here("DNA", "Samples.csv") %>% 
  read_csv() %T>%
  print()

sequences %<>%
  left_join(samples) %T>%
  print()

# 4. Check for primers ####
primer_341F <- DNAString("CCTAYGGGRBGCASCAG") %T>% print()
primer_806R <- DNAString("GGACTACNNGGGTATCTAAT") %T>% print()

# Check forward sequence with forward primer
example_forward %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_341F, fixed = FALSE) 
# Only one match that starts at the beginning of the DNA sequence.
# Aligned this looks like this:
example_forward %>%
  sample(1) %>% .[[1]] %>%
  pwalign::pairwiseAlignment(pattern = primer_341F)

# Check forward sequence with reverse primer
example_forward %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_806R, fixed = FALSE) 
# As expected, no match.

# Check reverse sequence with reverse primer
example_reverse %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_806R, fixed = FALSE) 
# Only one match that starts at the beginning of the DNA sequence.
example_reverse %>%
  sample(1) %>% .[[1]] %>%
  pwalign::pairwiseAlignment(pattern = primer_806R)

# Check reverse sequence with forward primer
example_reverse %>%
  sample(1) %>% .[[1]] %>%
  matchPattern(pattern = primer_341F, fixed = FALSE) 
# No match.

# So I can simply remove the first 17 bases from the forward sequence
# and the first 20 bases from the reverse sequence to remove primers.

# 5. Quality control ####
require(glue)
sequences %<>%
  rowwise() %>% # rowwise makes plotting from many variables easier
  mutate( 
    # Here I include the sample name and read direction along with
    # a conditional guessed trim for reference.
    Quality_Plot = 
      list(
        Path %>%
          plotQualityProfile() +
            geom_hline(yintercept = 30) + # 99.9% accuracy
            # My guesses for quality drop below 30 are 260 for 
            # forward and 200 for reverse.
            geom_vline(
              xintercept = if (Direction == "Forward") 260 else 200
            ) +
            ggtitle(glue("{Sample} ({Direction})")) 
            # Alternatively, ggtitle(Sample, Direction) creates 
            # title and subtitle.
      )
  ) %>%
  ungroup() %T>% # rowwise needs to be undone
  print()

require(patchwork)
sequences %$%
  wrap_plots(Quality_Plot) %>%
  ggsave(filename = "quality.pdf", path = here("DNA", "Plots"),
         device = cairo_pdf, height = 100, width = 100, units = "cm")
# Trim positions 290 and 220 generally look good.

# 6. Filter and trim ####
# Create paths for filtered sequences
sequences %<>%
  pivot_wider(names_from = Direction, values_from = c(Path, Quality_Plot)) %>%
  mutate(
    Path_Forward_Filtered = here("DNA", "Sequences", "Filtered") %>%
      file.path(str_c(Sample, "_F_filt.fastq.gz")),
    Path_Reverse_Filtered = here("DNA", "Sequences", "Filtered") %>%
      file.path(str_c(Sample, "_R_filt.fastq.gz"))
  ) %T>%
  print()

# My ideal trimming lengths are 260 and 200 bp for quality control, but 
# I need to ensure enough overlap for merging. I tried 260 and 200 and
# this lead to a large drop in read length at the merging step. 

# My amplicon length is
length(341:806) # 466 bp
# DADA2 needs 20 bp overlap
260 + 200 - 466 # My planned truncation allows no overlap, so I'll go
# with 280 and 210 bp instead:
280 + 210 - 466 # 24 bp overlap

filter <- sequences %$%
  filterAndTrim( # The paths need to be named to be matched correctly.
    fwd = Path_Forward %>% setNames(Sample), filt = Path_Forward_Filtered %>% setNames(Sample),
    rev = Path_Reverse %>% setNames(Sample), filt.rev = Path_Reverse_Filtered %>% setNames(Sample), 
    truncLen = c(280, 210), # Trim forward to 280 bases and reverse to 210
    trimLeft = c(length(primer_341F), length(primer_806R)), # Remove primers
    maxN = 0, maxEE = c(2, 2), truncQ = 2, 
    rm.phix = TRUE, compress = TRUE, multithread = TRUE
  )

# Watch files appear in ~/DNA/Sequences/Filtered to track progress

sequences %<>%
  mutate(Reads = filter[,1],
         Reads_Filtered = filter[,2],
         Prop_Filtered = Reads_Filtered / Reads) %T>%
  print()

rm(filter)

sequences %>%
  summarise( across(Prop_Filtered, list( mean = mean, sd = sd )) )
# 68 ± 4.2 % of reads survived filtering.

( sequences %>%
  ggplot(aes(Reads, Reads_Filtered)) +
      geom_abline(slope = c(1, sequences %$% mean(Prop_Filtered)),
                  colour = c("black", "grey")) +
      geom_label(aes(label = Sample),) +
      theme_minimal() ) %>%
  ggsave(filename = "filtering.pdf", path = here("DNA", "Plots"),
         device = cairo_pdf, height = 30, width = 30, units = "cm")

# 7. Post-filter quality control ####
sequences %<>%
  rowwise() %>%
  mutate(
    Quality_Plot_Forward_Filtered = 
      list(
        Path_Forward_Filtered %>%
          plotQualityProfile() +
            ggtitle(glue("{Sample} (Forward)"))
      ),
    Quality_Plot_Reverse_Filtered = 
      list(
        Path_Reverse_Filtered %>%
          plotQualityProfile() +
            ggtitle(glue("{Sample} (Reverse)"))
      )
  ) %>%
  ungroup()

sequences %$%
  wrap_plots(c(Quality_Plot_Forward_Filtered, 
               Quality_Plot_Reverse_Filtered)) %>%
  ggsave(filename = "quality_filtered.pdf", path = here("DNA", "Plots"),
         device = cairo_pdf, height = 100, width = 100, units = "cm")
# Looks good.

# 8. Learn error rates ####
# Samples are ordered by treatment to randomisation is important
error_forward <- sequences %$% 
  learnErrors(Path_Forward_Filtered,
              multithread = TRUE,
              randomize = TRUE) 

error_reverse <- sequences %$% 
  learnErrors(Path_Reverse_Filtered,
              multithread = TRUE,
              randomize = TRUE)

# Plot
( plotErrors(error_forward, nominalQ = TRUE) +
  ggtitle("Forward") | 
  plotErrors(error_reverse, nominalQ = TRUE) +
  ggtitle("Reverse") ) %>%
  ggsave(filename = "errors.pdf", path = here("DNA", "Plots"),
         device = cairo_pdf, height = 30, width = 60, units = "cm")

# 9. Denoise and infer ASVs ####
sequences %<>%
  mutate(
    DADA_Forward = dada(
      derep = Path_Forward_Filtered,
      err = error_forward,
      pool = "pseudo",
      multithread = TRUE
    ),
    DADA_Reverse = dada(
      derep = Path_Reverse_Filtered,
      err = error_reverse,
      pool = "pseudo",
      multithread = TRUE
    )
  ) %T>%
  print()
  
# 10. Merge forward and reverse ####
sequences %<>%
  mutate(
    DADA_Merged = mergePairs(
      dadaF = DADA_Forward, derepF = Path_Forward_Filtered,
      dadaR = DADA_Reverse, derepR = Path_Reverse_Filtered,
      verbose = TRUE
    )
  ) %T>%
  print()

# 11. Sequence table ####
# Create naive sequence table
sequence_table <- sequences %$% 
  makeSequenceTable(DADA_Merged)
dim(sequence_table) # 65 samples with 48964 ASVs

( sequence_table %>%
  getSequences() %>%
  str_length() %>%
  tibble(Basepairs = .) %>%
  count(Basepairs, name = "ASVs") %T>%
  print(n = 60) %>%
  ggplot(aes(Basepairs, ASVs)) +
    geom_point() +
    geom_line() +
    scale_y_log10() +
    theme_minimal() ) %>%
  ggsave(filename = "ASV_length_naive.pdf", path = here("DNA", "Plots"),
         device = cairo_pdf, height = 10, width = 10, units = "cm")
# Most ASVs are 404 and 429 bp long

sequence_table %>%
  as_tibble() %>%
  print(n = 65)

# Remove chimeras
sequence_table_clean <- removeBimeraDenovo(
  sequence_table, method = "consensus", 
  multithread = TRUE, verbose = TRUE
)
# Identified 42552 chimeras out of 48964 ASVs.

( sequence_table_clean %>%
  getSequences() %>%
  str_length() %>%
  tibble(Basepairs = .) %>%
  count(Basepairs, name = "ASVs") %T>%
  print(n = 60) %>%
  ggplot(aes(Basepairs, ASVs)) +
    geom_point() +
    geom_line() +
    scale_y_log10() +
    theme_minimal() ) %>%
  ggsave(filename = "ASV_length_clean.pdf", path = here("DNA", "Plots"),
         device = cairo_pdf, height = 10, width = 10, units = "cm")
# Most ASVs are still 404 and 429 bp long but the spikes are less extreme.
# There is a strange tail of a few short ASVs less than 380 bp long, which
# I'll remove:
sequence_table_clean <- 
  sequence_table_clean[ , nchar(colnames(sequence_table_clean)) > 380 ]

sequence_table_clean %>%
  as_tibble() %>%
  print(n = 65)

# Calculate proportion of ASVs lost during cleaning
1 - ncol(sequence_table_clean) / ncol(sequence_table)
# 87% of inferred ASVs are chimeras or too short
1 - sum(sequence_table_clean) / sum(sequence_table)
# But removed ASVs only account for 18% of merged reads

# 12. Track reads and ASVs ####
# Count reads at various stages after filtering
sequences %<>%
  mutate(Reads_Forward_Denoised = DADA_Forward %>% 
           map_dbl(~ .x %>% getUniques() %>% sum()),
         Reads_Reverse_Denoised = DADA_Reverse %>% 
           map_dbl(~ .x %>% getUniques() %>% sum()),
         Reads_Merged = DADA_Merged %>% 
           map_dbl(~ .x %>% getUniques() %>% sum()),
         Reads_Sequence_Table = sequence_table %>%
           rowSums(),
         Reads_Sequence_Table_Clean = sequence_table_clean %>%
           rowSums()) %T>%
  print()

sequences %>%
  select(Sample, starts_with("Reads")) %>%
  print(n = 65)

sequences %<>%
  mutate(Prop_Chimera = 1 - Reads_Sequence_Table_Clean / Reads_Sequence_Table,
         Prop_Retained = Reads_Sequence_Table_Clean / Reads) %T>%
  print()

sequences %>%
  select(Sample, starts_with("Prop")) %>%
  print(n = 65)

# Count ASVs per sample before and after chimera removal
sequences %<>%
  mutate(ASVs = rowSums(sequence_table > 0),
         ASVs_Clean = rowSums(sequence_table_clean > 0),
         Prop_ASVs = ASVs_Clean / ASVs) %T>%
  print()

sequences %>%
  select(Sample, contains("ASV")) %>%
  print(n = 65)

sequences %>%
  summarise( across(starts_with("Prop"), list( mean = mean, sd = sd )) )
# 18 ± 12 % of reads from chimeras and short ASVs
# 52 ± 11 % of initial reads retained
# 32 ± 15 % of initial ASVs retained

sequences %>%
  summarise( 
    across(
      c(Reads_Sequence_Table_Clean, ASVs_Clean), 
      list( mean = mean, sd = sd )
    ) 
  )
# 40046 ± 12919 final reads per sample
# 466 ± 176 final ASVs per sample

# Save sequence table
sequence_table_clean %>%
  write_rds(here("DNA", "RDS", "sequence_table.rds"))

# 13. Assign taxonomy ####
# 13.1 Genus-level naive Bayesian classification ####
taxa_table <- assignTaxonomy(
  seqs = sequence_table_clean,
  refFasta = here("DNA", "SILVA", "silva_nr99_v138.2_toGenus_trainset.fa.gz"),
  multithread = TRUE
)

taxa_table %>% as_tibble()

# 13.2 Species exact matching ####
# Add species to genus-level table
taxa_table %<>%
  addSpecies(refFasta = here("DNA", "SILVA", "silva_v138.2_assignSpecies.fa.gz"),
             allowMultiple = TRUE)

taxa_table %>% as_tibble()

# Assign species directly, skipping genus
species_table <- assignSpecies(
  seqs = sequence_table_clean,
  refFasta = here("DNA", "SILVA", "silva_v138.2_assignSpecies.fa.gz"),
  allowMultiple = TRUE
)

species_table %>% as_tibble()

# Compare to taxa table to species table
taxa_table %>%
  as_tibble() %>%
  select(Genus, Species) %>%
  bind_cols(species_table %>% 
              as_tibble()) %>%
  filter(Genus...1 != Genus...3 | Species...2 != Species...4 )

# Clearly addSpecies() is not independent because mistakes made by
# assignTaxonomy() at the genus level are carried through. It is now
# possible to run naive Bayesian classification to species level.

# 13.3 Species-level naive Bayesian classification ####
taxa_table_species <- assignTaxonomy(
  seqs = sequence_table_clean,
  refFasta = here("DNA", "SILVA", "silva_nr99_v138.2_toSpecies_trainset.fa.gz"),
  multithread = TRUE
)

taxa_table_species %>% as_tibble()

# Compare to taxa table to previous taxa table
taxa_table_species %>%
  as_tibble() %>%
  select(Genus, Species) %>%
  bind_cols(taxa_table %>% 
              as_tibble() %>%
              select(Genus, Species)) %>%
  filter(Genus...1 != Genus...3 | Species...2 != Species...4 ) %>%
  print(n = 200)
# 156 don't match

# Compare to taxa table to species table
taxa_table_species %>%
  as_tibble() %>%
  select(Genus, Species) %>%
  bind_cols(species_table %>% 
              as_tibble()) %>%
  filter(Genus...1 != Genus...3 | Species...2 != Species...4 ) %>%
  print(n = 200)
# 71 don't match

# 13.4 Resolve species-level classification ####
# The best native DADA2 taxonomy assignment seems to be to run the 
# naive Bayesian classifier to species level and independently run 
# exact species matching. The two separate estimates can then be 
# cross-validated to get the most complete taxa table. When only
# one approach yields species, that species is accepted. Similarly,
# when both approaches agree, there is no conflict. However,
# when both approaches yield different species, the choice is made 
# according to these rules:
# 1. If genus and/or species clash, exact matching is trusted
#    over naive Bayesian classification.
# 2. If genus aligns and exact matching identifies several species
#    including the one identified by naive Bayesian classification
#    the single species identified by the latter method is trusted.
# 3. If none of the species identified by exact matching match the
#    species identified by naive Bayesian classification, the list 
#    of species is favoured over the single species.

taxa_table_resolved <- 
  taxa_table_species %>%
  as_tibble() %>%
  bind_cols(
    species_table %>%
      as_tibble() %>% # rename doesn't work here
      mutate(Genus_exact = Genus, Species_exact = Species) %>%
      select(Genus_exact, Species_exact)
  ) %>%
  mutate(
    Genus_resolved =
      case_when(
        !is.na(Genus) & is.na(Genus_exact) ~ Genus,
        is.na(Genus) & !is.na(Genus_exact) ~ Genus_exact,
        Genus == Genus_exact ~ Genus,
        Genus != Genus_exact ~ Genus_exact
      ),
    Species_resolved =
      case_when(
        !is.na(Species) & is.na(Species_exact) ~ Species,
        is.na(Species) & !is.na(Species_exact) ~ Species_exact,
        Species_exact %>% str_detect(fixed(Species)) ~ Species,
        !Species_exact %>% str_detect("/") | 
          !Species_exact %>% str_detect(fixed(Species)) ~ Species_exact
      )
  ) %T>%
  print()

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact, 
         Genus_resolved, Species_resolved) %>%
  filter(!if_all(contains("Species"), is.na)) %>%
  print(n = 1e3)
# Looks fine but better double-check specific cases.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus == Genus_exact & Species == Species_exact) %>%
  print(n = 100)
# 55 exact genus and species matches.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus == Genus_exact & Species != Species_exact) %>%
  mutate(Corrected = Species != Species_resolved) %>%
  print(n = 100)
# 65 cases where genus matches but species clashes. Most cases
# have several options from exact matching that include the naive
# Bayesian classifier species. In 6 cases there is a species 
# mismatch despite genus match and the exact match is favoured. 
# Everything is resolved properly.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(Genus != Genus_exact) %>%
  print(n = 10)
# 6 cases where genus is mismatched. In the two cases where
# species is also assigned it's just spelled differently.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(is.na(Species) & !is.na(Species_exact)) %>%
  print(n = 100)
# 74 cases where naive Bayesian classifier didn't get the species
# but exact matching did.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(!is.na(Species) & is.na(Species_exact)) %>%
  print(n = 500)
# 431 cases where naive Bayesian classifier got the species but
# exact matching got neither genus nor species.

taxa_table_resolved %>%
  select(Genus, Species, Genus_exact, Species_exact,
         Genus_resolved, Species_resolved) %>%
  filter(is.na(Genus) & !is.na(Species_exact)) %>%
  print(n = 100)
# Only 1 case where naive Bayesian classifier got neither genus nor
# species but exact matching did.

# 13.5 DECIPHER classification algorithm ####
require(DECIPHER)
# The pipe doesn't work with load().
load(here("DNA", "SILVA", "SILVA_SSU_r138.2.RData"))

taxa_table_decipher <- IdTaxa(
  test = sequence_table_clean %>%
    getSequences() %>%
    DNAStringSet(),
  trainingSet = trainingSet,
  strand = "top", processors = NULL
)

# Reformat to match taxa_table.
taxa_table_decipher %<>%
  map(
    ~ {
      m <- match(c("domain", "phylum", "class", "order", 
                   "family", "genus", "species"), .x$rank)
      taxa <- .x$taxon[m]
      taxa[startsWith(taxa, "unclassified_")] <- NA
      taxa
    }
  ) %>% 
  simplify2array() %>%
  t() %>%
  `colnames<-`( c("domain", "phylum", "class", "order", 
                  "family", "genus", "species") ) %>%
  `rownames<-`( getSequences(sequence_table_clean) )

taxa_table_decipher %>% as_tibble()
# Generally DECIPHER seems to be worse at classifying species.
# Let's have a closer look.

taxa_table_resolved %>%
  select(Genus_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus)
  ) %>%
  filter(Genus_resolved != genus)
# Genus differences seem to be mostly due to spelling.

taxa_table_resolved %>%
  select(Genus_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus)
  ) %>%
  filter(Genus_resolved != genus & 
           !genus %>% str_detect(fixed(Genus_resolved))) %>%
  print(n = 50)
# Only 47 genera (4.1% of mismatches) are actually differently 
# classified and DECIPHER mostly assigns Incertae Sedis, 
# so placeholders.

taxa_table_decipher %>%
  as_tibble() %>%
  drop_na(species)
# DECIPHER didn't get a single species.

taxa_table_resolved %>%
  select(Genus_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus)
  ) %>%
  filter(is.na(Genus_resolved) & !is.na(genus))
# 766 times DECIPHER got the genus where the previous method didn't,
# but most are Incertae Sedis.

taxa_table_resolved %>%
  select(Genus_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus)
  ) %>%
  filter(!is.na(Genus_resolved) & is.na(genus))
# 1811 times the previous method got the genus but DECIPHER didn't.

taxa_table_resolved %>%
  select(Genus_resolved) %>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>%
      select(genus)
  ) %>%
  filter(genus %>% str_detect(fixed(Genus_resolved)))
# Genus matched or contained 1844 times.

# 13.6 Resolve all methods ####
# The only part that would improve my current taxa table would be 
# adding the 264 missing genera, but only if the other taxa match
# or are absent in the older method.
taxa_table_resolved %<>%
  bind_cols(
    taxa_table_decipher %>%
      as_tibble() %>% # rename also doesn't work here
      mutate(Phylum_decipher = phylum, Class_decipher = class, 
             Order_decipher = order, Family_decipher = family, 
             Genus_decipher = genus) %>%
      select(Phylum_decipher, Class_decipher, Order_decipher, 
             Family_decipher, Genus_decipher)
  ) %>%
  mutate(
    Genus_resolved_decipher =
      if_else( # Only add for missing genera if all else matches or is missing
        is.na(Genus_resolved) &
          ( Phylum_decipher %>%
              str_detect(Phylum %>% fixed()) |
              is.na(Phylum) ) &
          ( Class_decipher %>%
                str_detect(Class %>% fixed()) |
              is.na(Class) ) &
          ( Order_decipher %>%
                str_detect(Order %>% fixed()) |
              is.na(Order) ) &
          ( Family_decipher %>%
                str_detect(Family %>% fixed()) |
              is.na(Family) ),
        Genus_decipher, Genus_resolved
      )
  ) %T>%
  print()

# Filter cases where DECIPHER resolved a missing genus, but the 
# remaining taxonomy doesn't match, leading to rejection
taxa_table_resolved %>%
  filter(
    is.na(Genus_resolved_decipher) & !is.na(Genus_decipher) &
      (
        !Family_decipher %>% str_detect(Family %>% fixed()) | is.na(Family) |
          !Order_decipher %>% str_detect(Order %>% fixed()) | is.na(Order) |
          !Class_decipher %>% str_detect(Class %>% fixed()) | is.na(Class) |
          !Phylum_decipher %>% str_detect(Phylum %>% fixed()) | is.na(Phylum)
      )
  ) %>%
  select(Phylum, Class, Order, Family, Phylum_decipher, 
         Class_decipher, Order_decipher, Family_decipher)
# Only 14 cases. Everything matches up to class, but in
# order there are a few mismatches and more in family
# where most are Incertae Sedis for decipher, so I trust
# the previous method on these.
  
# Clean taxa_table_resolved  
taxa_table_resolved %<>%
  select(Kingdom, Phylum, Class, Order, Family, 
         Genus_resolved_decipher, Species_resolved) %>%
  mutate(Genus = Genus_resolved_decipher, Species = Species_resolved) %>%
  select(-c(Genus_resolved_decipher, Species_resolved)) %T>%
  print()

# Calculate final resolution
taxa_table_resolved %>%
  summarise(Prop_Species = sum(!is.na(Species)) / n(),
            Prop_Genus = sum(!is.na(Genus)) / n())
# 9.8% of ASVs are classified to species.
# 70% of ASVs are classified to genus.

# Convert back to matrix
taxa_table %>% str()
taxa_table_species %>% str()
taxa_table_decipher %>% str()
# This is the matrix structure I need.

taxa_table_resolved %<>%
  as.matrix() %>%
  `rownames<-`( rownames(taxa_table_species) ) %T>%
  # Alternatively, get sequences from sequence_table_clean.
  str()
# Looks fine.

# Save taxa table
taxa_table_resolved %>%
  write_rds(here("DNA", "RDS", "taxa_table.rds"))

# 14. Save ASVs with taxonomy as phyloseq object ####
# Add metadata
treatments <- here("DNA", "Treatments.csv") %>%
  read_csv() %>%
  mutate(Date = Date %>% dmy(),
         Days = Date[1] %--% Date / ddays()) %T>%
  print(n = 85)

sequences %<>%
  left_join(treatments) %T>%
  print(n = 65)

# Build phyloseq object
require(phyloseq)
ASV <- phyloseq(
  otu_table(
    sequence_table_clean, 
    taxa_are_rows = FALSE
  ),
  sample_data(
    sequences %>% 
      select(Number, Sample, Date, Days, Tank, Temperature, PAR) %>%
      as.data.frame() %>%
      `rownames<-`( rownames(sequence_table_clean) )
  ),
  tax_table(
    taxa_table_resolved
  )
) %T>%
  print()

# Simplify ASV names and save DNA in refseq
ASV %<>%
  taxa_names() %>%
  DNAStringSet() %>%
  `names<-`( ASV %>% taxa_names() ) %>%
  merge_phyloseq(ASV) %>%
  `taxa_names<-`( str_c("ASV", ASV %>% ntaxa() %>% seq()) ) %T>%
  print()

otu_table(ASV) # ASVs are now aptly named in both tables
tax_table(ASV)
refseq(ASV) # All sequences are preserved

# Simple test
ASV %>%
  plot_richness("Sample", measures = c("Shannon", "Simpson"))

# Save phyloseq object
ASV %>%
  write_rds(here("DNA", "RDS", "ASV.rds"))

# Save DADA2 pipeline metadata
sequences %>%
  select(Number, Sample, contains("Reads"), contains("Prop")) %>%
  write_rds(here("DNA", "RDS", "pipeline.rds"))

# Clean up
rm(list = ls())