rm(list = ls())
# Script to load, merge, filter and tidy metabarcoding data and field metadata
# Filtering is based on:
# https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13780

# Loading libraries
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)


# 00. Reading inputs ####
# Original sample size (without controls)
Os <- 9

load("../data/01_PUF_PAS.Rdata")

# Tidying and exploring
insects %>% 
  bind_rows(.id = "Filter") %>% 
  mutate(Type = ifelse(Sample=="CTRL.neg.extraction.PUF" | Sample=="ctrl.PCR.neg.PUF","Control-","Sample")) %>% 
  select(Filter, Sample, Type, Feature.ID, Confidence, Reads, Taxon) -> insects

# 01. Raw data stats ####
# Estimating reads per denoising algorithm
insects %>%
  select(Filter, Feature.ID) %>% 
  unique() %>%  
  group_by(Filter) %>%
  summarise(n = n()) %>% 
  mutate(PostFilter = "0.Raw") %>% 
  select(Filter, PostFilter, n) -> counts

insects %>% 
  select(Filter, Sample, Type) %>%
  filter(Type != "Control-") %>% 
  unique() %>%
  group_by(Filter) %>% 
  summarise(inds = n()) %>%
  mutate(inds = round(inds/Os*100,1)) %>% 
  right_join(counts, by = join_by(Filter)) %>% 
  select(PostFilter, Filter, n, inds) -> counts

# 02. Negative controls filtering ####
insects %>% 
  filter(Type == "Control-") -> nc

for(i in unique(nc$Filter)) {
  tmp <- nc[nc$Filter == i,]
  for(r in 1:nrow(tmp)){
    Tax = tmp[r,"Taxon"]
    reads = tmp[r,"Reads"]
    
    insects %>% 
      filter(!(Filter==i & Taxon == Tax & Reads <= reads)) -> insects
  }
}
rm(tmp, i, r, reads, Tax)

# Negative controls shoould be clean
insects %>% 
  filter(Type == "Control-")

# Summarising negative controls filtering
insects %>%
  select(Filter, Feature.ID) %>% 
  unique() %>%
  group_by(Filter) %>%
  summarise(n = n()) %>% 
  mutate(PostFilter = "1.NC") %>% 
  select(Filter, PostFilter, n) %>% 
  arrange(Filter, PostFilter) -> counts.f

insects %>% 
  select(Filter, Sample, Type) %>%
  filter(Type != "Control-") %>% 
  unique() %>%
  group_by(Filter) %>% 
  summarise(inds = n()) %>%
  mutate(inds = round(inds/Os*100,1)) %>% 
  right_join(counts.f, by = join_by(Filter)) %>% 
  select(PostFilter, Filter, n, inds) %>%
  rbind(counts) %>% 
  arrange(PostFilter, Filter) -> counts

rm(counts.f, nc)

# 03. Exploring and plotting ####
filt <- "Deblur_t26_l10"
insects %>% 
  filter(Filter == filt) %>% 
  separate(Taxon, into = c("k","p","c","o","f","g","s"), sep = "; ", extra = "merge") -> insects.f

# Replacing empty filds by NA
as.data.frame(lapply(insects.f, function(x) gsub(".?__","", x))) -> insects.f
insects.f[insects.f == ""] <- NA

# Exploring
insects.f %>% 
  select(Sample, Feature.ID, Confidence, Reads, k, p, c, o, f, g, s) -> insects.f

# Plotting
insects.f %>%
  pivot_longer(cols = c("k", "p", "c", "o", "f", "g", "s"), names_to = "Level", values_to = "Taxon") %>% 
  mutate(Level = factor(Level, levels = c("k", "p", "c", "o", "f", "g", "s"))) %>% 
  select(Sample, Level, Taxon) %>% 
  unique() %>% 
  filter(Level == "k") %>%
  group_by(Sample, Taxon) %>% 
  summarise(count = n(), .groups = "drop_last") %>% # Count occurrences, then drop the last grouping level
  mutate(relative_frequency = count / sum(count)) %>% # Calculate relative frequency within each group_var
  ungroup() %>% 
    ggplot() +
    geom_col(aes(x=Sample, fill = Taxon, y=relative_frequency))

insects.f %>%
  filter(k == "Metazoa_33208") %>% 
  pivot_longer(cols = c("k", "p", "c", "o", "f", "g", "s"), names_to = "Level", values_to = "Taxon") %>% 
  mutate(Level = factor(Level, levels = c("k", "p", "c", "o", "f", "g", "s"))) %>% 
  select(Sample, Level, Taxon) %>% 
  unique() %>% 
  filter(Level == "p") %>%
  group_by(Sample, Taxon) %>% 
  summarise(count = n(), .groups = "drop_last") %>% # Count occurrences, then drop the last grouping level
  mutate(relative_frequency = count / sum(count)) %>% # Calculate relative frequency within each group_var
  ungroup() %>% 
  ggplot() +
  geom_col(aes(x=Sample, fill = Taxon, y=relative_frequency))

insects.f %>%
  filter(k == "Metazoa_33208") %>% 
  pivot_longer(cols = c("k", "p", "c", "o", "f", "g", "s"), names_to = "Level", values_to = "Taxon") %>% 
  mutate(Level = factor(Level, levels = c("k", "p", "c", "o", "f", "g", "s"))) %>% 
  select(Sample, Level, Taxon) %>% 
  unique() %>% 
  filter(Level == "c") %>%
  group_by(Sample, Taxon) %>% 
  summarise(count = n(), .groups = "drop_last") %>% # Count occurrences, then drop the last grouping level
  mutate(relative_frequency = count / sum(count)) %>% # Calculate relative frequency within each group_var
  ungroup() %>% 
  ggplot() +
  geom_col(aes(x=Sample, fill = Taxon, y=relative_frequency))

insects.f %>%
  filter(k == "Metazoa_33208") %>% 
  pivot_longer(cols = c("k", "p", "c", "o", "f", "g", "s"), names_to = "Level", values_to = "Taxon") %>% 
  mutate(Level = factor(Level, levels = c("k", "p", "c", "o", "f", "g", "s"))) %>% 
  select(Sample, Level, Taxon) %>% 
  unique() %>% 
  filter(Level == "o") %>%
  group_by(Sample, Taxon) %>% 
  summarise(count = n(), .groups = "drop_last") %>% # Count occurrences, then drop the last grouping level
  mutate(relative_frequency = count / sum(count)) %>% # Calculate relative frequency within each group_var
  ungroup() %>% 
  ggplot() +
  geom_col(aes(x=Sample, fill = Taxon, y=relative_frequency))
