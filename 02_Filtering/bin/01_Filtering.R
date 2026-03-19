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
# Reading metabarcoding data
meta <- list()
for(i in c("pt1","pt2")){
  load(paste0("../data/01_PUF_PAS_",i,".Rdata"))
  bind_rows(insects, .id = "Filter") -> meta[[i]]
}
bind_rows(meta, .id = "Library") -> insects
rm(meta)

# Removing mont.tremblant and insectarium
insects %>% 
  filter(!grepl("Ins", Sample)) %>% 
  filter(Sample!="Mont.Tremblant") -> insects

# Reading metadata and removing extra samples
meta <- read.csv("../data/metadata.csv")
meta %>% 
  filter(!grepl("Ins", ID_R)) %>% 
  filter(ID_R!="Mont.Tremblant") -> meta

read.csv("../data/mock_PUFs.csv") %>% 
  select(Species) %>% 
  mutate(Species = gsub(" ","_",Species)) -> positive.species

positive.species$Species -> positive.species

# Estimating sample size (without controls)
meta %>% 
  filter(Type == "Sample") %>%
  select(ID_R) %>% 
  reframe(samples = unique(ID_R)) %>% 
  nrow() -> Os

# Checking string compatibility between the metadata and the metabarcoding data
identical(sort(unique(meta$ID_R)), sort(unique(insects$Sample)))

# Tidying and exploring
insects %>%
  left_join(unique(meta[,c("ID_R","Type")]), by = join_by(Sample == ID_R)) %>% 
  select(Library, Filter, Sample, Type, Feature.ID, Confidence, Reads, Taxon) -> insects

# 01. Raw data stats ####
# Estimating reads per denoising algorithm
insects %>%
  select(Filter, Feature.ID) %>% 
  unique() %>%  
  group_by(Filter) %>%
  summarise(n = n()) %>% 
  mutate(PostFilter = "0.Raw") %>% 
  select(Filter, PostFilter, n) -> counts

# Estimating number of samples per denoising algorithm
insects %>% 
  select(Filter, Sample, Type) %>%
  filter(Type == "Sample") %>%
  unique() %>%
  group_by(Filter) %>% 
  summarise(inds = n()) %>%
  mutate(inds = round(inds/Os*100,1)) %>% 
  right_join(counts, by = join_by(Filter)) %>% 
  select(PostFilter, Filter, n, inds) -> counts

# Viewing positive controls
insects %>% 
  filter(Type == "C+") %>% 
  select(Library, Filter, Sample, Taxon) %>% 
  unique() %>% 
  mutate(PS=str_detect(Taxon, paste(positive.species, collapse = "|"))) %>% 
  group_by(Filter, Library, Sample, PS) %>% 
  summarise(n = n()) %>% 
  mutate(PostFilter = "0.Raw") %>% 
  select(Filter, PostFilter, Library, Sample, PS, n) -> pc

# Checking presence of positive species in other samples
insects %>% 
  filter(Type != "C+") %>%
  filter(!str_detect(Sample,"Ins")) %>% 
  mutate(PS=str_detect(Taxon, paste(positive.species, collapse = "|"))) %>%
  filter(PS == T) %>% 
  mutate(PostFilter = "0.Raw") %>% 
  select(Filter, PostFilter, Sample, Feature.ID, Confidence, Reads, Taxon, PS) -> pc.cont

# 02. Negative controls filtering ####
insects %>% 
  filter(Type == "C-") -> nc

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
  filter(Type == "C-")

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
  filter(Type == "Sample") %>%
  unique() %>%
  group_by(Filter) %>% 
  summarise(inds = n()) %>%
  mutate(inds = round(inds/Os*100,1)) %>% 
  right_join(counts.f, by = join_by(Filter)) %>% 
  select(PostFilter, Filter, n, inds) %>%
  rbind(counts) %>% 
  arrange(PostFilter, Filter) -> counts

rm(counts.f, nc)

# Viewing positive controls
insects %>% 
  filter(Type == "C+") %>%
  select(Library, Filter, Sample, Taxon) %>% 
  unique() %>% 
  mutate(PS=str_detect(Taxon, paste(positive.species, collapse = "|"))) %>% 
  group_by(Filter, Library, Sample, PS) %>% 
  summarise(n = n()) %>% 
  mutate(PostFilter = "1.NC") %>% 
  select(Filter, PostFilter, Library, Sample, PS, n) %>% 
  rbind(pc) %>% 
  arrange(Filter, Library, Sample, PostFilter, desc(PS)) -> pc

# Checking presence of positive species in other samples
insects %>% 
  filter(Type != "C+") %>%
  filter(!str_detect(Sample,"Ins")) %>% 
  mutate(PS=str_detect(Taxon, paste(positive.species, collapse = "|"))) %>%
  filter(PS == T) %>% 
  mutate(PostFilter = "1.NC") %>% 
  select(Filter, PostFilter, Sample, Feature.ID, Confidence, Reads, Taxon, PS) %>% 
  rbind(pc.cont) %>% 
  arrange(Filter, Sample, PostFilter) -> pc.cont

# 03. Sample % filtering ####
insects %>%
  group_by(Library, Filter, Sample) %>% 
  mutate(sample.per=Reads/sum(Reads)*100) -> insects

thresholds <- c("0.000","0.025","0.050", "0.075", "0.100", "0.250","0.500","0.750","1.000")

for(i in thresholds) { # Sample percent filtering
  name <- paste0("2.sample.per.",i)
  
  i <- as.numeric(i)
  
  # Summarising filtering
  insects %>%
    filter(sample.per > i) %>%
    ungroup() %>% 
    select(Filter, Feature.ID) %>% 
    unique() %>%
    group_by(Filter) %>%
    summarise(n = n()) %>% 
    mutate(PostFilter = name) %>% 
    select(Filter, PostFilter, n) %>% 
    arrange(Filter, PostFilter) -> counts.f
  
  insects %>%
    filter(sample.per > i) %>%
    filter(Type == "Sample") %>%
    select(Library, Filter, Sample) %>%
    unique() %>%
    group_by(Filter) %>% 
    summarise(inds = n()) %>% 
    mutate(inds = round(inds/Os*100,1)) %>% 
    right_join(counts.f, by = join_by(Filter)) %>% 
    select(PostFilter, Filter, n, inds) %>%
    rbind(counts) %>% 
    arrange(PostFilter, Filter) -> counts
  
  rm(counts.f)
  
  # Viewing positive controls
  insects %>%
    filter(sample.per > i) %>%  
    filter(Type == "C+") %>% 
    select(Library, Filter, Sample, Taxon) %>% 
    unique() %>% 
    mutate(PS=str_detect(Taxon, paste(positive.species, collapse = "|"))) %>% 
    group_by(Filter, Library, Sample, PS) %>% 
    summarise(n = n()) %>% 
    mutate(PostFilter = name) %>% 
    select(Filter, PostFilter, Library, Sample, PS, n) %>% 
    rbind(pc) %>% 
    arrange(Filter, Library, Sample, PostFilter, desc(PS)) -> pc
  
  # Checking presence of positive species in other samples
  insects %>%
    ungroup() %>% 
    filter(sample.per > i) %>%  
    filter(Type != "C+") %>%
    filter(!str_detect(Sample,"Ins")) %>% 
    mutate(PS=str_detect(Taxon, paste(positive.species, collapse = "|"))) %>%
    filter(PS == T) %>% 
    mutate(PostFilter = name) %>% 
    select(Filter, PostFilter, Sample, Feature.ID, Confidence, Reads, Taxon, PS) %>% 
    rbind(pc.cont) %>% 
    arrange(Filter, Sample, PostFilter) -> pc.cont
}

# Plotting filtering schemes
for(i in thresholds) {
  name <- paste0("Deblur_t20_l10_sample.per.",i)
  
  i <- as.numeric(i)
  
  insects %>%
    filter(Filter == "Deblur_t20_l10") %>%
    filter(sample.per > i) %>%
    select(Library, Filter, Sample, Feature.ID, Confidence, Reads, Taxon, sample.per) -> insects.f
  
  # Finally, we will merge both dataets to explore the database in detail
  insects.f %>%
    separate(Taxon, into = c("k","p","c","o","f","g","s"), sep = "; ", extra = "merge") %>% 
    left_join(na.omit(meta[,c("ID_R","ID_order","site","from","to")]), by = join_by(Sample == ID_R)) -> insects.f # Joining with the fild database
  
  # Replacing empty filds by NA
  as.data.frame(lapply(insects.f, function(x) gsub(".?__","", x))) -> insects.f
  insects.f[insects.f == ""] <- NA
  
  # Tidying
  insects.f %>% 
    select(Library, Sample, ID_order, site, from, to, Feature.ID, Confidence, Reads, k, p, c, o, f, g, s) %>% 
    mutate(Reads = as.numeric(Reads))-> insects.f
  
  
  write.csv(insects.f, paste0("../results/A_Database_",name,".csv"), row.names = F)
  
  # Summary tables of fraction of reads per taxonomic level
  insects.f %>%
    filter(k == "Metazoa_33208") %>% 
    pivot_longer(cols = c("k", "p", "c", "o", "f", "g", "s"), names_to = "Level", values_to = "Taxon") %>% 
    mutate(Level = factor(Level, levels = c("k", "p", "c", "o", "f", "g", "s"))) %>% 
    select(Level, site, ID_order, from, to, Taxon, Reads) %>% 
    group_by(Level, site, ID_order, from, to, Taxon) %>% 
    summarise(Reads = sum(Reads)) %>% 
    mutate(Fx = Reads / sum(Reads)) %>%
    arrange(Level, site, ID_order, desc(Fx)) %>% 
    write.csv(paste0("../results/B_Summary_",name,".csv"), row.names = F)
}

# 04. Summarising filtering schemes ####
# Merging dataframes into a final summary one that can be added to the manuscript
counts %>% 
  mutate(FilFil = paste(PostFilter, Filter, sep = "-")) %>% 
  select(FilFil, PostFilter, Filter, n, inds) %>% 
  rename(ASVs = n) -> counts

pc %>%
  ungroup() %>% 
  mutate(Sample = paste(Sample,PS,sep = "_")) %>% 
  select(PostFilter, Filter, Sample, n) %>% 
  pivot_wider(names_from = Sample, values_from = n, values_fill = 0) %>% 
  mutate(FilFil = paste(PostFilter, Filter, sep = "-")) %>% 
  select(!c(PostFilter, Filter)) -> pc

# Merging
counts %>% 
  left_join(pc, by = join_by(FilFil)) -> counts
rm(pc)

# Now the number of positive control seqs in the samples
pc.cont %>% 
  select(PostFilter, Filter, Sample, Feature.ID) %>%
  group_by(PostFilter, Filter) %>% 
  summarise(n = n()) %>% 
  mutate(FilFil = paste(PostFilter, Filter, sep = "-")) %>%
  ungroup() %>% 
  select(FilFil, n) %>% 
  rename(MockASVs = n) -> pc.cont

#  Merging
counts %>% 
  left_join(pc.cont, by = join_by(FilFil)) %>% 
  mutate(MockASVs = replace_na(MockASVs, 0)) %>% 
  select(!FilFil) %>% 
  filter(!str_detect(PostFilter,"3.")) %>% 
  mutate(PostFilter = gsub("2\\.","",PostFilter)) %>% 
  mutate(PostFilter = gsub("4\\.","",PostFilter)) %>%
  relocate(Filter, PostFilter) %>% 
  arrange(Filter, PostFilter) -> counts

# Write csv
write.csv(counts, "../results/00_Filtering_summaries.csv", row.names = F)

# Creating a README
sink("../results/README.txt")
print("I suggest to use Deblur_t20_l10 after cleaning with negative controls (sample filtering 0.000), to recover the maximum number of taxa.")
print("And, 0.100 to clean the positive controls as much as possible without removing their known taxa")
sink()