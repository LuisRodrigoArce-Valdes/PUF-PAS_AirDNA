# This is an auxiliatry script to read and tidy up the final results matrices
rm(list = ls())

# Calling up libraries
library(dplyr)
library(tidyr)
library(stringr)

insects <- list()
for(i in c("Deblur", "DADA2")){
  for(n in c("t0_l10","t20_l10","t20_l50","t26_l10")){
    name <- paste0(i,"_",n)
    meta <- read.delim(paste0("../results/03_Matrices/",name,"/metadata.tsv"), comment.char = "#")
    matrix <- read.delim(paste0("../results/02_feature-table_",name,".tsv"), skip = 1)
    
    # Merging
    meta %>% 
      left_join(matrix, join_by(Feature.ID == X.OTU.ID)) -> meta
    meta %>% 
      gather("Sample","Reads",c(4:ncol(meta))) %>% 
      select(Sample, Feature.ID, Confidence, Reads, Taxon) %>% 
      filter(Reads > 0) -> insects[[name]]
  }
}
rm(meta, matrix)

# Saving list
save(insects, file = "../../../02_Filtering/data/01_PUF_PAS_pt1.Rdata")
