## table name on LONI: ApoE Genotyping - Results [ADNI1,GO,2,3]

library(tidyverse)

apoeres <- read.csv("~/Downloads/APOERES_31Jul2024.csv") 

## creates two additional utility variables for APOE information
## apoe4_status: whether a subject has at least one APOE-e4 allele
## apoe4_status_text: a leveled factor for each subject's APOE-e4 genotype
apoeres <- apoeres %>% 
  dplyr::mutate(apoe4_status=GENOTYPE %in% c("2/4","3/4","4/4"),
                                             apoe4_status_text = factor(case_when(
                                               GENOTYPE=="4/4" ~ "Homozygotes",
                                               GENOTYPE %in% c("2/4","3/4") ~ "Heterozygotes",
                                               TRUE ~ "Non-Carrier"),
                                               levels = c("Non-Carrier","Heterozygotes","Homozygotes"))) %>%
  dplyr::distinct_at(vars(RID),.keep_all = TRUE)

apoeres <- apoeres %>%
  dplyr::select(RID,GENOTYPE,apoe4_status_text,apoe4_status)
