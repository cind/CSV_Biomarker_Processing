# table name: Diagnostic Summary [ADNI1,GO,2,3,4]

library(tidyverse)

adni_diagnoses <-
  readr::read_delim("~/Projects/Amprion Project/Source Data/DXSUM_PDXCONV_21May2024.csv") %>% # change to your file name and location here 
  dplyr::select(RID,DIAGNOSIS,EXAMDATE,VISCODE2) %>%
  dplyr::mutate(DX = case_when(
    (DIAGNOSIS == 1) ~ "CU",
    (DIAGNOSIS == 2) ~ "MCI",
    (DIAGNOSIS == 3) ~ "Dementia")
  ) %>%
  dplyr::select(RID,DX,VISCODE2) %>%
  dplyr::filter(!is.na(DX))
