## table name: 	Mini-Mental State Examination (MMSE) [ADNI1,GO,2,3,4]

library(tidyverse)

mmse <- readr::read_delim("~/Projects/Amprion Project/Source Data/MMSE_22May2024.csv") %>% # edit to your source file name + location 
  dplyr::mutate(MMSE=MMSCORE) %>%
  dplyr::mutate(VISCODE2=case_when(
    VISCODE2 == "sc" ~ "bl",
    TRUE ~ VISCODE2
  )) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)
