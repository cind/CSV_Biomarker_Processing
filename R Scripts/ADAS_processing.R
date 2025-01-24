## 2 tables

## ADAS Scores, ADNI1: ADAS Sub-Scores and Total Scores [ADNI1]
## ADAS Scores, ADNIGO/2/3: Alzheimer's Disease Assessment Scale (ADAS) [ADNIGO,2,3]
library(tidyverse)

adas_1_score <- readr::read_delim("~/Projects/Amprion Project/Source Data/ADASSCORES_22May2024.csv") %>% # edit to file name and directory for ADAS Scores, ADNI 1
  dplyr::rename(VISCODE2 = VISCODE) 
adas_2_3_go <- readr::read_delim("~/Projects/Amprion Project/Source Data/ADAS_ADNIGO23_22May2024.csv") %>% # edit to file name and directory for ADAS Scores, ADNI GO/2/3
  dplyr::mutate(ADAS13=TOTAL13,EXAMDATE=VISDATE)

## recodes missing data for total ADAS scores in ADNI1
adas_1_score$ADAS13 <- ifelse(adas_1_score$TOTALMOD==-4,NA,adas_1_score$TOTALMOD)

## joins two ADAS tables
adas_scores <- dplyr::bind_rows(adas_2_3_go,adas_1_score) %>%
  dplyr::select(RID,VISCODE2,EXAMDATE,ADAS13) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)
