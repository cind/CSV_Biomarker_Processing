## ADAS-Cognitive Behavior [ADNI1,GO,2,3,4]

library(tidyverse)

adas_scores <- readr::read_delim("~/Downloads/ADAS_26Jun2025.csv")

adas_scores <- adas_scores %>%
  dplyr::rename(EXAMDATE = VISDATE,
                ADAS10 = TOTSCORE,
                ADAS13 = TOTAL13) %>%
  dplyr::filter(!(is.na(ADAS10)&is.na(ADAS13))) %>%
  dplyr::select(RID,VISCODE2,EXAMDATE,ADAS10,ADAS13)

adas_scores <- adas_scores %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)