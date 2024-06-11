## table name: Blennow Lab ADNI1-2 Plasma neurofilament light (NFL) longitudinal [ADNI1,GO,2]
## data is longitudinal

library(tidyverse)

plasma_nfl <-
  readr::read_delim("~/Projects/Amprion Project/Source Data/ADNI_BLENNOWPLASMANFLLONG_10_03_18.csv") %>% 
  dplyr::select(RID,VISCODE2,EXAMDATE,PLASMA_NFL)

## Some subjects had two measurements at a given visit; averaged two measurements for those subjects
plasma_nfl <- plasma_nfl %>% 
  dplyr::group_by(RID,EXAMDATE) %>% 
  dplyr::mutate(PLASMA_NFL=mean(PLASMA_NFL)) %>% 
  dplyr::distinct_at(vars(RID,EXAMDATE),.keep_all = TRUE) 

plasma_nfl <- plasma_nfl %>% 
  dplyr::arrange(RID,EXAMDATE) %>% 
  dplyr::group_by(RID) %>% 
  dplyr::mutate(plasma_nfl_change_from_bl=PLASMA_NFL-PLASMA_NFL[1L],
                baseline_plasma_nfl=PLASMA_NFL[1L],
                plasma_nfl_time=EXAMDATE-EXAMDATE[1L]) %>% 
  dplyr::ungroup()

## recode reccord with incorrect VISCODE2
plasma_nfl$VISCODE2[which(plasma_nfl$RID=="1097"&plasma_nfl$EXAMDATE=="2013-01-31")]<-"m72"
