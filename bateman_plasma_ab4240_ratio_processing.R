## 3 tables

## Bateman 2022: 	Bateman Lab Plasma Abeta42/Abeta40 Ratio as a Predictor of Brain Amyloidosis [ADNI1,2,GO], Version 2022-11-18
## Bateman 2019: 	Bateman Lab Plasma Abeta42/Abeta40 Ratio as a Predictor of Brain Amyloidosis [ADNIGO,2], Version 2019-06-21
## FNIH: FNIH Biomarkers Consortium Plasma Abeta Project: Wash U[ADNI1,GO,2]


library(tidyverse)
library(lubridate)
library(mcr)

## select runs at primary injection site with a Lumos instrument that passed QC
bateman_plasma_nov_22 <-readr::read_delim("~/Projects/Amprion Project/Source Data/batemanlab_20221118.csv") # edit to your Bateman 2022 file location

colnames(bateman_plasma_nov_22) <- stringr::str_to_upper(colnames(bateman_plasma_nov_22))

bateman_plasma_nov_22 <- bateman_plasma_nov_22 %>%
  dplyr::rename(abeta_ratio=ABETA_4240_STANDARDIZED) %>%
  dplyr::filter(!is.na(RID),QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a",RID!=999999) 

bateman_plasma_nov_22$EXAMDATE<-lubridate::mdy(bateman_plasma_nov_22$EXAMDATE)
bateman_plasma_nov_22$MS_RUN_DATE<-lubridate::ymd(bateman_plasma_nov_22$MS_RUN_DATE)

bateman_plasma_jun_19 <-
  readr::read_delim("~/Projects/Amprion Project/Source Data/batemanlab_20190621.csv") %>% # edit to your Bateman 2019 file location
  dplyr::rename(abeta_ratio = RATIO_ABETA42_40_BY_ISTD_TOUSE) %>% 
  dplyr::mutate(origin="Bateman") %>% 
  dplyr::filter(BATCH_N!="Batch #00 Trial",QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a",RID!=999999) 

washu_plasma <-
  readr::read_delim("~/Projects/Amprion Project/Source Data/PLASMA_ABETA_PROJECT_WASH_U_11_05_21.csv") # edit to your FNIH file location
washu_plasma <-
  washu_plasma %>% dplyr::rename(abeta_ratio = STANDARDIZED_PLASMAAB4240,abeta_ratio_non_std = PLASMAAB4240) %>% dplyr::mutate(origin="WashU") %>% 
  dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a")

## calculate adjustment factors among AB42/40 datasets using Deming regression
## 2022 Bateman used as reference DF

deming_df_old_bateman<- dplyr::left_join(bateman_plasma_nov_22,
                                         bateman_plasma_jun_19 ,
                                         by=c("RID","VISCODE2"),suffix=c(".new",".old"))
deming_df_old_bateman <- deming_df_old_bateman %>% 
  tidyr::drop_na(abeta_ratio.new,abeta_ratio.old)

deming_reg_old_bateman <-mcr::mcreg(y=deming_df_old_bateman$abeta_ratio.new,x=deming_df_old_bateman$abeta_ratio.old,method.reg = "Deming")

deming_df_washu<- dplyr::left_join(bateman_plasma_nov_22,
                                   washu_plasma,
                                   by=c("RID","VISCODE2"),
                                   suffix=c(".new",".washu"))

deming_df_washu <- deming_df_washu %>% 
  dplyr::filter(RID!=999999) %>% 
  tidyr::drop_na(abeta_ratio.new,abeta_ratio.washu) %>% 
  dplyr::mutate(diff=abeta_ratio.new-abeta_ratio.washu)

deming_reg_washu <-mcr::mcreg(y=deming_df_washu$abeta_ratio.new,x=deming_df_washu$abeta_ratio.washu,method.reg = "Deming")

## assign adjusted AB42/40 values to FNIH and Bateman 2019 DFs
bateman_plasma_jun_19 <- bateman_plasma_jun_19 %>% 
  dplyr::rename(abeta_ratio_unadj=abeta_ratio) %>% 
  dplyr::mutate(abeta_ratio=deming_reg_old_bateman@para[1,1]+(deming_reg_old_bateman@para[2,1]*abeta_ratio_unadj))

washu_plasma <- washu_plasma %>% 
  dplyr::rename(abeta_ratio_unadj=abeta_ratio) %>% 
  dplyr::mutate(abeta_ratio=deming_reg_washu@para[1,1]+(deming_reg_washu@para[2,1]*abeta_ratio_unadj))

## merge datasets
bateman_merged <-dplyr::bind_rows(bateman_plasma_nov_22 %>% dplyr::mutate(origin="New Bateman"),
                                  bateman_plasma_jun_19 %>% dplyr::mutate(origin="Old Bateman"))

plasma_abeta_merged <- dplyr::bind_rows(
  washu_plasma %>% dplyr::mutate(MS_RUN_DATE=lubridate::ymd(MS_RUN_DATE)),
  bateman_merged %>% dplyr::mutate(VISCODE=VISCODE2))

## remove older replicates of QC records replicated in multiple runs
plasma_abeta_merged <- plasma_abeta_merged %>% dplyr::mutate(is_qc_record=((RID %in% deming_df_old_bateman$RID)|(RID %in% deming_df_washu$RID)))
plasma_abeta_merged <- plasma_abeta_merged %>% dplyr::filter(is_qc_record==FALSE|origin=="New Bateman")
plasma_abeta_merged <- plasma_abeta_merged %>% dplyr::select(RID,VISCODE,VISCODE2,EXAMDATE,origin,abeta_ratio,abeta_ratio_unadj,MS_RUN_DATE,is_qc_record)
plasma_abeta_merged <- plasma_abeta_merged %>% dplyr::mutate(sample_age = MS_RUN_DATE - EXAMDATE)
