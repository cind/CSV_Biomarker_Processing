## Clinical Utility of Blood Biomarker Trajectories Study [ADNI1,2,GO,3]

library(tidyverse)

cubbt_biomarkers <- readr::read_delim("~/Downloads/FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES_16Jul2024.csv") # change path to your path

## split biomarker string into variables with analyte and manufacturer information
cubbt_biomarkers <- cubbt_biomarkers %>%
  dplyr::mutate(manufacturer = stringr::str_split_i(PLASMA_BIOMARKER,"_plasma_",1),
                analyte = stringr::str_split_i(PLASMA_BIOMARKER,"_plasma_",2),
                analyte_clean = case_when(
                  analyte %in% c("Abeta40","Ab40") ~ "AB40",
                  analyte %in% c("Abeta42","Ab42") ~ "AB42",
                  TRUE ~ analyte
                )
  )

## drop records with no corresponding ADNI RID
cubbt_biomarkers <- cubbt_biomarkers %>%
  dplyr::filter(!is.na(RID))

## drop records with unusually high coefficient of variation
cubbt_biomarkers <- cubbt_biomarkers %>%
  dplyr::filter(is.na(COMMENTS)|COMMENTS!="High CV (>25%)")

## preserve original dataframe
cubbt_biomarkers_long <- cubbt_biomarkers

## pivot data format to wide (one observation per subject visit)
cubbt_biomarkers <- cubbt_biomarkers_long %>%
  tidyr::pivot_wider(.,id_cols = c(RID,VISCODE2,EXAMDATE,manufacturer),
                     names_from = analyte_clean,
                     values_from = TESTVALUE)

## add "_plasma_cubbt" tag to variable names
cubbt_biomarkers <- cubbt_biomarkers %>%
  dplyr::rename_with(~paste0(.,"_plasma_cubbt"),
                     c(colnames(cubbt_biomarkers)[4:length(colnames(cubbt_biomarkers))]))

## generate AB42/40 ratio where not calculated but AB42 and AB40 are available
cubbt_biomarkers <- cubbt_biomarkers %>%
  dplyr::mutate(Abeta42_Abeta40_plasma_cubbt = case_when(
    !is.na(Abeta42_Abeta40_plasma_cubbt) ~ Abeta42_Abeta40_plasma_cubbt,
              !is.na(AB42_plasma_cubbt) & !is.na(AB40_plasma_cubbt) ~ AB42_plasma_cubbt/AB40_plasma_cubbt))

## merge duplicates within manufacturer
cubbt_biomarkers <- cubbt_biomarkers %>%
  dplyr::group_by(RID,VISCODE2,manufacturer_plasma_cubbt) %>%
  dplyr::mutate(across(ends_with("_cubbt") & where(is.numeric),.fns = ~mean(.,na.rm=TRUE))) %>%
  dplyr::ungroup() %>%
  dplyr::distinct_at(vars(RID,VISCODE2,manufacturer_plasma_cubbt),.keep_all = TRUE)