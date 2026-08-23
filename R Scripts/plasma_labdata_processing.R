library(tidyverse)
library(transplantr)

source("~/Downloads/demographics_processing.R")
source("~/Downloads/adni_generate_VISCODE3_fn.R")

lab_data <- readr::read_delim("~/Downloads/LABDATA_10Aug2026.csv")

urmc_df <- readr::read_delim("~/Downloads/URMC_LABDATA_10Aug2026.csv")

urmc_df <- urmc_df %>%
  dplyr::filter(TestID %in% c("CRE2","GLU","TP","CHOL"),
                TestName %in% c("Creatinine","Glucose","Total Protein","Cholesterol")) %>%
  dplyr::select(RID,VISCODE2,SampleDate,TestID,TestName,ResultValue) %>%
  tidyr::pivot_wider(names_from = TestName,
                     values_from = ResultValue,
                     id_cols = c(RID,VISCODE2,SampleDate)) %>%
  dplyr::mutate(
    cholesterol = as.numeric(Cholesterol),
    creatinine = as.numeric(Creatinine),
    total_protein = as.numeric('Total Protein'),
    glucose = as.numeric(Glucose)
  ) %>%
  dplyr::rename(lab_date = SampleDate) %>%
  dplyr::select(RID,VISCODE2,lab_date,
                cholesterol,creatinine,total_protein,glucose)

lab_data <- lab_data %>%
  dplyr::mutate(glucose = as.numeric(RCT11),
                glucose = case_when(
                  glucose < 0 ~ NA,
                  TRUE ~ glucose
                ),
                total_protein = as.numeric(RCT12),
                total_protein = case_when(
                  total_protein < 0 ~ NA,
                  TRUE ~ total_protein
                ),
                creatinine = as.numeric(RCT392),
                creatinine = case_when(
                  creatinine < 0 ~ NA,
                  TRUE ~ creatinine
                ),
                cholesterol = as.numeric(RCT20),
                cholesterol = case_when(
                  cholesterol < 0 ~ NA,
                  TRUE ~ cholesterol
                )
  )

lab_data <- lab_data %>%
  dplyr::group_by(RID) %>%
  dplyr::summarize(glucose = mean(glucose,na.rm=TRUE),
                   total_protein = mean(total_protein,na.rm=TRUE),
                   creatinine = mean(creatinine,na.rm=TRUE),
                   cholesterol = mean(cholesterol,na.rm=TRUE),
                   lab_date = as.Date(mean(EXAMDATE,na.rm=TRUE)))

lab_data <- dplyr::left_join(lab_data,demographics,by="RID")

lab_data$age_at_lab <-
  round(as.numeric(lubridate::year(
    lab_data$lab_date
  )) -
    lab_data$birth_year +
    ((
      as.numeric(
        lubridate::month(lab_data$lab_date)
      ) - as.numeric(lab_data$birth_month)
    ) / 12),
  1)

lab_data$race_for_egfr <-
  ifelse(
    lab_data$race == "Black or African American",
    "black",
    "non-black"
  )

lab_data$gender_for_egfr <-
  ifelse(
    lab_data$sex == "Male",
    "M",
    "F"
  )

lab_data$eGFR <-
  transplantr::ckd_epi(
    lab_data$creatinine,
    lab_data$age_at_lab,
    lab_data$gender_for_egfr,
    lab_data$race_for_egfr,
    units = "US"
  )

lab_data <- lab_data %>%
  dplyr::select(RID,glucose,total_protein,creatinine,cholesterol,
                lab_date,eGFR
  )

## read in registry to calculate alternate VISCODE (VISCODE3)
registry <- readr::read_delim("~/Downloads/REGISTRY_10Aug2026.csv") # change path to your path

lab_data <- adni_generate_VISCODE3(df = lab_data,registry_df = registry,
                                 value_cols = c("glucose","total_protein",
                                                "creatinine","cholesterol"),
                                 date_col = lab_data$lab_date)

lab_data <- lab_data %>%
  dplyr::mutate(diabetes_bl_ind = as.factor(case_when(
    glucose > 125 ~ 1,
    glucose <= 125 ~ 0
  )),
  hyperlipidemia_bl_ind = as.factor(case_when(
    cholesterol > 200 ~ 1,
    cholesterol <= 200 ~ 0
  )),
  ckd_ind = as.factor(case_when(
    (eGFR <= 60) ~ 1,
    eGFR > 60 ~ 0)))