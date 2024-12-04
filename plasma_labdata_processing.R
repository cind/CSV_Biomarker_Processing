library(tidyverse)

lab_data <- readr::read_delim("~/Downloads/All_Subjects_LABDATA_10Oct2024.csv")

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
                   lab_date = as.Date(mean(EXAMDATE,na.rm=TRUE))) %>%
  dplyr::mutate(diabetes_bl_ind = as.factor(case_when(
    glucose > 125 ~ 1,
    glucose <= 125 ~ 0
  )),
  hyperlipidemia_bl_ind = as.factor(case_when(
    cholesterol > 200 ~ 1,
    cholesterol <= 200 ~ 0
  )))