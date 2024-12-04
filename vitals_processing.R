## table name: Vital Signs [ADNI1,GO,2,3,4]

library(tidyverse)

vitals <- readr::read_delim("~/Downloads/All_Subjects_VITALS_15Oct2024.csv") %>% ## change this file path to your own path
  dplyr::mutate(height_metric = case_when(
    VSHEIGHT < 0 ~ NA,
    VSHTUNIT == 2 ~ VSHEIGHT,
    VSHTUNIT == 1 ~ VSHEIGHT*2.54,
    TRUE ~ NA),
    weight_metric = case_when(
      VSWTUNIT == 2 ~ VSWEIGHT,
      VSWTUNIT == 1 ~ VSWEIGHT/2.205,
      TRUE ~ NA),
    height_imperial = case_when(
      VSHEIGHT < 0 ~ NA,
      VSHTUNIT == 1 ~ VSHEIGHT,
      VSHTUNIT == 2 ~ VSHEIGHT/2.54,
      TRUE ~ NA),
    weight_imperial = case_when(
      VSWTUNIT == 1 ~ VSWEIGHT,
      VSWTUNIT == 2 ~ VSWEIGHT*2.205,
      TRUE ~ NA)) %>%
  dplyr::mutate(height_metric = case_when(
    height_metric < 121 ~ NA,
    height_metric > 214 ~ NA,
    TRUE ~ height_metric 
  ),
  height_imperial = case_when(
    height_imperial > 84 ~ NA,
    height_imperial < 48 ~ NA,
    TRUE ~ height_imperial),
  diastolic_bp = as.numeric(VSBPDIA),
  systolic_bp = as.numeric(VSBPSYS),
  diastolic_bp = case_when(
    diastolic_bp <= 0 ~ NA,
    TRUE ~ diastolic_bp
  ),
  sysstolic_bp = case_when(
    systolic_bp <= 0 ~ NA,
    TRUE ~ systolic_bp
  ),
  hypertension_cat = case_when(
    diastolic_bp >= 90 | systolic_bp >= 140 ~ "Stage 2 Hypertension",
    diastolic_bp >= 80 | systolic_bp >= 130 ~ "Stage 1 Hypertension",
    diastolic_bp < 60 | systolic_bp < 90 ~ "Hypotension",
    diastolic_bp >= 60 & systolic_bp >= 90 ~ "Normal Blood Pressure"
  )
  )

vitals <- vitals %>%
  dplyr::group_by(RID) %>%
  dplyr::arrange(VISDATE) %>%
  tidyr::fill(c(height_metric,height_imperial),.direction="down") %>%
  dplyr::mutate(bmi = weight_metric/((height_metric/100)^2),
  abnormal_bmi_ind = (bmi>25|bmi<18.5),
  bmi_cat = case_when(
    bmi>30 ~ "Obese",
    bmi>25 ~ "Overweight",
    bmi>18.5 ~ "Normal Weight",
    bmi<=18.5 ~ "Underweight"
  ))

vitals <- vitals %>%
  dplyr::select(RID,VISCODE2,VISDATE,
                height_metric,weight_metric,height_imperial,weight_imperial,
                diastolic_bp,systolic_bp,hypertension_cat,
                bmi,abnormal_bmi_ind,bmi_cat)
