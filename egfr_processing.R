library(tidyverse)
library(transplantr)

source("~/Projects/demographics_processing.R")
source("~/Downloads/plasma_labdata_processing.R")

egfr_df <- dplyr::left_join(demographics,lab_data,by="RID")

egfr_df$age_at_lab <-
  round(as.numeric(lubridate::year(
    egfr_df$lab_date
  )) -
    egfr_df$birth_year +
    ((
      as.numeric(
        lubridate::month(egfr_df$lab_date)
      ) - as.numeric(egfr_df$birth_month)
    ) / 12),
  1)

egfr_df$race_for_egfr <-
  ifelse(
    egfr_df$race == "Black or African American",
    "black",
    "non-black"
  )

egfr_df$gender_for_egfr <-
  ifelse(
    egfr_df$sex == "Male",
    "M",
    "F"
  )

egfr_df$eGFR <-
  transplantr::ckd_epi(
    egfr_df$creatinine,
    egfr_df$age_at_lab,
    egfr_df$gender_for_egfr,
    egfr_df$race_for_egfr,
    units = "US"
  )

egfr_df <- egfr_df %>%
  mutate(ckd_ind = as.factor(case_when(
                           (eGFR <= 60) ~ 1,
                           eGFR > 60 ~ 0)))

egfr_df <- egfr_df %>%
  dplyr::select(RID,eGFR,ckd_ind)