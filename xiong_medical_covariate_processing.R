library(tidyverse)

hachinski <- readr::read_delim("~/Downloads/All_Subjects_MODHACH_15Oct2024.csv") %>%
  dplyr::mutate(hypertension_ind = case_when(
    HMHYPERT == 1 ~ 1,
    TRUE ~ 0
  ),
  stroke_ind = case_when(
    HMSTROKE == 2 ~ 1,
    TRUE ~ 0
  )
  )

vitals <- readr::read_delim("~/Downloads/All_Subjects_VITALS_15Oct2024.csv") %>% 
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
    bmi = weight_metric/((height_metric/100)^2),
    bmi_ind = (bmi>25|bmi<18.5))

lab_data <- readr::read_delim("~/Downloads/All_Subjects_LABDATA_10Oct2024.csv")
lab_data <- lab_data %>%
  dplyr::mutate(glucose = as.numeric(RCT11),
                glucose = case_when(
                  glucose < 0 ~ NA,
                  TRUE ~ glucose
                ))
                
lab_data <- lab_data %>%
                  dplyr::group_by(RID) %>%
                  dplyr::summarize(glucose = mean(glucose,na.rm=TRUE)) %>%
  dplyr::mutate(diabetes_ind = as.factor(case_when(
                  glucose > 125 ~ 1,
                  TRUE ~ 0
                )))

medical_covariates <- dplyr::full_join(hachinski %>%
                                         dplyr::select(RID,hypertension_ind,stroke_ind),
                                       vitals %>%
                                         dplyr::filter(VISCODE2 %in% c("sc","f")) %>%
                                         dplyr::select(RID,bmi,bmi_ind),
                                       by = "RID"
                                         )

medical_covariates <- dplyr::full_join(medical_covariates,
                                       lab_data,
                                       by = "RID") 

medical_covariates <- medical_covariates %>%
  dplyr::mutate(VISCODE2="bl")

