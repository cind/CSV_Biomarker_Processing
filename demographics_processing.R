## 	Subject Demographics [ADNI1,GO,2,3,4]

library(tidyverse)

demographics <- readr::read_delim("~/Downloads/All_Subjects_PTDEMOG_10Oct2024.csv") # change path to your path

demographics <- demographics %>% dplyr::select(RID,VISCODE,VISCODE2,VISDATE,
                                               PTGENDER,PTDOB,PTDOBYY,PTEDUCAT,
                                               PTETHCAT,PTRACCAT)

## replace numeric values with factor levels from data dictionary
demographics <- demographics %>%
  dplyr::mutate(
    PTGENDER = case_when(
      PTGENDER==1 ~ "Male",
      PTGENDER==2 ~ "Female"),
    PTETHCAT = case_when(
      PTETHCAT==1 ~ "Hispanic or Latino",
      PTETHCAT==2 ~ "Not Hispanic or Latino",
      PTETHCAT==3 ~ "Unknown"),
    PTRACCAT = case_when(
      PTRACCAT==1 ~ "American Indian or Alaskan Native",
      PTRACCAT==2 ~ "Asian",
      PTRACCAT==3 ~ "Native Hawaiian or Other Pacific Islander",
      PTRACCAT==4 ~ "Black or African American",
      PTRACCAT==5 ~ "White",
      PTRACCAT==6 ~ "More than one race",
      PTRACCAT==7 ~ "Unknown"),
    PTDOBMM = stringr::str_split_i(PTDOB,"/",1)
      )

## remove years of education below 0
demographics$PTEDUCAT <- ifelse(demographics$PTEDUCAT>0,demographics$PTEDUCAT,NA)

## calculate age based on the difference in the month and year terms between the subject's birthdate and the visit date
demographics$age_bl <-round(as.numeric(lubridate::year(demographics$VISDATE)) -
                              as.numeric(demographics$PTDOBYY) +
          ((as.numeric(lubridate::month(demographics$VISDATE)) - 
              as.numeric(demographics$PTDOBMM)) / 12), digits=1)

demographics <- demographics %>%
  dplyr::rename(sex = PTGENDER,
               birth_year = PTDOBYY,
               birth_month = PTDOBMM,
               yrs_education = PTEDUCAT,
               ethnicity = PTETHCAT,
               race = PTRACCAT,
               demog_date = VISDATE)

demographics <- demographics %>%
  dplyr::arrange(demog_date) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

demographics <- demographics %>%
  dplyr::select(RID,VISCODE2,sex,age_bl,yrs_education,ethnicity,race,birth_year,birth_month)