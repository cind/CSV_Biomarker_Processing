## 	Subject Demographics [ADNI1,GO,2,3,4]

library(tidyverse)

demographics <- readr::read_delim("~/Downloads/PTDEMOG_26Apr2026.csv") # change path to your path

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
    education_cat = case_when(
      PTEDUCAT > 16 ~ "Any postgraduate education",
      PTEDUCAT == 16 ~ "Bachelor's",
      PTEDUCAT == 12 ~ "High school or GED",
      PTEDUCAT > 12 ~ "Associate's/Some college",
      PTEDUCAT > 0 ~ "Less than high school"
    ),
    education_cat = factor(education_cat,levels = c("Less than high school",
                                                    "High school or GED",
                                                    "Associate's/Some college",
                                                    "Bachelor's",
                                                    "Any postgraduate education")),
    PTDOBMM = stringr::str_split_i(PTDOB,"/",1)
      )

## calculate age based on the difference in the month and year terms between the subject's birthdate and the visit date
demographics$age_bl <-round(as.numeric(lubridate::year(demographics$VISDATE)) -
                              as.numeric(lubridate::year(demographics$PTDOBYY)) +
          ((as.numeric(lubridate::month(demographics$VISDATE)) - 
              as.numeric(demographics$PTDOBMM)) / 12), digits=1)

demographics <- demographics %>%
  dplyr::rename(sex = PTGENDER,
               birth_year = PTDOBYY,
               birth_month = PTDOBMM,
               ethnicity = PTETHCAT,
               race = PTRACCAT,
               demog_date = VISDATE)

demographics <- demographics %>%
  dplyr::mutate(birth_year = lubridate::year(birth_year))

demographics <- demographics %>%
  dplyr::arrange(demog_date) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

demographics <- demographics %>%
  dplyr::select(RID,VISCODE2,sex,age_bl,education_cat,PTEDUCAT,ethnicity,race,birth_year,birth_month,demog_date)
