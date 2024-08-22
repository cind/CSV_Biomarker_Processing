## table name on LONI - ADNI1 Protocol 1.5T: 	UCSF - Cross-Sectional FreeSurfer (FreeSurfer Version 4.3) [ADNI1,GO,2]
## table name on LONI - ADNIGO/2 Protocol 3T: 	UCSF - Cross-Sectional FreeSurfer (5.1) [ADNI1,GO,2]
## table name on LONI - ADNI3 Protocol 3T: 	UCSF - Cross-Sectional FreeSurfer (6.0) [ADNI3]
## table name on LONI - Demographics:  	Subject Demographics [ADNI1,GO,2,3,4]

library(dplyr)
library(reshape)

# ComBat harmonization
source('~/projects/ComGamPackage-master/ComGamFunctionHelpers.R')
source('~/projects/ComGamPackage-master/ComGamHarmFunction.R')

#####################################################################################################################
# processing steps for ICV-adjusted, age-adjusted, harmonized ROI values (also creating cortical thickness meta-ROI)
#####################################################################################################################

# Diagnoses as CN, MCI, or AD
adni_diagnoses <- readr::read_delim("~/data/DXSUM_PDXCONV_23Jul2024.csv") %>% # change to your file name and location here 
  dplyr::select(RID,DIAGNOSIS,EXAMDATE, VISCODE2) %>%
  dplyr::rename(DX.DATE = EXAMDATE) %>%
  dplyr::mutate(DX = case_when(
    (DIAGNOSIS == 1) ~ "CU",
    (DIAGNOSIS == 2) ~ "MCI",
    (DIAGNOSIS == 3) ~ "Dementia")
  ) %>%
  dplyr::filter(!is.na(DX))

# Demographics for Gender, Education, Age
dem <- read.csv("~/data/demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))
dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))

# Amyloid PET data in order to do ICV-adjustment
amyloid_pet <- readr::read_delim("~/data/UCBERKELEY_AMY_6MM_31Jul2024.csv") [,1:16]
amyloid_pet <- amyloid_pet %>% dplyr::rename(suvr_summary=SUMMARY_SUVR,
                                             Centiloid=CENTILOIDS,
                                             AmyloidPosPET=AMYLOID_STATUS,
                                             EXAMDATE_pet = SCANDATE) %>%
  dplyr::mutate(EXAMDATE_pet = as.Date(EXAMDATE_pet),
                RID = as.character(RID),
                AmyloidPosPET_Centiloid = case_when(Centiloid >= 20 ~ 1,
                                                    TRUE ~ 0)) %>%
  dplyr::select(RID, EXAMDATE_pet, AmyloidPosPET, AmyloidPosPET_Centiloid, Centiloid, suvr_summary)

# getting amyloid negative RID's that were consistently negative for every scan they had
amyloid_negs <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(first_a_neg_date_pet = min(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, first_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

# ADNI 1 data
adni1_1.5 <- read.csv("~/data/freesurfer_4.3_ADNI1_protocol_1.5T.csv") %>%
  dplyr::mutate(STUDY = "ADNI1",
                Field_Strength = "1.5T")
adni1_1.5 <- adni1_1.5 %>%
  dplyr::filter(OVERALLQC == "Pass")

# getting demographics and closest diagnosis into the data
adni1_1.5 <- merge(adni1_1.5, dem, all.x = TRUE) %>% 
  dplyr::mutate(age = round(lubridate::time_length(difftime(EXAMDATE, birthdate), "years"), digits = 1))
adni1_1.5$RID <- as.factor(adni1_1.5$RID)
adni_diagnoses$RID <- as.factor(adni_diagnoses$RID)

adni1_1.5 <- merge(adni1_1.5, adni_diagnoses %>% 
                  dplyr::select(-VISCODE2) %>%
                  dplyr::filter(!is.na(DX.DATE)), 
                by = "RID") %>%
  dplyr::mutate(date_diff = abs(lubridate::time_length(difftime(EXAMDATE, DX.DATE), "years"))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::mutate(min_match_date = min(date_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(min_match_date == date_diff) %>%
  dplyr::distinct() %>%
  dplyr::rename(diags = DX)

adni1_1.5 <- adni1_1.5 %>%
  dplyr::rename(ICV = ST10CV)

# getting rid of variables with all NA
adni1_1.5 <- Filter(function(x)!all(is.na(x)), adni1_1.5)

# ADNI 2 data
adni2_3 <- read.csv("~/data/freesurfer_5.1_ADNI2GO_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI2",
                Field_Strength = "3T")
adni2_3 <- adni2_3 %>%
  dplyr::filter(OVERALLQC == "Pass")

# getting demographics and closest diagnosis into the data
adni2_3 <- merge(adni2_3, dem, all.x = TRUE) %>% 
  dplyr::mutate(age = round(lubridate::time_length(difftime(EXAMDATE, birthdate), "years"), digits = 1))
adni2_3$RID <- as.factor(adni2_3$RID)
adni_diagnoses$RID <- as.factor(adni_diagnoses$RID)

adni2_3 <- merge(adni2_3, adni_diagnoses %>% 
                        dplyr::select(-VISCODE2) %>%
                        dplyr::filter(!is.na(DX.DATE)), 
                      by = "RID") %>%
  dplyr::mutate(date_diff = abs(lubridate::time_length(difftime(EXAMDATE, DX.DATE), "years"))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::mutate(min_match_date = min(date_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(min_match_date == date_diff) %>%
  dplyr::distinct() %>%
  dplyr::rename(diags = DX)

adni2_3 <- adni2_3 %>%
  dplyr::rename(ICV = ST10CV) 

# getting rid of variables with all NA
adni2_3 <- Filter(function(x)!all(is.na(x)), adni2_3)
  
# ADNI 3 data
adni3_3 <- read.csv("~/data/freesurfer_6.0_ADNI3_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI3",
                Field_Strength = "3T")
adni3_3 <- adni3_3 %>%
  dplyr::filter(OVERALLQC == "Pass")

# getting demographics and closest diagnosis into the data
adni3_3 <- merge(adni3_3, dem, all.x = TRUE) %>% 
  dplyr::mutate(age = round(lubridate::time_length(difftime(EXAMDATE, birthdate), "years"), digits = 1))
adni3_3$RID <- as.factor(adni3_3$RID)
adni_diagnoses$RID <- as.factor(adni_diagnoses$RID)

adni3_3 <- merge(adni3_3, adni_diagnoses %>% 
                  dplyr::select(-VISCODE2) %>%
                  dplyr::filter(!is.na(DX.DATE)), 
                by = "RID") %>%
  dplyr::mutate(date_diff = abs(lubridate::time_length(difftime(EXAMDATE, DX.DATE), "years"))) %>%
  dplyr::group_by(RID, EXAMDATE) %>%
  dplyr::mutate(min_match_date = min(date_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(min_match_date == date_diff) %>%
  dplyr::distinct() %>%
  dplyr::rename(diags = DX)

adni3_3 <- adni3_3 %>%
  dplyr::rename(ICV = ST10CV)

# getting rid of variables with all NA
adni3_3 <- Filter(function(x)!all(is.na(x)), adni3_3)

#############################
# ICV-adjustments: getting CU, amyloid-negative scans
#############################

# getting earliest, amyloid-negative scan per RID
# ADNI 1
adni1_1.5_temp <- adni1_1.5 %>%
  dplyr::filter(diags == "CU")
adni1_1.5_temp <- adni1_1.5_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
adni1_1.5_temp <- adni1_1.5_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

# ADNI GO/2
adni2_3_temp <- adni2_3 %>%
  dplyr::filter(diags == "CU")
adni2_3_temp <- adni2_3_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
adni2_3_temp <- adni2_3_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

# ADNI 3
adni3_3_temp <- adni3_3 %>%
  dplyr::filter(diags == "CU")
adni3_3_temp <- adni3_3_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
adni3_3_temp <- adni3_3_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

#############################
# ICV-adjustments: setting up regression
#############################
# ADNI 1
# creating list of variables that need to be ICV-adjusted
volumes_adni1_1.5 <- adni1_1.5_temp %>%
  dplyr::select(ICV | contains("CV") | contains("SV"))

variables_adni1_1.5 <- c()

for (roi_num in 1:length(names(volumes_adni1_1.5))){
  current_roi <- names(volumes_adni1_1.5)[roi_num]
  variables_adni1_1.5 <- append(variables_adni1_1.5, current_roi)
}

for (i in variables_adni1_1.5[-1]){
  fit <- lm(adni1_1.5_temp[[i]] ~ ICV, data = adni1_1.5_temp, na.action=na.exclude)
  adni1_1.5[paste0(i,'_adjusted')] <- adni1_1.5[paste0(i)] - (adni1_1.5$ICV - mean(adni1_1.5_temp$ICV, na.rm = T)) * fit$coefficients[2]
}

adni1_1.5 <- adni1_1.5 %>%
  dplyr::select(-names(volumes_adni1_1.5))

# ADNI GO/2
# creating list of variables that need to be ICV-adjusted
volumes_adni2_3 <- adni2_3_temp %>%
  dplyr::select(ICV | contains("CV") | contains("SV"))

variables_adni2_3 <- c()

for (roi_num in 1:length(names(volumes_adni2_3))){
  current_roi <- names(volumes_adni2_3)[roi_num]
  variables_adni2_3 <- append(variables_adni2_3, current_roi)
}

for (i in variables_adni2_3[-1]){
  fit <- lm(adni2_3_temp[[i]] ~ ICV, data = adni2_3_temp, na.action=na.exclude)
  adni2_3[paste0(i,'_adjusted')] <- adni2_3[paste0(i)] - (adni2_3$ICV - mean(adni2_3_temp$ICV, na.rm = T)) * fit$coefficients[2]
}

adni2_3 <- adni2_3 %>%
  dplyr::select(-names(volumes_adni2_3))

# ADNI 3
# creating list of variables that need to be ICV-adjusted
volumes_adni3_3 <- adni3_3_temp %>%
  dplyr::select(ICV | contains("CV") | contains("SV"))

variables_adni3_3 <- c()

for (roi_num in 1:length(names(volumes_adni3_3))){
  current_roi <- names(volumes_adni3_3)[roi_num]
  variables_adni3_3 <- append(variables_adni3_3, current_roi)
}

for (i in variables_adni3_3[-1]){
  fit <- lm(adni3_3_temp[[i]] ~ ICV, data = adni3_3_temp, na.action=na.exclude)
  adni3_3[paste0(i,'_adjusted')] <- adni3_3[paste0(i)] - (adni3_3$ICV - mean(adni3_3_temp$ICV, na.rm = T)) * fit$coefficients[2]
}

adni3_3 <- adni3_3 %>%
  dplyr::select(-names(volumes_adni3_3))

#############################
# merging variables together
#############################

mri_list <- list(adni1_1.5, adni2_3, adni3_3)

mri_data <- merge_recurse(mri_list)

# ridding mri_data of columns with NA 
mri_data <- mri_data[ , apply(mri_data, 2, function(x) !any(is.na(x)))]

#############################
# ComBat Harmonization
#############################

all_features <- mri_data[, c(17:220, 232:336)] 

all_covariates <- mri_data[, c(1:5, 222:226, 229)] 

# now creating dataset for harmonization
all_data_combined_CN <- cbind(all_features, all_covariates) %>%
  dplyr::filter(diags == "CU")

all_features_CN <- all_data_combined_CN[, 1:309] 
all_covariates_CN <- all_data_combined_CN[, c(315, 317, 319)] 
extra_covariates_CN <- all_data_combined_CN[, c(310:314, 316, 318, 320)]
all_covariates_CN$STUDY <- as.factor(all_covariates_CN$STUDY)
all_covariates_CN$age <- as.numeric(all_covariates_CN$age)
all_covariates_CN$PTGENDER <- as.factor(all_covariates_CN$PTGENDER)

CN_data_harmonized <- ComGamHarm(feature.data = all_features_CN,
                                 covar.data = all_covariates_CN, 
                                 eb                = TRUE,
                                 parametric        = TRUE,
                                 smooth.terms      = c("age"),  # c("Examdate_Age"),
                                 k.val             = 5,
                                 verbose           = TRUE,
                                 model.diagnostics = FALSE)

CN_harmonized_data <- as.data.frame(t(CN_data_harmonized$harm.results))

CN_harmonized_data <- cbind(extra_covariates_CN, all_covariates_CN, CN_harmonized_data)

all_data_combined_MCI_AD <- cbind(all_features, all_covariates) %>%
  dplyr::filter(!diags == "CN")

all_features_MCI_AD <- all_data_combined_MCI_AD[, 1:309] 
all_covariates_MCI_AD <- all_data_combined_MCI_AD[, c(315, 317, 319)]
extra_covariates_MCI_AD <- all_data_combined_MCI_AD[, c(310:314, 316, 318, 320)]
all_covariates_MCI_AD$STUDY <- as.factor(all_covariates_MCI_AD$STUDY)
all_covariates_MCI_AD$age <- as.numeric(all_covariates_MCI_AD$age)
all_covariates_MCI_AD$PTGENDER <- as.factor(all_covariates_MCI_AD$PTGENDER)

MCI_AD_harmonized <- t(ApplyHarm(feature.data   = all_features_MCI_AD,
                                 covariate.data = all_covariates_MCI_AD,
                                 comgam.out     = CN_data_harmonized))

MCI_AD_harmonized <- as.data.frame(MCI_AD_harmonized)

MCI_AD_harmonized <- cbind(extra_covariates_MCI_AD, all_covariates_MCI_AD, MCI_AD_harmonized)

harmonized_data_from_CN <- rbind(CN_harmonized_data, MCI_AD_harmonized)

#############################
# Age Adjustment
#############################
mri_plot_data <- harmonized_data_from_CN %>%
  dplyr::distinct()

#getting amyloid-negative CN cases
mri_plot_data_temp <- mri_plot_data %>%
  dplyr::filter(diags == "CU")
mri_plot_data_temp <- mri_plot_data_temp %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
mri_plot_data_temp <- mri_plot_data_temp %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

# creating list of variables that need to be age-adjusted - for now just volumes and cortical thickness
features_mri_plot_data <- mri_plot_data %>%
  dplyr::select(contains("CV") | contains("TS") | contains("SV"))

variables_mri_plot_data <- c()

for (roi_num in 1:length(names(features_mri_plot_data))){
  current_roi <- names(features_mri_plot_data)[roi_num]
  variables_mri_plot_data <- append(variables_mri_plot_data, current_roi)
}

# regression for age-adjustment
for (i in variables_mri_plot_data){
  fit <- lm(mri_plot_data_temp[[i]] ~ age + poly(age, 2, raw = TRUE)[,"2"], data = mri_plot_data_temp, na.action=na.exclude)
  mri_plot_data[i] <- mri_plot_data[paste0(i)] - (mri_plot_data$age - mean(mri_plot_data_temp$age, na.rm = T)) * fit$coefficients[2]
}
