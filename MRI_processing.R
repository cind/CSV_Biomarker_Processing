## table name on LONI - ADNI1 Protocol 1.5T: 	UCSF - Cross-Sectional FreeSurfer (FreeSurfer Version 4.3) [ADNI1,GO,2]
## table name on LONI - ADNIGO/2 Protocol 3T: 	UCSF - Cross-Sectional FreeSurfer (5.1) [ADNI1,GO,2]
## table name on LONI - ADNI3 Protocol 3T: 	UCSF - Cross-Sectional FreeSurfer (6.0) [ADNI3]
## table name on LONI - Demographics:  	Subject Demographics [ADNI1,GO,2,3,4]

library(dplyr)

# ComBat harmonization
source('~/ComGamPackage-master/ComGamFunctionHelpers.R')
source('~/ComGamPackage-master/ComGamHarmFunction.R')

#####################################################################################################################
# processing steps for ICV-adjusted, age-adjusted, harmonized ROI values (also creating cortical thickness meta-ROI)
#####################################################################################################################

#first getting demographics (getting information to calculate age)
dem <- read.csv("~/demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))
dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))
dem <- merge(dem, apoeres, all = TRUE) %>% # apoe processing- steps are in other file
  dplyr::mutate(apoe = case_when(apoe == "E2/E3" | apoe == "E2/E2" ~ "E2",
                                 apoe == "E3/E3" ~ "E3",
                                 apoe == "E3/E4" | apoe == "E4/E4" ~ "E4"),
                RID = as.character(RID))

# calling in all freesurfer
adni1_1.5 <- read.csv("~/freesurfer_4.3_ADNI1_protocol_1.5T.csv") %>%
  dplyr::mutate(STUDY = "ADNI1",
                Field_Strength = "1.5T")
adni2_3 <- read.csv("~/freesurfer_5.1_ADNI2GO_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI2",
                Field_Strength = "3T")
adni3_3 <- read.csv("~/freesurfer_6.0_ADNI3_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI3",
                Field_Strength = "3T")

# creating one dataset with freesurfer data
mri_list <- list(adni1_1.5, adni2_3, adni3_3)
mri_data <- adni1_1.5
for ( df in mri_list ) {
  mri_data <-merge(mri_data, df, all=T)
}
mri_data <- mri_data %>%
  dplyr::filter(!(OVERALLQC == "Fail" | OVERALLQC == "" | OVERALLQC == " ")) %>%
  dplyr::distinct(across(c(-TEMPQC, -FRONTQC, -VENTQC, -PARQC, -INSULAQC, -OCCQC, -BGQC, -CWMQC, -Field_Strength, -HIPPOQC, -VISCODE, -IMAGETYPE, -update_stamp, -COLPROT, -VISCODE2, -LONISID, -FLDSTRENG, -RUNDATE, -STATUS, -VERSION, -LHIPQC, -RHIPQC)))

#merging demographics and freesurfer data
mri_data <- merge(mri_data, dem, all.x = TRUE) %>%
  dplyr::mutate(age = round(lubridate::time_length(difftime(EXAMDATE, birthdate), "years"), digits = 1)) %>%
  dplyr::rename(ST88CV = ST88SV, #these are hippocampal volumes - relabeling as CV regions
                ST29CV = ST29SV) %>%
  dplyr::filter(OVERALLQC == "Pass")

mri_data <- mri_data %>%
  dplyr::select(RID, EXAMDATE, LONIUID, IMAGEUID, STUDY, PTGENDER, age, ends_with("TA") | ends_with("CV"))

##########################################################################################################################
# ICV-adjustment
##########################################################################################################################
#making sure that the adjusted volumes are back in
volumes <- mri_data %>%
  dplyr::select(contains("CV"))

# creating list of variables that need to be ICV-adjusted
variables <- c()

for (roi_num in 1:length(names(volumes))){
  
  current_roi <- names(volumes)[roi_num]
  variables <- append(variables, current_roi)
  
}

# now doing the ICV adjustment
for (i in variables[-9]){
  #regressing out ICV
  lm_mod <- lm(mri_data[[i]] ~ ST10CV, data = mri_data, na.action=na.exclude)
  
  #creating a new variable for each i
  new_var_name <- paste("adjusted", i, sep = "_")
  mri_data[[new_var_name]] <- coef(lm_mod)[1] + resid(lm_mod)} #intercept + residuals, most often you just see people using the residuals as the adjusted value (see also "adjust" function from datawizard package)

#getting rid of the previous columns that are now replaced with ICV-adjusted values
mri_data <- mri_data %>%
  dplyr::select(-names(volumes))

##########################################################################################################################
# Hatmonization and pre-processing for harmonization
##########################################################################################################################

# getting diagnosis into the dataset by matching the closest diagnosis
mri_data <- merge(mri_data, diagnoses, all.x = TRUE) %>% # diagnoses processing in other processing step file
  dplyr::mutate(date_diff = abs(lubridate::time_length(difftime(EXAMDATE, DX.DATE), "years"))) %>%
  dplyr::group_by(RID, EXAMDATE, LONIUID, IMAGEUID, STUDY, PTGENDER, age) %>%
  dplyr::mutate(min_match_date = min(date_diff)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(min_match_date == date_diff) %>%
  dplyr::distinct()

# getting mri data split into features and covariates
all_features <- mri_data[, c(8:150)]
all_features <- all_features %>%
  dplyr::select_if(~mean(is.na(.)) < 0.2)
all_covariates <- mri_data[, c(1:7, 155)]

# now creating dataset for harmonization
all_data_combined_CN <- cbind(all_features, all_covariates) %>%
  dplyr::filter(diags == "CN")

all_features_CN <- all_data_combined_CN[, 1:138]
all_covariates_CN <- all_data_combined_CN[, c(143:145)]
extra_covariates_CN <- all_data_combined_CN[, c(139:142, 146)]

all_covariates_CN$STUDY <- as.factor(all_covariates_CN$STUDY)
all_covariates_CN$age <- as.numeric(all_covariates_CN$age)
all_covariates_CN$PTGENDER <- as.factor(all_covariates_CN$PTGENDER)

#harmonization
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

all_features_MCI_AD <- all_data_combined_MCI_AD[, 1:138] 
all_covariates_MCI_AD <- all_data_combined_MCI_AD[, c(143:145)]  
extra_covariates_MCI_AD <- all_data_combined_MCI_AD[, c(139:142, 146)] 

all_covariates_MCI_AD$STUDY <- as.factor(all_covariates_MCI_AD$STUDY)
all_covariates_MCI_AD$age <- as.numeric(all_covariates_MCI_AD$age)
all_covariates_MCI_AD$PTGENDER <- as.factor(all_covariates_MCI_AD$PTGENDER)

MCI_AD_harmonized <- t(ApplyHarm(feature.data   = all_features_MCI_AD,
                                 covariate.data = all_covariates_MCI_AD,
                                 comgam.out     = CN_data_harmonized))

MCI_AD_harmonized <- as.data.frame(MCI_AD_harmonized)
MCI_AD_harmonized <- cbind(extra_covariates_MCI_AD, all_covariates_MCI_AD, MCI_AD_harmonized)

harmonized_data_from_CN <- rbind(CN_harmonized_data, MCI_AD_harmonized) #%>%

##########################################################################################################################
# Age-adjustment
##########################################################################################################################
# creating list of variables that need to be ICV-adjusted
roi_variables <- c()

rois <- harmonized_data_from_CN %>%
  dplyr::select(ends_with("SV") | ends_with("TA") | ends_with("CV") | ends_with("SA"))

for (roi_num in 1:length(names(rois))){
  
  current_roi <- names(rois)[roi_num]
  roi_variables <- append(roi_variables, current_roi)
  
}

mri_plot_data <- harmonized_data_from_CN

mri_data_cn <- mri_plot_data %>%
  dplyr::filter(diags == "CN") %>%
  dplyr::mutate(RID = as.character(RID))

#getting amyloid negative cases
amyloid_negs <- amyloid_pet %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(AmyloidPosPET == 0) %>%
  dplyr::mutate(first_a_neg_date_pet = min(EXAMDATE_pet),
                RID = as.character(RID)) %>%
  dplyr::select(RID, first_a_neg_date_pet) %>%
  dplyr::ungroup() %>%
  dplyr::distinct()

mri_data_cn <- mri_data_cn %>%
  dplyr::left_join(amyloid_negs) %>%
  dplyr::mutate(diff_time = abs(lubridate::time_length(difftime(EXAMDATE, first_a_neg_date_pet), "years"))) %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(diff_time == min(diff_time)) %>%
  dplyr::ungroup()
mri_data_cn <- mri_data_cn %>%
  dplyr::group_by(RID) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()


# now doing the age adjustment
for (i in roi_variables){
  #regressing out age in cn
  lm_mod_age <- lm(mri_data_cn[[i]] ~ age + poly(age, 2, raw = TRUE)[,"2"], data = mri_data_cn, na.action=na.exclude)
  
  #creating a new variable for each i
  predictions_age <- predict(lm_mod_age, mri_plot_data)
  mri_plot_data[[i]] <- mri_plot_data[[i]] - predictions_age
  mri_plot_data[[i]] <- coef(lm_mod_age)[1] + mri_plot_data[i] } #intercept + residuals, most often you just see people using the residuals as the adjusted value (see also "adjust" function from datawizard package)

##########################################################################################################################
# calculating meta-ROI from regions according to Jack Jr. et al
##########################################################################################################################

mri_plot_data <- mri_plot_data %>%  
  dplyr::mutate(meta_ROI_left = (ST24TA[[1]] + ST32TA[[1]] + ST40TA[[1]] + ST26TA[[1]])/4,  # meta_ROI is mean cortical thickness in the following individual ROIs: entorhinal, inferior temporal, middle temporal, and fusiform.
                meta_ROI_right = (ST83TA[[1]] + ST91TA[[1]] + ST99TA[[1]] + ST85TA[[1]])/4)
