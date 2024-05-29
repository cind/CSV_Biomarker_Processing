# processing steps for ICV-adjusted, age-adjusted, harmonized ROI values (also creating cortical thickness meta-ROI)
library(dplyr)
# ComBat harmonization
source('~/Documents/Projects/MixedPathologies/ComGamPackage-master/ComGamFunctionHelpers.R')
source('~/Documents/Projects/MixedPathologies/ComGamPackage-master/ComGamHarmFunction.R')

#first getting demographics (getting information to calculate age)
dem <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\demographics.csv") %>%
  dplyr::mutate(day = 1,
                birthdate = as.character(paste(day, PTDOB, sep = "/")))
dem <- dem %>%
  dplyr::mutate(birthdate = as.Date(birthdate, format = "%d/%m/%Y")) %>%
  dplyr::select(RID, PTGENDER, PTEDUCAT, birthdate) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(birthdate))
dem <- merge(dem, apoeres, all = TRUE) %>% #(apoe processing- steps are in other file)
  dplyr::mutate(apoe = case_when(apoe == "E2/E3" | apoe == "E2/E2" ~ "E2",
                                 apoe == "E3/E3" ~ "E3",
                                 apoe == "E3/E4" | apoe == "E4/E4" ~ "E4"),
                RID = as.character(RID))

# calling in all freesurfer
adni1_1.5 <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\freesurfer_4.3_ADNI1_protocol_1.5T.csv") %>%
  dplyr::mutate(STUDY = "ADNI1",
                Field_Strength = "1.5T")
adni2_3 <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\freesurfer_5.1_ADNI2GO_protocol_3T.csv") %>%
  dplyr::mutate(STUDY = "ADNI2",
                Field_Strength = "3T")
adni3_3 <- read.csv("C:\\Work Folder\\paper data longitudinal phases\\freesurfer_6.0_ADNI3_protocol_3T.csv") %>%
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
  dplyr::rename(ST88CV = ST88SV, #these are labelled SV but are actually CV
                ST29CV = ST29SV) %>%
  dplyr::filter(OVERALLQC == "Pass")

mri_data <- mri_data[,colSums(is.na(mri_data))<nrow(mri_data)]

#putting aside cortical volumes to ICV-adjust
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
  dplyr::select(-names(volumes)) %>%
  dplyr::filter(!STUDY == "ADNI1") #currently getting rid of ADNI1 values due to weird age association issues with ADNI1 freesurfer data
