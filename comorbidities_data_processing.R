library(ADNIMERGE)
library(tidyverse)
library(transplantr) ## used to calculate eGFR

## temp code to get DXs in as same type 
adnimerge_unlabelled<-sjlabelled::unlabel(ADNIMERGE::adnimerge %>% 
                                            dplyr::mutate(ORIGPROT=as.character(ORIGPROT),
                                                          DX=as.character(DX),
                                                          DX.bl=as.character(DX.bl)))

## would recommend downloading fresh copies of any listed CSVs
## init_health is health reported at initial visit
init_health<-readr::read_delim("~/Projects/Comorbidities Project/INITHEALTH.csv")
## rec_hist contains updates in subsequent visits
rec_hist<- readr::read_delim("~/Projects/Comorbidities Project/RECMHIST.csv")

## creates dataframes with cleaned date formats for initial health and recent medical history

init_health <- init_health %>% dplyr::mutate(
  onset=case_when(
    (is.na(IHDTONSET)) ~ USERDATE,
    (lubridate::month(IHDTONSET)==7 & (lubridate::day(IHDTONSET)==2|lubridate::day(IHDTONSET)==1)) ~ lubridate::ceiling_date(init_health$IHDTONSET,unit="years"), ## Initial health incorrectly has most records as having onset/cease dates of 7/1 or 7/2
    (!is.na(IHDTONSET)) ~ lubridate::ymd(IHDTONSET)
  ),
  cease=case_when(
    (lubridate::month(IHCEASE)==7 & (lubridate::day(IHCEASE)==2|lubridate::day(IHCEASE)==1)) ~ lubridate::floor_date(init_health$IHCEASE,unit="years"), ## Initial health incorrectly has most records as having onset/cease dates of 7/1 or 7/2
    (!is.na(IHCEASE)) ~ lubridate::ymd(IHCEASE)
  )
)

rec_hist$MHDTONSET<-ifelse(rec_hist$MHDTONSET=="x1",NA,rec_hist$MHDTONSET)
rec_hist$MHDTONSET<-ifelse(rec_hist$MHDTONSET=="x4",NA,rec_hist$MHDTONSET)

rec_hist$onset_year<-lubridate::year(strptime(stringr::str_sub(start=-4,rec_hist$MHDTONSET),format="%Y"))
rec_hist$onset_month<-lubridate::month(strptime(paste(stringr::str_sub(end=2,rec_hist$MHDTONSET),"01",stringr::str_sub(start=-4,rec_hist$MHDTONSET),sep="-"),format="%m-%d-%Y"))
rec_hist$onset_day<-lubridate::day(strptime(stringr::str_sub(start=4,end=5,rec_hist$MHDTONSET),format="%d"))

rec_hist <- rec_hist %>% dplyr::mutate(onset=case_when(
  (is.na(onset_year)) ~ USERDATE,
  (is.na(onset_month)) ~ lubridate::ceiling_date(lubridate::ymd(paste(onset_year,"01","01",sep="-")),unit="years"),
  (is.na(onset_day)) ~ lubridate::ceiling_date(lubridate::ymd(paste(onset_year,onset_month,"01",sep="-")),unit="months"),
  (!is.na(onset_day)) ~ lubridate::mdy(MHDTONSET)
))

rec_hist_scraper<- rec_hist %>%
  dplyr::filter(MHCUR == 1) %>%
  dplyr::select(RID, VISCODE2, onset, MHDESC,USERDATE) %>%
  dplyr::rename(viscode=VISCODE2,desc=MHDESC)

init_health_scraper <- init_health %>%
  dplyr::select(RID, VISCODE2, onset, cease, IHDESC,USERDATE) %>%
  dplyr::rename(viscode=VISCODE2,desc=IHDESC)

merged_health_scraper<- dplyr::bind_rows(rec_hist_scraper,init_health_scraper)

## AB42/40 data load-in
## Removed records w/ BATCH_N = "Batch #00 Trial" from Bateman because they were duplicates
bateman_plasma_nov_22 <-readr::read_delim("~/Downloads/batemanlab_20221118.csv") 
colnames(bateman_plasma_nov_22) <- str_to_upper(colnames(bateman_plasma_nov_22))
bateman_plasma_nov_22 <- bateman_plasma_nov_22 %>%
  dplyr::rename(abeta_ratio=ABETA_4240_STANDARDIZED) %>%
  ## some plasma records don't have RIDs and some have them = 999999; these have to be thrown out because they can't be matched to other tables
  dplyr::filter(!is.na(RID))
bateman_plasma_nov_22$EXAMDATE<-lubridate::mdy(bateman_plasma_nov_22$EXAMDATE)
bateman_plasma_nov_22$MS_RUN_DATE<-lubridate::ymd(bateman_plasma_nov_22$MS_RUN_DATE)

bateman_plasma <-
  readr::read_delim("~/Projects/Comorbidities Project/batemanlab_20190621.csv") %>% dplyr::rename(abeta_ratio = RATIO_ABETA42_40_BY_ISTD_TOUSE) %>% dplyr::mutate(origin="Bateman") %>% dplyr::filter(BATCH_N!="Batch #00 Trial")

deming_df_old_bateman<- dplyr::left_join(bateman_plasma_nov_22 %>% dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a"),
                                         bateman_plasma %>% dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a"),
                                         by=c("RID","VISCODE2"),suffix=c(".new",".old"))
deming_df_old_bateman <- deming_df_old_bateman %>% dplyr::filter(RID!=999999) %>% tidyr::drop_na(abeta_ratio.new,abeta_ratio.old) %>% dplyr::mutate(diff=abeta_ratio.new-abeta_ratio.old)
deming_reg_old_bateman <-mcr::mcreg(y=deming_df_old_bateman$abeta_ratio.new,x=deming_df_old_bateman$abeta_ratio.old,method.reg = "Deming")

## coefficients for adjusted abeta ratio were calculated using Deming regression predicting Bateman 2022 values with Bateman 2019 values
bateman_plasma <- bateman_plasma %>% dplyr::rename(abeta_ratio_unadj=abeta_ratio) %>% dplyr::mutate(abeta_ratio=deming_reg_old_bateman@para[1,1]+(deming_reg_old_bateman@para[2,1]*abeta_ratio_unadj))

bateman_merged <-dplyr::bind_rows(bateman_plasma_nov_22 %>% dplyr::mutate(origin="New Bateman"),bateman_plasma %>% dplyr::mutate(origin="Old Bateman"))
bateman_merged <- bateman_merged %>% dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a",RID!=999999)

washu_plasma <-
  readr::read_delim("~/Projects/Comorbidities Project/PLASMA_ABETA_PROJECT_WASH_U_11_05_21.csv")
washu_plasma <-
  washu_plasma %>% dplyr::rename(abeta_ratio = STANDARDIZED_PLASMAAB4240,abeta_ratio_non_std = PLASMAAB4240) %>% dplyr::mutate(origin="WashU") %>% 
  dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a")

deming_df_washu<- dplyr::left_join(bateman_plasma_nov_22 %>% dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a"),
                                   washu_plasma %>% dplyr::filter(QC_STATUS == "Passed", (INSTRUMENT == "OTF Lumos"|INSTRUMENT == "Lumos"), INJECTION == "a"),
                                   by=c("RID","VISCODE2"),suffix=c(".new",".washu"))
deming_df_washu <- deming_df_washu %>% dplyr::filter(RID!=999999) %>% tidyr::drop_na(abeta_ratio.new,abeta_ratio.washu) %>% dplyr::mutate(diff=abeta_ratio.new-abeta_ratio.washu)
deming_reg_washu <-mcr::mcreg(y=deming_df_washu$abeta_ratio.new,x=deming_df_washu$abeta_ratio.washu,method.reg = "Deming")

## coefficients for adjusted abeta ratio were calculated using Deming regression predicting Bateman 2022 values with WashU 2021 values
washu_plasma <- washu_plasma %>% dplyr::rename(abeta_ratio_unadj=abeta_ratio) %>% dplyr::mutate(abeta_ratio=deming_reg_washu@para[1,1]+(deming_reg_washu@para[2,1]*abeta_ratio_unadj))

plasma_merged <- dplyr::bind_rows(
  washu_plasma %>% dplyr::mutate(MS_RUN_DATE=lubridate::ymd(MS_RUN_DATE)),
  bateman_merged %>% dplyr::filter(RID != 999999) %>% dplyr::mutate(VISCODE=VISCODE2))

plasma_merged <- plasma_merged %>% dplyr::mutate(is_qc_record=((RID %in% deming_df_old_bateman$RID)|(RID %in% deming_df_washu$RID)))
plasma_merged <- plasma_merged %>% dplyr::filter(is_qc_record==FALSE|origin=="New Bateman")
plasma_merged <- plasma_merged %>% dplyr::select(RID,VISCODE,VISCODE2,EXAMDATE,origin,abeta_ratio,abeta_ratio_unadj,MS_RUN_DATE,is_qc_record)
plasma_merged <- plasma_merged %>% dplyr::mutate(sample_age = MS_RUN_DATE - EXAMDATE)

## Plasma ptau-181 load-in
## Longitudinal data
## Had some subjects with two measurements at a given visit; averaged two measurements for those subjects
plasma_ptau <- 
  readr::read_delim("~/Projects/Comorbidities Project/UGOTPTAU181_06_18_20.csv") %>% dplyr::select(RID,VISCODE2,EXAMDATE,PLASMAPTAU181) %>% dplyr::rename(VISCODE=VISCODE2,ptau_181=PLASMAPTAU181)
plasma_ptau <- plasma_ptau %>% dplyr::group_by(RID,EXAMDATE) %>% dplyr::mutate(ptau_181=mean(ptau_181)) %>% dplyr::distinct_at(vars(RID,EXAMDATE),.keep_all = TRUE) 

plasma_ptau <- plasma_ptau %>% dplyr::arrange(RID,EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::mutate(ptau_change_from_bl=ptau_181-ptau_181[1L],baseline=ptau_181[1L],time=EXAMDATE-EXAMDATE[1L]) %>% dplyr::ungroup()

plasma_ptau %>% dplyr::group_by(RID) %>% dplyr::summarize(count=n()) %>% dplyr::group_by(count) %>% dplyr::summarize(metacount=n())

## Plasma NfL load-in
## Had some subjects with two measurements at a given visit; averaged two measurements for those subjects
plasma_nfl <-
  readr::read_delim("~/Projects/Comorbidities Project/ADNI_BLENNOWPLASMANFLLONG_10_03_18.csv") %>% dplyr::select(RID,VISCODE2,EXAMDATE,PLASMA_NFL) %>% dplyr::rename(VISCODE=VISCODE2,nfl=PLASMA_NFL)
plasma_nfl <- plasma_nfl %>% dplyr::group_by(RID,EXAMDATE) %>% dplyr::mutate(nfl=mean(nfl)) %>% dplyr::distinct_at(vars(RID,EXAMDATE),.keep_all = TRUE) 

plasma_nfl <- plasma_nfl %>% dplyr::arrange(RID,EXAMDATE) %>% dplyr::group_by(RID) %>% dplyr::mutate(nfl_change_from_bl=nfl-nfl[1L],baseline=nfl[1L],time=EXAMDATE-EXAMDATE[1L]) %>% dplyr::ungroup()

plasma_nfl %>% dplyr::group_by(RID) %>% dplyr::summarize(count=n()) %>% dplyr::group_by(count) %>% dplyr::summarize(metacount=n())

## PET data load-in
## Adapted from Adam's code
av45 <- readr::read_delim("~/Projects/Comorbidities Project/UCBERKELEYAV45_04_26_22.csv") %>% 
  dplyr::select(RID,EXAMDATE,VISCODE2,SUMMARYSUVR_WHOLECEREBNORM,SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF) %>% 
  dplyr::rename(suvr_summary=SUMMARYSUVR_WHOLECEREBNORM,AmyloidPos=SUMMARYSUVR_WHOLECEREBNORM_1.11CUTOFF) %>%
  dplyr::mutate(Centiloid = (188.22 * suvr_summary) - 189.16,tracer="av45") ## reflects most recent values for AV45 from "Instructions for converting ADNI processing results to Centiloids (PDF)" on LONI
fbb <- readr::read_delim("~/Projects/Comorbidities Project/UCBERKELEYFBB_04_26_22.csv") %>% 
  dplyr::select(RID,EXAMDATE,VISCODE2,SUMMARYSUVR_WHOLECEREBNORM,SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF) %>% 
  dplyr::rename(suvr_summary=SUMMARYSUVR_WHOLECEREBNORM,AmyloidPos=SUMMARYSUVR_WHOLECEREBNORM_1.08CUTOFF) %>%
  dplyr::mutate(Centiloid = (157.15 * suvr_summary) - 151.87,tracer="fbb") ## reflects most recent values for FBB from "Instructions for converting ADNI processing results to Centiloids (PDF)" on LONI
amyloid_pet<-dplyr::bind_rows(av45,fbb) %>% dplyr::mutate(tracer=as.factor(tracer)) %>% dplyr::rename(VISCODE=VISCODE2)

## CSF data load-in
upenn_mk_12 <-readr::read_delim("~/Projects/Comorbidities Project/UPENNBIOMK12_01_04_21.csv") %>% dplyr::mutate(origin="MK 12")
upenn_mk_10 <-readr::read_delim("~/Downloads/UPENNBIOMK10_07_29_19.csv") %>% dplyr::rename(EXAMDATE=DRAWDATE,ABETA=ABETA42) %>% dplyr::mutate(origin="MK 10")

upenn_mk_9 <-readr::read_delim("~/Projects/Comorbidities Project/UPENNBIOMK9_04_19_17.csv")
upenn_mk_9 <- upenn_mk_9 %>% dplyr::mutate(ABETA=as.numeric(ABETA),TAU=as.numeric(TAU),PTAU=as.numeric(PTAU),origin="MK 9")
upenn_mk_9$ABETA <- ifelse(is.na(upenn_mk_9$ABETA),readr::parse_number(upenn_mk_9$COMMENT),upenn_mk_9$ABETA)

upenn_merged_csf_biomarkers <- dplyr::bind_rows(upenn_mk_12,upenn_mk_10,upenn_mk_9)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::select(-VISCODE) %>% dplyr::rename(VISCODE=VISCODE2)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::mutate(ABETA=round(ABETA),ptau_pos=PTAU>24)
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::arrange(RUNDATE) %>% dplyr::distinct_at(vars(RID,VISCODE),.keep_all=TRUE)

## split ABETA values here; values <200 or >1700 are usable for status assignment but not for numerical analysis
## don't need to split CSF PTAU - all values are within technical limits
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::mutate(ABETA_for_status=ABETA,
  ABETA=case_when(
  ABETA>=200 & ABETA<=1700 ~ ABETA
))
upenn_merged_csf_biomarkers <- upenn_merged_csf_biomarkers %>% dplyr::group_by(RID,EXAMDATE) %>% dplyr::filter(RUNDATE==max(RUNDATE)) %>% dplyr::ungroup()

##U-p53 load-in

alzosure <-readr::read_delim("~/Downloads/DIAMEM_V2_08_30_22.csv")
alzosure <- alzosure %>% dplyr::select(RID,VISCODE2,EXAMDATE,AZ284) %>% dplyr::rename(VISCODE=VISCODE2) %>%
  dplyr::group_by(RID,EXAMDATE) %>% dplyr::mutate(up53=mean(AZ284))
alzosure <- alzosure %>% dplyr::distinct_at(vars(RID,EXAMDATE),.keep_all=TRUE) %>% dplyr::select(-AZ284)

## Hachinski and vitals
## Hachinski refers to a questionnare for Hachinski's ischemic score, which uses information on known vascular comorbidities to differentiate types of dementia
adni_mod_hach <- ADNIMERGE::modhach %>% dplyr::filter(VISCODE!="nv")
adni_vitals <- ADNIMERGE::vitals %>% dplyr::filter(VISCODE!="nv")

## Sleep load-in
adni_npi<-ADNIMERGE::npi
adni_sleep_npi<-adni_npi %>% dplyr::select(RID,USERDATE,NPIK,NPIKTOT) %>% 
  dplyr::mutate(sleep_quest=case_when(
    is.na(NPIK) ~ 0,
    NPIK=="No" ~ 0,
    NPIK=="Yes" ~ 1)) %>% 
  dplyr::rename(sleep_score=NPIKTOT,sleep_date=USERDATE)

## Demographics
## I use ADNIMERGE here because the demographic categories are well-labeled in ADNIMERGE but not LONI
## Still need web data because ADNIMERGE doesn't have month and year of birth
adni_merge_demog <- ADNIMERGE::ptdemog
adni_web_demog <-
  readr::read_delim("~/Projects/Comorbidities Project/PTDEMOG.csv")
adni_merge_demog_uniques <- adni_merge_demog %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adni_web_demog_uniques <- adni_web_demog %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adni_joined_demog_uniques <-
  dplyr::left_join(
    adni_merge_demog_uniques,
    adni_web_demog_uniques %>%
      dplyr::select(PTDOBYY, PTDOBMM, RID),
    by = "RID"
  )

adni_joined_demog_uniques_reduced <-
  adni_joined_demog_uniques %>% dplyr::select(RID,PTGENDER,PTETHCAT,PTRACCAT,PTDOBYY,PTDOBMM,PTEDUCAT,ORIGPROT)

## Blood lab data load-in
## This is somewhat involved because I needed to calculate eGFR for CKD, which depends on age, gender, and race
## adni_lab_data <- ADNIMERGE::labdata %>% dplyr::filter(VISCODE!="nv")

adni_lab_data <- readr::read_delim("~/Projects/Comorbidities Project/LABDATA.csv") %>% dplyr::filter(VISCODE!="nv") %>% dplyr::mutate(across(everything(), ~replace(., . ==  -1 , NA)))

adni_lab_data_demog_reduced <-
  dplyr::left_join(adni_lab_data,adni_joined_demog_uniques_reduced,by = "RID",suffix = c(".labs", ".demos"))

adni_lab_data_demog_reduced <-
  adni_lab_data_demog_reduced %>% dplyr::mutate(
    glucose = as.numeric(RCT11),
    creatinine = as.numeric(RCT392),
    cholesterol = as.numeric(RCT20),
    triglycerides = as.numeric(RCT19),
    ast = as.numeric(RCT5),
    alt = as.numeric(RCT4),
    platelets = as.numeric(HMT13)
  )

adni_lab_data_demog_reduced <-adni_lab_data_demog_reduced %>% dplyr::filter(!(is.na(creatinine)|is.na(glucose)|is.na(cholesterol)|is.na(triglycerides)))

adni_lab_data_demog_reduced$age_at_lab <-round(as.numeric(lubridate::year(adni_lab_data_demog_reduced$EXAMDATE)) -
                                                 adni_lab_data_demog_reduced$PTDOBYY +
                                                 ((as.numeric(lubridate::month(adni_lab_data_demog_reduced$EXAMDATE)) - 
                                                     adni_lab_data_demog_reduced$PTDOBMM) / 12), digits=1)

adni_lab_data_demog_reduced$race_for_egfr <-
  ifelse(
    adni_lab_data_demog_reduced$PTRACCAT == "Black or African American",
    "black",
    "non-black"
  )

adni_lab_data_demog_reduced$gender_for_egfr <-
  ifelse(
    adni_lab_data_demog_reduced$PTGENDER == "Male",
    "M",
    "F"
  )

adni_lab_data_demog_reduced$eGFR <-
  transplantr::ckd_epi(
    adni_lab_data_demog_reduced$creatinine,
    adni_lab_data_demog_reduced$age_at_lab,
    adni_lab_data_demog_reduced$gender_for_egfr,
    adni_lab_data_demog_reduced$race_for_egfr,
    units = "US"
  )

adni_lab_data_demog_reduced <- adni_lab_data_demog_reduced %>% dplyr::mutate(
  fibrosis_4_score = ((age_at_lab*ast)/(platelets*sqrt(alt))),
  ast_alt_ratio = ast/alt
)

## Text strings for NAs are to avoid missing labs causing an error in regression analyses
## Should consider improving by either imputing lab values (if subjects have some lab results but not all) or just NAing subjects without labs (if some subjects have no lab results)

adni_lab_data_demog_reduced <- adni_lab_data_demog_reduced %>% dplyr::mutate(eGFR_cat = with(.,case_when(
  (eGFR>90)~"Greater than 90",
  (eGFR>60 & eGFR<=90)~"61-90",
  (eGFR>45 & eGFR <=60)~"45-60",
  (eGFR<=45)~"Less than 45",
  (is.na(eGFR))~"No eGFR Available")))

adni_lab_data_demog_reduced <- adni_lab_data_demog_reduced %>%
  mutate(., diabetes_cat = with(
    .,
    case_when(
      (glucose >= 125) ~ 'Diabetes',
      (glucose >= 100 & glucose < 125) ~ 'Impaired Glycemia',
      (glucose < 100) ~ 'Normoglycemic',
      is.na(glucose) ~ 'Failed Draw',
    )
  )) %>%
  mutate(., ckd_cat = with(
    .,
    case_when(
      (eGFR <= 60) ~ 'Suspected CKD',
      (eGFR > 60) ~ 'No Suspected CKD',
      is.na(eGFR) ~ 'Could Not Calculate eGFR',
    )
  )) %>%
  mutate(., dyslipidemia_cat = with(
    .,
    case_when(
      (cholesterol >= 240) ~ 'Hypercholesteremia',
      (cholesterol >= 200) ~ 'Borderline High',
      (cholesterol < 200) ~ 'Normal Cholesterol',
      is.na(cholesterol) ~ 'Failed Draw',
    )),
    hypertriglyceridemia_cat = 
      case_when(
        (triglycerides >= 150) ~ 'Hypertriglyceridemia',
        (triglycerides < 150) ~ 'Normal Triglycerides',
        is.na(cholesterol) ~ 'Failed Draw',
      ),
    fibrosis_cat = 
      case_when(
        (fibrosis_4_score > 3.25) ~ 'Likely Fibrosis',
        (fibrosis_4_score < 3.25 & fibrosis_4_score >= 1.45) ~ 'No Predictive Value',
        (fibrosis_4_score < 1.45) ~ 'Unlikely Fibrosis',
        is.na(fibrosis_4_score) ~ 'Failed Draw'
      ),
    alcoholic_liver_disease_cat = 
      case_when(
        (ast_alt_ratio > 2) ~ 'Likely ALD',
        (ast_alt_ratio > 1) ~ 'Possible Cirrhosis',
        (ast_alt_ratio <= 1) ~ 'Unlikely ALD',
        is.na(ast_alt_ratio) ~ 'Failed Draw'
      ),
    elevated_alt = 
  case_when(
    ((PTGENDER=="Male" & alt<=42)|(PTGENDER=="Female" & alt<=30)) ~ "ALT Below ULN",
    ((PTGENDER=="Male" & alt<=42*2)|(PTGENDER=="Female" & alt<=30*2)) ~ "Borderline ALT Elevation",
    ((PTGENDER=="Male" & alt<=42*5)|(PTGENDER=="Female" & alt<=30*5)) ~ "Mild ALT Elevation",
    ((PTGENDER=="Male" & alt<=42*15)|(PTGENDER=="Female" & alt<=30*15)) ~ "Moderate ALT Elevation",
    ((PTGENDER=="Male" & alt>42*15)|(PTGENDER=="Female" & alt>30*15)) ~ "Severe ALT Elevation"
  )) %>%
  mutate(., ckd_cat_alt = with(
    .,
    case_when(
      (eGFR <= 40) ~ 'Suspected CKD',
      (eGFR > 40) ~ 'No Suspected CKD',
      is.na(eGFR) ~ 'Could Not Calculate eGFR',
    )
  )) %>%
  dplyr::select(RID,age_at_lab,glucose,diabetes_cat,creatinine,eGFR,ckd_cat,ckd_cat_alt,cholesterol,dyslipidemia_cat,race_for_egfr,gender_for_egfr,eGFR_cat,hypertriglyceridemia_cat,fibrosis_4_score,fibrosis_cat,ast_alt_ratio,alcoholic_liver_disease_cat,ast,alt,elevated_alt,EXAMDATE
  ) %>% dplyr::rename(lab_date=EXAMDATE)

adni_lab_data_use<- adni_lab_data_demog_reduced %>%
  dplyr::distinct_at(vars(RID), .keep_all = TRUE)

## Hachinski & Vitals load-in

adni_mod_hach <- ADNIMERGE::modhach
adni_mod_hach$htn_hach <-
  ifelse(adni_mod_hach$HMHYPERT == "Present - 1 point", 1, 0)

adni_vitals <- ADNIMERGE::vitals

adni_vitals <-
  adni_vitals %>% dplyr::mutate(
    diastolic_bp = as.numeric(VSBPDIA),
    systolic_bp = as.numeric(VSBPSYS),
    vitals_date = USERDATE) 

adni_vitals <- adni_vitals %>% dplyr::filter(!(is.na(systolic_bp)|is.na(diastolic_bp)))

adni_vitals$htn_vitals <-
  ifelse(adni_vitals$diastolic_bp >= 90 |
           adni_vitals$systolic_bp >= 140,
         1,
         0)

test<-adni_vitals %>% dplyr::distinct_at(vars(USERDATE,RID))

## Diagnostic info load-in
adni_diagnoses<-adnimerge_unlabelled %>% 
  dplyr::mutate(RID=as.numeric(RID)) %>%
  dplyr::select(RID,DX,DX.bl,EXAMDATE,VISCODE) %>% dplyr::rename(dx_date=EXAMDATE) %>%
  dplyr::mutate(DX.baseline = case_when(
    (DX.bl %in% c("EMCI","LMCI")) ~ "MCI",
    (DX.bl %in% c("SMC","CN")) ~ "CN",
    (DX.bl == "Dementia" ~ "AD")
  )) %>%
  dplyr::filter(!is.na(DX))

## APOE info load-in
adni_apoe<-ADNIMERGE::adnimerge %>% 
  dplyr::mutate(RID=as.numeric(RID)) %>%
  dplyr::select(RID,APOE4) %>% dplyr::mutate(apoe_status=(APOE4>0)) %>%
  dplyr::filter(!is.na(APOE4)) %>%
  dplyr::distinct_at(vars(RID),.keep_all = TRUE)

## merged, all biomarker dataset
all_biomarker_data<-dplyr::full_join(plasma_merged,plasma_nfl %>% dplyr::select(RID,nfl),by="RID")
all_biomarker_data<- all_biomarker_data %>% dplyr::full_join(plasma_ptau %>% dplyr::select(RID,ptau_181),by="RID")
all_biomarker_data<- all_biomarker_data %>% dplyr::full_join(amyloid_pet %>% dplyr::select(RID,AmyloidPos),by="RID")
all_biomarker_data<- all_biomarker_data %>% dplyr::full_join(upenn_merged_csf_biomarkers %>% dplyr::select(RID,ptau_pos),by="RID")
## all_biomarker_data<- all_biomarker_data %>% dplyr::full_join(alzosure %>% dplyr::select(RID,up53),by="RID")

all_biomarker_data<- all_biomarker_data %>% dplyr::mutate(has_abeta=!is.na(abeta_ratio),has_ptau=!is.na(ptau_181),has_nfl=!is.na(nfl),
                                                          has_pet=!is.na(AmyloidPos),has_csf=!is.na(ptau_pos))
all_biomarker_data <- all_biomarker_data %>% dplyr::filter(has_abeta|has_ptau|has_nfl) %>% dplyr::filter(has_pet|has_csf)
all_biomarker_data<- all_biomarker_data %>% dplyr::distinct_at(vars(RID),.keep_all = TRUE)
all_biomarker_data <- all_biomarker_data %>% dplyr::left_join(adni_vitals %>% dplyr::filter(VISCODE=="sc"|VISCODE=="f") %>% dplyr::select(RID,USERDATE),by="RID")
all_biomarker_data <- all_biomarker_data %>% dplyr::mutate(EXAMDATE=USERDATE)
