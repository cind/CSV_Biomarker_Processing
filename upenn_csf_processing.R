## table name: UPENN CSF Biomarkers Roche Elecsys [ADNI1,GO,2,3]

library(tidyverse)

upenn_csf_biomarkers <-readr::read_delim("~/Projects/Amprion Project/Source Data/UPENNBIOMK_ROCHE_ELECSYS_21May2024.csv") # edit to your source file name + location 

## keep only most recent replicate for replicated samples
upenn_csf_biomarkers <- upenn_csf_biomarkers %>% 
  dplyr::arrange(desc(RUNDATE)) %>% 
  dplyr::distinct_at(vars(RID,VISCODE2),.keep_all=TRUE)

## rename biomarkers with CSF tags
upenn_csf_biomarkers <- upenn_csf_biomarkers %>%
  dplyr::rename(ABETA40_csf=ABETA40,
                ABETA42_csf=ABETA42,
                TAU_csf=TAU,
                PTAU_csf=PTAU)

## create positivity status variables for CSF AB42 < 980, CSF p-tau181 > 24, and AD pathology: p-tau181/AB42 > 0.025
upenn_csf_biomarkers <- upenn_csf_biomarkers %>% dplyr::mutate(ABETA42_csf=round(ABETA42_csf),
                                                               ptau_pos_csf=PTAU_csf>24,
                                                               amyloid_pos_csf=ABETA42_csf<980,
                                                               ptau_ab_ratio_csf=(as.numeric(PTAU_csf)/as.numeric(ABETA42_csf)),
                                                               ad_pathology_pos_csf=ptau_ab_ratio_csf>0.025)

## convert all biomarker vars to numeric
upenn_csf_biomarkers <- upenn_csf_biomarkers %>%
  dplyr::mutate(across(.cols=c(ABETA42_csf,ABETA40_csf,TAU_csf,PTAU_csf,ptau_ab_ratio_csf),.fns=as.numeric))

## reassign records with VISCODE2 = "UNK" to correct VISCODE2s
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==89 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m162"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==232 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m36"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==726 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m24"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==790 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m18"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==1200 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m36"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==2373 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m120"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==4216 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m102"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==4393 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m96"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==4643 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m72"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==5140 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"m84"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==6161 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"bl"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==6880 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"bl"
upenn_csf_biomarkers$VISCODE2[upenn_csf_biomarkers$RID==6906 & upenn_csf_biomarkers$VISCODE2=="UNK"]<-"bl"
