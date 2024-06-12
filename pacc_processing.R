## table name: Neuropsychological Battery [ADNI1,GO,2,3]	

library(tidyverse)

source("~/ADAS_processing.R") # change to name and location of your ADAS processing file
source("~/mmse_processing.R") # change to name and location of your MMSE processing file
source("~/Downloads/diagnoses_processing.R") # change to name and location of your diagnosis processing file

neurobat_source <- readr::read_delim("~/Projects/Amprion Project/Source Data/NEUROBAT_22May2024.csv") # change to name and location of Neuropsychological Battery

## debugged version of the pacc function from ADNIMERGE
## dd is a dataset that contains cols DX.bl, ADASQ4, LDELTOTAL, DIGITSCOR, TRABSCOR, MMSE

debugged_pacc_fn<-function (dd, keepComponents = FALSE) 
{
  require(dplyr)
  require(tidyr)
  holdNames <- colnames(dd)
  dd <- mutate(dd, log.TRABSCOR = log(TRABSCOR + 1))
  dd$log.TRABSCOR <- ifelse(is.finite(dd$log.TRABSCOR),dd$log.TRABSCOR,NA)
  bl.summary <- dplyr::filter(dd, VISCODE == "bl") %>% dplyr::select(DX.bl, 
                                                                     ADASQ4, LDELTOTAL, DIGITSCOR, log.TRABSCOR, MMSE) %>% 
    tidyr::gather(VARIABLE, SCORE, -DX.bl) %>% dplyr::filter(!is.na(SCORE)) %>% 
    dplyr::group_by(DX.bl, VARIABLE) %>% dplyr::summarize(N = n(), 
                                                          mean = mean(SCORE), sd = sd(SCORE)) 
  zscore <- function(x, var) {
    (x - dplyr::filter(bl.summary, DX.bl == "CN" & VARIABLE == var)$mean)/dplyr::filter(bl.summary, 
                                                                                        DX.bl == "CN" & VARIABLE == var)$sd
  }
  dd <- dplyr::mutate(dd, ADASQ4.z = -zscore(ADASQ4, "ADASQ4"), LDELTOTAL.z = zscore(LDELTOTAL, 
                                                                                     "LDELTOTAL"), DIGITSCOR.z = zscore(DIGITSCOR, "DIGITSCOR"), 
                      log.TRABSCOR.z = -zscore(log.TRABSCOR, "log.TRABSCOR"), 
                      MMSE.z = zscore(MMSE, "MMSE"))
  corTest <- cor(dplyr::select(dd, ADASQ4.z, LDELTOTAL.z, DIGITSCOR.z, 
                               log.TRABSCOR.z, MMSE.z), use = "pairwise.complete.obs")
  if (any(corTest < 0)) 
    stop("Some PACC z scores are negatively correlated!")
  compscore <- function(x, n.components = 4, n.missing = 2) {
    ifelse(sum(is.na(x)) > n.missing, NA, mean(x, na.rm = TRUE)) * 
      n.components
  }
  dd$mPACCdigit <- apply(dd[, c("ADASQ4.z", "LDELTOTAL.z", 
                                "DIGITSCOR.z", "MMSE.z")], 1, compscore)
  dd$mPACCtrailsB <- apply(dd[, c("ADASQ4.z", "LDELTOTAL.z", 
                                  "log.TRABSCOR.z", "MMSE.z")], 1, compscore)
  if (!keepComponents) 
    dd <- dd[, c(holdNames, "mPACCdigit", "mPACCtrailsB")]
  dd <- as.data.frame(dd) %>%
    dplyr::distinct_at(vars(RID,VISCODE,mPACCdigit,mPACCtrailsB))
}

## LDEL and DIGITSCOR/TRABSCOR are recorded at separate baseline/screening visits - below code aligns to one visit
neurobat_non_bl <- neurobat_source %>%
  dplyr::filter(!(VISCODE2 %in% c("bl","sc","f"))) %>%
  dplyr::select(RID,VISCODE2,LDELTOTAL,DIGITSCOR,TRABSCOR)

neurobat_ldel_bl <- neurobat_source %>%
  dplyr::filter(VISCODE2 %in% c("bl","sc","f")) %>%
  dplyr::select(RID,VISCODE2,LDELTOTAL) %>%
  dplyr::mutate(VISCODE2="bl") %>%
  dplyr::arrange(RID)

neurobat_scores_bl <- neurobat_source %>%
  dplyr::filter(VISCODE2 %in% c("bl","sc","f")) %>%
  dplyr::select(RID,VISCODE2,DIGITSCOR,TRABSCOR) %>%
  dplyr::mutate(VISCODE2="bl") %>%
  dplyr::arrange(RID)

## selects correct merged rows
neurobat_bl <- dplyr::full_join(neurobat_ldel_bl,neurobat_scores_bl,by=c("RID","VISCODE2")) %>%
  dplyr::mutate(total_na = is.na(DIGITSCOR)+is.na(TRABSCOR)+is.na(LDELTOTAL)) %>%
  dplyr::arrange(total_na) %>%
  dplyr::distinct_at(vars(RID,VISCODE2),.keep_all=TRUE)

## averages records for subjects with multiple batteries evaluated
neurobat <- dplyr::bind_rows(neurobat_bl,neurobat_non_bl) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adas_q4_for_pacc_1 <- adas_1_score %>%
  dplyr::select(RID,VISCODE2,Q4) %>%
  dplyr::rename(ADASQ4=Q4)

adas_q4_for_pacc_2_3_go <- adas_2_3_go %>%
  dplyr::select(RID,VISCODE2,Q4SCORE) %>%
  dplyr::rename(ADASQ4=Q4SCORE)

adas_q4_for_pacc <- dplyr::bind_rows(adas_q4_for_pacc_1,adas_q4_for_pacc_2_3_go) %>%
  dplyr::group_by(RID,VISCODE2) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

pacc_df <- dplyr::left_join(neurobat,
                            adas_q4_for_pacc,
                            by=c("RID","VISCODE2")) %>%
  dplyr::left_join(.,
                   mmse %>% dplyr::select(RID,VISCODE2,MMSE),
                   by=c("RID","VISCODE2")) %>%
  dplyr::left_join(.,
                   adni_diagnoses %>%
                     dplyr::rename(DX.bl=DX) %>%
                     dplyr::filter(VISCODE2=="bl") %>% 
                     dplyr::select(RID,DX.bl),
                   by=c("RID")) %>%
  dplyr::mutate(log.TRABSCOR=log(TRABSCOR+1)) %>%
  dplyr::select(RID,VISCODE2,DX.bl,ADASQ4, LDELTOTAL, DIGITSCOR,TRABSCOR,MMSE) %>%
  tidyr::drop_na(DX.bl) %>%
  dplyr::mutate(DX.bl = 
                  case_when(DX.bl == "CU" ~ "CN",
                            TRUE ~ DX.bl)) %>%
  dplyr::rename(VISCODE=VISCODE2)

pacc <- debugged_pacc_fn(dd=data.frame(pacc_df)) %>%
  dplyr::group_by(RID,VISCODE) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

