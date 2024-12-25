## last updated 6/23/24 ZH

## SEARCH TABLE NAME BELOW IN LONI
## table name: Neuropsychological Battery [ADNI1,GO,2,3]

library(tidyverse)
library(DescTools)

source("~/ADAS_processing.R") # change to name and location of your ADAS processing file
source("~/mmse_processing.R") # change to name and location of your MMSE processing file
source("~/Downloads/diagnoses_processing.R") # change to name and location of your diagnosis processing file

mmse_baseline_dates <- mmse %>%
  dplyr::ungroup() %>%
  dplyr::filter(VISCODE2 %in% c("bl","f","sc")) %>%
  dplyr::select(RID,VISDATE) %>%
  dplyr::rename(ref_date = VISDATE)

mmse <- dplyr::left_join(mmse,mmse_baseline_dates,by="RID") %>%
  dplyr::rename(mmse_date=VISDATE)

mmse <- mmse %>%
  dplyr::mutate(VISCODE_mmse = round(as.numeric(mmse_date - ref_date)/182.125))


neurobat_source <- readr::read_delim("~/Downloads/NEUROBAT_for_review_18Jun2024.csv") # change to name and location of Neuropsychological Battery

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
  dplyr::select(RID,VISCODE2,LDELTOTAL,DIGITSCOR,TRABSCOR,VISDATE) %>%
  dplyr::rename(neurobat_date = VISDATE)

neurobat_non_bl <- dplyr::left_join(neurobat_non_bl,mmse_baseline_dates,by="RID")

neurobat_non_bl <- neurobat_non_bl %>%
  dplyr::mutate(VISCODE_mmse = round(as.numeric(neurobat_date - ref_date)/182.125))

neurobat_ldel_bl <- neurobat_source %>%
  dplyr::filter(VISCODE2 %in% c("sc","f")) %>%
  dplyr::select(RID,VISCODE2,LDELTOTAL,VISDATE) %>%
  dplyr::mutate(VISCODE2="bl") %>%
  dplyr::rename(ldel_date=VISDATE) %>%
  dplyr::arrange(RID) 

neurobat_ldel_bl <- dplyr::left_join(neurobat_ldel_bl,mmse_baseline_dates,by="RID")

neurobat_ldel_bl <- neurobat_ldel_bl %>%
  dplyr::mutate(VISCODE_mmse = round(as.numeric(ldel_date - ref_date)/182.125)) %>%
  dplyr::select(-VISCODE2)

neurobat_scores_bl <- neurobat_source %>%
  dplyr::filter(VISCODE2 %in% c("bl")) %>%
  dplyr::select(RID,VISCODE2,DIGITSCOR,TRABSCOR,VISDATE) %>%
  dplyr::rename(neurobat_scores_date=VISDATE) %>%
  dplyr::arrange(RID)

neurobat_scores_bl <- dplyr::left_join(neurobat_scores_bl,mmse_baseline_dates,by="RID")

neurobat_scores_bl <- neurobat_scores_bl %>%
  dplyr::mutate(VISCODE_mmse = round(as.numeric(neurobat_scores_date - ref_date)/182.125)) %>%
  dplyr::select(-VISCODE2)

## merges LDEL and score battery records
neurobat_bl <- dplyr::full_join(neurobat_ldel_bl,neurobat_scores_bl,by=c("RID","VISCODE_mmse"))

## averages records for subjects with multiple batteries evaluated
neurobat <- dplyr::bind_rows(neurobat_bl,neurobat_non_bl) %>%
  dplyr::group_by(RID,VISCODE_mmse) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE)

adas_q4_for_pacc_1 <- adas_1_score %>%
  dplyr::select(RID,VISCODE2,Q4,EXAMDATE) %>%
  dplyr::rename(ADASQ4=Q4)

adas_q4_for_pacc_2_3_go <- adas_2_3_go %>%
  dplyr::select(RID,VISCODE2,Q4SCORE,EXAMDATE) %>%
  dplyr::rename(ADASQ4=Q4SCORE)

adas_q4_for_pacc <- dplyr::bind_rows(adas_q4_for_pacc_1,adas_q4_for_pacc_2_3_go) %>%
  dplyr::rename(adas_date = EXAMDATE)

adas_q4_for_pacc <- dplyr::left_join(adas_q4_for_pacc,mmse_baseline_dates,by="RID")

adas_q4_for_pacc <- adas_q4_for_pacc %>%
  dplyr::mutate(VISCODE_mmse = round(as.numeric(adas_date - ref_date)/182.125))

adas_q4_for_pacc <- adas_q4_for_pacc %>%
  dplyr::group_by(RID,VISCODE_mmse) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE) 

pacc_df <- dplyr::left_join(mmse %>% 
                            dplyr::select(RID,VISCODE2,VISCODE_mmse,MMSE),
                            adas_q4_for_pacc %>%
                              dplyr::select(-VISCODE2),
                            by=c("RID","VISCODE_mmse")) %>%
  dplyr::left_join(.,
                   neurobat %>%
                     dplyr::select(-VISCODE2),
                   by=c("RID","VISCODE_mmse")) %>%
  dplyr::left_join(.,
                   adni_diagnoses %>%
                     dplyr::rename(DX.bl=DX) %>%
                     dplyr::filter(VISCODE2=="bl") %>% 
                     dplyr::select(RID,DX.bl),
                   by=c("RID")) %>%
  dplyr::mutate(log.TRABSCOR=log(TRABSCOR+1)) %>%
  dplyr::select(RID,VISCODE2,VISCODE_mmse,DX.bl,ADASQ4, LDELTOTAL,DIGITSCOR,TRABSCOR,MMSE) %>%
  tidyr::drop_na(DX.bl) %>%
  dplyr::mutate(DX.bl = 
                  case_when(DX.bl == "CU" ~ "CN",
                            TRUE ~ DX.bl)) %>%
  dplyr::rename(VISCODE=VISCODE2)

pacc <- debugged_pacc_fn(dd=data.frame(pacc_df)) %>%
  dplyr::group_by(RID,VISCODE) %>%
  dplyr::mutate(across(.cols=dplyr::where(is.numeric),.fns=~mean(.,na.rm=TRUE))) %>%
  dplyr::distinct_at(vars(RID),.keep_all=TRUE) %>%
  dplyr::rename(VISCODE2=VISCODE)

date_table_pacc <- dplyr::full_join(adas_q4_for_pacc %>%
                                      dplyr::select(RID,VISCODE_mmse,adas_date),
                                    dplyr::full_join(neurobat %>%
                                                       dplyr::select(RID,VISCODE_mmse,ldel_date,neurobat_scores_date,neurobat_date),
                                                     mmse %>%
                                                       dplyr::select(RID,VISCODE2,VISCODE_mmse,mmse_date),
                                                     by=c("RID","VISCODE_mmse")
                                    ),by=c("RID","VISCODE_mmse")
)

date_table_pacc <-
  date_table_pacc %>%
  dplyr::mutate(
    date_na_count = sum(is.na(c(adas_date,neurobat_date,neurobat_scores_date,ldel_date,mmse_date))),
    pacc_date = median(c(adas_date,neurobat_date,neurobat_scores_date,ldel_date,mmse_date),
                      na.rm=TRUE))

date_table_pacc <- date_table_pacc %>%
  dplyr::mutate(pacc_date = as.Date(pacc_date,origin="1970-01-01")) %>%
  dplyr::ungroup()

pacc <- dplyr::left_join(pacc,date_table_pacc %>%
                           dplyr::select(RID,VISCODE2,pacc_date),
                         by=c("RID","VISCODE2"))