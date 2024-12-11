library(tidyverse)
library(mediation)

source("~/Downloads/cubbt_biomarkers_processing.R")
source("~/egfr_processing.R")
source("~/Downloads/adni_generate_VISCODE3_fn.R")
source("upenn_csf_processing.R") ## updated files
source("~/Downloads/CDR_processing.R") ## updated files
source("~/Projects/demographics_processing.R") ## updated files
source("~/xiong_medical_covariate_processing.R") ## updated files
source("~/Downloads/cubbt_biomarkers_processing.R") ## updated files
source("~/Downloads/apoe_processing.R") ## updated files
source("~/Downloads/diagnoses_processing.R")

## minimal theme for graphs
theme_blank <- function()
{
  theme_minimal() %+replace%
    theme(
      panel.background = element_rect(fill='transparent'),
      plot.background = element_rect(fill='transparent', color=NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill='transparent'),
      legend.box.background = element_rect(fill='transparent')
    )
}

## read in registry to calculate alternate VISCODE (VISCODE3)
registry <- readr::read_delim("~/Downloads/REGISTRY_04Nov2024.csv") # change path to your path

## read in AB PET
amyloid_pet <- readr::read_delim("~/Downloads/UCBERKELEY_AMY_6MM_31Jul2024.csv") %>%
  dplyr::filter(qc_flag %in% c(-1,1,2)) %>%
  dplyr::select(RID,VISCODE,SCANDATE,CENTILOIDS) %>%
  dplyr::rename(pet_date = SCANDATE)

amyloid_pet <- adni_generate_VISCODE3(df = amyloid_pet,registry_df = registry,value_cols = c("CENTILOIDS"),
                       date_col = amyloid_pet$pet_date)

amyloid_pet <- amyloid_pet %>%
  dplyr::mutate(amyloid_pos_pet = case_when(
    CENTILOIDS > 20 ~ 1,
    CENTILOIDS <= 20 ~ 0
  ))

## read in Tau PET
tau_pet <- readr::read_delim("~/Downloads/UCBERKELEY_TAU_6MM_26Sep2024.csv")

tau_pet <- tau_pet %>%
  dplyr::mutate(early_tau_pet_metaroi = 
                  ((CTX_ENTORHINAL_SUVR*CTX_ENTORHINAL_VOLUME) +
                     (CTX_PARAHIPPOCAMPAL_SUVR*CTX_PARAHIPPOCAMPAL_VOLUME) +
                     (AMYGDALA_SUVR*AMYGDALA_VOLUME))
                /(
                  CTX_ENTORHINAL_VOLUME +
                    CTX_PARAHIPPOCAMPAL_VOLUME +
                    AMYGDALA_VOLUME
                ))

tau_pet <- tau_pet %>%
  dplyr::select(RID,VISCODE,SCANDATE,early_tau_pet_metaroi,qc_flag) %>%
  dplyr::filter(qc_flag %in% c(1,2)) %>%
  dplyr::rename(tau_pet_date = SCANDATE)

tau_pet <- adni_generate_VISCODE3(df = tau_pet,registry_df = registry,value_cols = c("early_tau_pet_metaroi"),
                                      date_col = tau_pet$tau_pet_date)

tau_pet <- tau_pet %>%
  dplyr::mutate(early_tau_pet_pos = case_when(
  early_tau_pet_metaroi > 1.328 ~ 1,
  early_tau_pet_metaroi <= 1.328 ~ 0
)
)

## add VISCODE3 to other fields
upenn_csf_biomarkers <- upenn_csf_biomarkers %>%
  dplyr::rename(csf_date = EXAMDATE,
                VISCODE2_csf = VISCODE2) %>%
  dplyr::select(RID,VISCODE2_csf,csf_date,ABETA42_csf)

upenn_csf_biomarkers <- adni_generate_VISCODE3(df = upenn_csf_biomarkers,registry_df = registry,value_cols = c("ABETA42_csf"),
                                  date_col = upenn_csf_biomarkers$csf_date)

upenn_csf_biomarkers <- upenn_csf_biomarkers %>%
  dplyr::mutate(amyloid_pos_csf = case_when(
    ABETA42_csf < 980 ~ 1,
    ABETA42_csf >= 980 ~ 0
  ))

cdr <- cdr %>%
  dplyr::rename(cdr_date = VISDATE,
                VISCODE2_cdr = VISCODE2) %>%
  dplyr::select(RID,VISCODE2_cdr,cdr_date,CDRSB,CDGLOBAL) %>%
  dplyr::filter(cdr_date != "")

cdr <- adni_generate_VISCODE3(df = cdr,registry_df = registry,value_cols = c("CDRSB","CDGLOBAL"),
                                               date_col = cdr$cdr_date)

adni_diagnoses <- adni_diagnoses %>%
  dplyr::rename(dx_date = EXAMDATE,
                VISCODE2_dx = VISCODE2)

adni_diagnoses <- adni_generate_VISCODE3(df = adni_diagnoses,registry_df = registry,,
                              date_col = adni_diagnoses$dx_date)

## multiple diagnoses per VISCODE3; more conservatively take earlier DX
adni_diagnoses <- adni_diagnoses %>%
  dplyr::arrange(dx_date) %>%
  dplyr::distinct_at(vars(RID,VISCODE3),.keep_all=TRUE)

## add cutoffs from "Head-to-head comparison of leading blood tests for Alzheimer’s disease pathology"
cubbt_biomarkers_c2n <- 
  cubbt_biomarkers %>%
  dplyr::filter(manufacturer_plasma_cubbt == "C2N") %>%
  dplyr::rename(EXAMDATE_c2n = EXAMDATE)

cubbt_biomarkers_c2n <- adni_generate_VISCODE3(df = cubbt_biomarkers_c2n,registry_df = registry,value_cols = c("Abeta42_Abeta40_plasma_cubbt",
                                                                                                               "ptau217_plasma_cubbt",
                                                                                                               "ptau217_ratio_plasma_cubbt"),
                                               date_col = cubbt_biomarkers_c2n$EXAMDATE_c2n)

cubbt_biomarkers_c2n <- cubbt_biomarkers_c2n %>%
  dplyr::select(RID,VISCODE2,VISCODE3,EXAMDATE_c2n,Abeta42_Abeta40_plasma_cubbt,
                ptau217_plasma_cubbt,ptau217_ratio_plasma_cubbt) %>%
  dplyr::mutate(c2n_ab_ratio_ab_pet_pos = Abeta42_Abeta40_plasma_cubbt<0.0924,
                c2n_ptau_ratio_ab_pet_pos = ptau217_ratio_plasma_cubbt>4.06,
                c2n_ptau_217_ab_pet_pos = ptau217_plasma_cubbt>2.34,
                c2n_ab_ratio_tau_pet_pos = Abeta42_Abeta40_plasma_cubbt<0.0951,
                c2n_ptau_ratio_tau_pet_pos = ptau217_ratio_plasma_cubbt>5.60,
                c2n_ptau_217_tau_pet_pos = ptau217_plasma_cubbt>2.94) %>%
  dplyr::rename_with(~paste0(.,"_c2n"),(c(Abeta42_Abeta40_plasma_cubbt,
                                          ptau217_plasma_cubbt,ptau217_ratio_plasma_cubbt)))

cubbt_biomarkers_fuji <- 
  cubbt_biomarkers %>%
  dplyr::filter(manufacturer_plasma_cubbt == "Fuji") %>%
  dplyr::rename(EXAMDATE_fuji = EXAMDATE)

cubbt_biomarkers_fuji <- adni_generate_VISCODE3(df = cubbt_biomarkers_fuji,registry_df = registry,value_cols = c("Abeta42_Abeta40_plasma_cubbt",
                                                                                                               "ptau217_plasma_cubbt"),
                                               date_col = cubbt_biomarkers_fuji$EXAMDATE_fuji)

cubbt_biomarkers_fuji <- 
  cubbt_biomarkers_fuji %>%
  dplyr::select(RID,VISCODE2,VISCODE3,EXAMDATE_fuji,Abeta42_Abeta40_plasma_cubbt,
                ptau217_plasma_cubbt) %>%
  dplyr::mutate(fuji_ab_ratio_ab_pet_pos = Abeta42_Abeta40_plasma_cubbt<0.0869,
                fuji_ptau_217_ab_pet_pos = ptau217_plasma_cubbt>0.158,
                fuji_ab_ratio_tau_pet_pos = Abeta42_Abeta40_plasma_cubbt<0.0869, ## double checked - early Tau PET and AB PET values are actually the same in preprint
                fuji_ptau_217_tau_pet_pos = ptau217_plasma_cubbt>0.177) %>%
  dplyr::rename_with(~paste0(.,"_fuji"),(c(EXAMDATE_fuji,Abeta42_Abeta40_plasma_cubbt,
                                           ptau217_plasma_cubbt)))

cubbt_biomarkers_alzpath <- 
  cubbt_biomarkers %>%
  dplyr::filter(manufacturer_plasma_cubbt == "AlzPath") %>%
  dplyr::rename(EXAMDATE_AlzPath = EXAMDATE)

cubbt_biomarkers_alzpath <- adni_generate_VISCODE3(df = cubbt_biomarkers_alzpath,registry_df = registry,value_cols = c("ptau217_plasma_cubbt"),
                                                date_col = cubbt_biomarkers_alzpath$EXAMDATE_AlzPath)

cubbt_biomarkers_alzpath <- 
  cubbt_biomarkers_alzpath %>%
  dplyr::select(RID,VISCODE2,VISCODE3,EXAMDATE_AlzPath,
                ptau217_plasma_cubbt) %>%
  dplyr::mutate(AlzPath_ptau_217_ab_pet_pos = ptau217_plasma_cubbt>0.444,
                AlzPath_ptau_217_tau_pet_pos = ptau217_plasma_cubbt>0.559) %>%
  dplyr::rename_with(~paste0(.,"_AlzPath"),(c(ptau217_plasma_cubbt)))

cubbt_biomarkers_janssen <- 
  cubbt_biomarkers %>%
  dplyr::filter(manufacturer_plasma_cubbt == "Janssen") %>%
  dplyr::rename(EXAMDATE_Janssen = EXAMDATE)

cubbt_biomarkers_janssen <- adni_generate_VISCODE3(df = cubbt_biomarkers_janssen,registry_df = registry,value_cols = c("ptau217_plasma_cubbt"),
                                                   date_col = cubbt_biomarkers_janssen$EXAMDATE_Janssen)

cubbt_biomarkers_janssen <- 
  cubbt_biomarkers_janssen %>%
  dplyr::select(RID,VISCODE2,VISCODE3,EXAMDATE_Janssen,
                ptau217_plasma_cubbt) %>%
  dplyr::mutate(Janssen_ptau_217_ab_pet_pos = ptau217_plasma_cubbt>0.0615,
                Janssen_ptau_217_tau_pet_pos = ptau217_plasma_cubbt>0.0715) %>%
  dplyr::rename_with(~paste0(.,"_Janssen"),(c(ptau217_plasma_cubbt)))

cubbt_biomarkers_roche <- 
  cubbt_biomarkers %>%
  dplyr::filter(manufacturer_plasma_cubbt == "Roche") %>%
  dplyr::rename(EXAMDATE_roche = EXAMDATE)

cubbt_biomarkers_roche <- adni_generate_VISCODE3(df = cubbt_biomarkers_roche,registry_df = registry,value_cols = c("Abeta42_Abeta40_plasma_cubbt"),
                                                   date_col = cubbt_biomarkers_roche$EXAMDATE_roche)

cubbt_biomarkers_roche <- 
  cubbt_biomarkers_roche %>%
  dplyr::select(RID,VISCODE2,VISCODE3,EXAMDATE_roche,
                Abeta42_Abeta40_plasma_cubbt) %>%
  dplyr::mutate(roche_ab_ratio_ab_pet_pos = Abeta42_Abeta40_plasma_cubbt<0.126,
                ) %>%
  dplyr::rename_with(~paste0(.,"_roche"),(c(Abeta42_Abeta40_plasma_cubbt)))

cubbt_biomarkers_quanterix <- 
  cubbt_biomarkers %>%
  dplyr::filter(manufacturer_plasma_cubbt == "QX") %>%
  dplyr::rename(EXAMDATE_quanterix = EXAMDATE)

cubbt_biomarkers_quanterix <- adni_generate_VISCODE3(df = cubbt_biomarkers_quanterix,registry_df = registry,value_cols = c("Abeta42_Abeta40_plasma_cubbt"),
                                                 date_col = cubbt_biomarkers_quanterix$EXAMDATE_quanterix)

cubbt_biomarkers_quanterix <- 
  cubbt_biomarkers_quanterix %>%
  dplyr::select(RID,VISCODE2,VISCODE3,EXAMDATE_quanterix,
                Abeta42_Abeta40_plasma_cubbt) %>%
  dplyr::mutate(quanterix_ab_ratio_ab_pet_pos = Abeta42_Abeta40_plasma_cubbt<0.0582) %>%
  dplyr::rename_with(~paste0(.,"_quanterix"),(c(Abeta42_Abeta40_plasma_cubbt)))

## join all platforms
cubbt_pos_comparison_df <- dplyr::full_join(cubbt_biomarkers_c2n,cubbt_biomarkers_fuji,by=c("RID","VISCODE3"))
cubbt_pos_comparison_df <- dplyr::full_join(cubbt_pos_comparison_df,cubbt_biomarkers_alzpath,by=c("RID","VISCODE3"))
cubbt_pos_comparison_df <- dplyr::full_join(cubbt_pos_comparison_df,cubbt_biomarkers_roche,by=c("RID","VISCODE3"))
cubbt_pos_comparison_df <- dplyr::full_join(cubbt_pos_comparison_df,cubbt_biomarkers_quanterix,by=c("RID","VISCODE3"))
cubbt_pos_comparison_df <- dplyr::full_join(cubbt_pos_comparison_df,cubbt_biomarkers_janssen,by=c("RID","VISCODE3"))

## take the earliest record with observations for all plasma biomarkers
cubbt_pos_comparison_cs <- cubbt_pos_comparison_df %>%
  tidyr::drop_na(ends_with("pos")) %>%
  dplyr::arrange(EXAMDATE_c2n) %>%
  dplyr::distinct_at(vars(RID),.keep_all = T) 

## add gold standard for AB and then select earliest complete observation
cubbt_pos_comparison_cs_plus_ab_gold <- dplyr::left_join(cubbt_pos_comparison_df,
                                                         amyloid_pet,
                                                             by=c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_ab_gold <- cubbt_pos_comparison_cs_plus_ab_gold %>%
  tidyr::drop_na(contains("pos")) %>%
  dplyr::arrange(as.numeric(VISCODE3)) %>%
  dplyr::distinct_at(vars(RID),.keep_all = T) 

## add gold standard for Tau and then select earliest complete observation
cubbt_pos_comparison_cs_plus_tau_gold <- dplyr::left_join(cubbt_pos_comparison_df,
                                                         tau_pet,
                                                         by=c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_tau_gold <- cubbt_pos_comparison_cs_plus_tau_gold %>%
  tidyr::drop_na(contains("pos")) %>%
  dplyr::arrange(as.numeric(VISCODE3)) %>%
  dplyr::distinct_at(vars(RID),.keep_all = T) 

## create agreement matrices with and without gold standard measures + true agreement matrices

## function to create agreement matrices
create_agreement_matrix <- function(matrix_df,gold_standard=NULL,percentage=TRUE){
  pos_comparison_mat <- as.matrix(matrix_df)
  pos_comparison_mat_transpose <- t(pos_comparison_mat)
  results_matrix <- matrix(nrow=ncol(pos_comparison_mat),ncol=ncol(pos_comparison_mat))
  if(missing(gold_standard)){
  for(i in 1:ncol(pos_comparison_mat)){
    for(j in 1:nrow(pos_comparison_mat_transpose)){
      results_matrix[i,j] <- sum(pos_comparison_mat[,i] == pos_comparison_mat_transpose[j,])
    }
  }
  }
  if(!missing(gold_standard)){
    for(i in 1:ncol(pos_comparison_mat)){
      for(j in 1:nrow(pos_comparison_mat_transpose)){
        results_matrix[i,j] <- sum(pos_comparison_mat[,i]== pos_comparison_mat_transpose[j,] &
                                     pos_comparison_mat[,i] == pos_comparison_mat[,gold_standard])
      }
    }
  }
  
  results_matrix_pct <- results_matrix/nrow(pos_comparison_mat)
  
  rownames(results_matrix) <- colnames(pos_comparison_mat)
  colnames(results_matrix) <- colnames(pos_comparison_mat)
  
  rownames(results_matrix_pct) <- colnames(pos_comparison_mat)
  colnames(results_matrix_pct) <- colnames(pos_comparison_mat)
  
  if(percentage){return (results_matrix_pct)} else {return (results_matrix)}
}

## create datasets for function
cubbt_ab_pet_pos_comparison_mat <- cubbt_pos_comparison_cs %>%
                                               dplyr::select(ends_with("ab_pet_pos"))

cubbt_ab_pet_pos_plus_ab_gold_comparison_mat <- cubbt_pos_comparison_cs_plus_ab_gold %>%
  dplyr::select(ends_with("ab_pet_pos"),"amyloid_pos_pet")

cubbt_tau_pet_pos_comparison_mat <- cubbt_pos_comparison_cs %>%
  dplyr::select(ends_with("tau_pet_pos"))

cubbt_tau_pet_pos_plus_tau_gold_comparison_mat <- cubbt_pos_comparison_cs_plus_tau_gold %>%
  dplyr::select(ends_with("tau_pet_pos"),"early_tau_pet_pos")

## create agreement matrices
cubbt_ab_pos_agreement_matrix <- create_agreement_matrix(matrix_df = cubbt_ab_pet_pos_comparison_mat,
                                                                            percentage=TRUE)

cubbt_ab_pos_agreement_matrix_plus_gold_standard <- create_agreement_matrix(matrix_df = cubbt_ab_pet_pos_plus_ab_gold_comparison_mat,
                                                              percentage=TRUE)

cubbt_ab_pos_true_agreement_matrix <- create_agreement_matrix(matrix_df = cubbt_ab_pet_pos_plus_ab_gold_comparison_mat,
                                gold_standard = "amyloid_pos_pet",
                                percentage=TRUE)

cubbt_tau_pos_agreement_matrix <- create_agreement_matrix(matrix_df = cubbt_tau_pet_pos_comparison_mat,
                                                         percentage=TRUE)

cubbt_tau_pos_agreement_matrix_plus_gold_standard <- create_agreement_matrix(matrix_df = cubbt_tau_pet_pos_plus_tau_gold_comparison_mat,
                                                                            percentage=TRUE)

cubbt_tau_pos_true_agreement_matrix <- create_agreement_matrix(matrix_df = cubbt_tau_pet_pos_plus_tau_gold_comparison_mat,
                                                              gold_standard = "early_tau_pet_pos",
                                                              percentage=TRUE)

## create agreement charts

## AB PET ~ AB42/40 agreement charts
cubbt_ab_pos_agreement_matrix_ab_ratio <- cubbt_ab_pos_agreement_matrix[which(rownames(cubbt_ab_pos_agreement_matrix) %in% c("c2n_ab_ratio_ab_pet_pos","fuji_ab_ratio_ab_pet_pos","roche_ab_ratio_ab_pet_pos","quanterix_ab_ratio_ab_pet_pos")),
                                                  which(rownames(cubbt_ab_pos_agreement_matrix) %in% c("c2n_ab_ratio_ab_pet_pos","fuji_ab_ratio_ab_pet_pos","roche_ab_ratio_ab_pet_pos","quanterix_ab_ratio_ab_pet_pos"))]

colnames(cubbt_ab_pos_agreement_matrix_ab_ratio) <- c("C2N","Fuji","Roche","Quanterix")
rownames(cubbt_ab_pos_agreement_matrix_ab_ratio) <- c("C2N","Fuji","Roche","Quanterix")
cubbt_ab_pos_agreement_matrix_ab_ratio_plot <- pheatmap::pheatmap(cubbt_ab_pos_agreement_matrix_ab_ratio,display_numbers=T,fontsize=12,
                   cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("AB PET + ~ AB42/40 Agreement, N = ",nrow(cubbt_ab_pet_pos_comparison_mat)))

cubbt_ab_pos_agreement_matrix_plus_gold_standard_ab_ratio <- cubbt_ab_pos_agreement_matrix_plus_gold_standard[which(rownames(cubbt_ab_pos_agreement_matrix_plus_gold_standard) %in% c("c2n_ab_ratio_ab_pet_pos","fuji_ab_ratio_ab_pet_pos","roche_ab_ratio_ab_pet_pos","quanterix_ab_ratio_ab_pet_pos","amyloid_pos_pet")),
                                                                        which(rownames(cubbt_ab_pos_agreement_matrix_plus_gold_standard) %in% c("c2n_ab_ratio_ab_pet_pos","fuji_ab_ratio_ab_pet_pos","roche_ab_ratio_ab_pet_pos","quanterix_ab_ratio_ab_pet_pos","amyloid_pos_pet"))]

colnames(cubbt_ab_pos_agreement_matrix_plus_gold_standard_ab_ratio) <- c("C2N","Fuji","Roche","Quanterix","AB PET")
rownames(cubbt_ab_pos_agreement_matrix_plus_gold_standard_ab_ratio) <- c("C2N","Fuji","Roche","Quanterix","AB PET")
cubbt_ab_pos_agreement_matrix_plus_gold_standard_ab_ratio_plot <- pheatmap::pheatmap(cubbt_ab_pos_agreement_matrix_plus_gold_standard_ab_ratio,display_numbers=T,fontsize=12,
                                                                  cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("AB PET + ~ AB42/40 Agreement, N = ",nrow(cubbt_ab_pet_pos_plus_ab_gold_comparison_mat)))

cubbt_ab_pos_true_agreement_ab_ratio <- cubbt_ab_pos_true_agreement_matrix[which(rownames(cubbt_ab_pos_true_agreement_matrix) %in% c("c2n_ab_ratio_ab_pet_pos","fuji_ab_ratio_ab_pet_pos","roche_ab_ratio_ab_pet_pos","quanterix_ab_ratio_ab_pet_pos","amyloid_pos_pet")),
                                                                                                              which(rownames(cubbt_ab_pos_true_agreement_matrix) %in% c("c2n_ab_ratio_ab_pet_pos","fuji_ab_ratio_ab_pet_pos","roche_ab_ratio_ab_pet_pos","quanterix_ab_ratio_ab_pet_pos","amyloid_pos_pet"))]

colnames(cubbt_ab_pos_true_agreement_ab_ratio) <- c("C2N","Fuji","Roche","Quanterix","AB PET")
rownames(cubbt_ab_pos_true_agreement_ab_ratio) <- c("C2N","Fuji","Roche","Quanterix","AB PET")
cubbt_ab_pos_true_agreement_ab_ratio_plot <- pheatmap::pheatmap(cubbt_ab_pos_true_agreement_ab_ratio,display_numbers=T,fontsize=12,
                                                                                     cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("AB PET + ~ AB42/40 Correct Agreement, N = ",nrow(cubbt_ab_pet_pos_plus_ab_gold_comparison_mat)))

## AB PET ~ p-tau217 agreement charts
cubbt_ab_pos_agreement_matrix_ptau217 <- cubbt_ab_pos_agreement_matrix[which(rownames(cubbt_ab_pos_agreement_matrix) %in% c("c2n_ptau_217_ab_pet_pos","fuji_ptau_217_ab_pet_pos","AlzPath_ptau_217_ab_pet_pos","Janssen_ptau_217_ab_pet_pos")),
                                                                        which(rownames(cubbt_ab_pos_agreement_matrix) %in% c("c2n_ptau_217_ab_pet_pos","fuji_ptau_217_ab_pet_pos","AlzPath_ptau_217_ab_pet_pos","Janssen_ptau_217_ab_pet_pos"))]

colnames(cubbt_ab_pos_agreement_matrix_ptau217) <- c("C2N","Fuji","AlzPath","Janssen")
rownames(cubbt_ab_pos_agreement_matrix_ptau217) <- c("C2N","Fuji","AlzPath","Janssen")
cubbt_ab_pos_agreement_matrix_ptau217_plot <- pheatmap::pheatmap(cubbt_ab_pos_agreement_matrix_ptau217,display_numbers=T,fontsize=12,
                                                                  cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("AB PET + ~ p-tau217 Agreement, N = ",nrow(cubbt_ab_pet_pos_comparison_mat)))

cubbt_ab_pos_agreement_matrix_plus_gold_standard_ptau217 <- cubbt_ab_pos_agreement_matrix_plus_gold_standard[which(rownames(cubbt_ab_pos_agreement_matrix_plus_gold_standard) %in% c("c2n_ptau_217_ab_pet_pos","fuji_ptau_217_ab_pet_pos","AlzPath_ptau_217_ab_pet_pos","Janssen_ptau_217_ab_pet_pos","amyloid_pos_pet")),
                                                                                                              which(rownames(cubbt_ab_pos_agreement_matrix_plus_gold_standard) %in% c("c2n_ptau_217_ab_pet_pos","fuji_ptau_217_ab_pet_pos","AlzPath_ptau_217_ab_pet_pos","Janssen_ptau_217_ab_pet_pos","amyloid_pos_pet"))]

colnames(cubbt_ab_pos_agreement_matrix_plus_gold_standard_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","AB PET")
rownames(cubbt_ab_pos_agreement_matrix_plus_gold_standard_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","AB PET")
cubbt_ab_pos_agreement_matrix_plus_gold_standard_ptau217_plot <- pheatmap::pheatmap(cubbt_ab_pos_agreement_matrix_plus_gold_standard_ptau217,display_numbers=T,fontsize=15,
                                                                                     cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("AB PET + ~ p-tau217 Agreement, N = ",nrow(cubbt_ab_pet_pos_plus_ab_gold_comparison_mat)))

cubbt_ab_pos_true_agreement_ptau217 <- cubbt_ab_pos_true_agreement_matrix[which(rownames(cubbt_ab_pos_true_agreement_matrix) %in% c("c2n_ptau_217_ab_pet_pos","fuji_ptau_217_ab_pet_pos","AlzPath_ptau_217_ab_pet_pos","Janssen_ptau_217_ab_pet_pos","amyloid_pos_pet")),
                                                                           which(rownames(cubbt_ab_pos_true_agreement_matrix) %in% c("c2n_ptau_217_ab_pet_pos","fuji_ptau_217_ab_pet_pos","AlzPath_ptau_217_ab_pet_pos","Janssen_ptau_217_ab_pet_pos","amyloid_pos_pet"))]

colnames(cubbt_ab_pos_true_agreement_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","AB PET")
rownames(cubbt_ab_pos_true_agreement_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","AB PET")
cubbt_ab_pos_true_agreement_ptau217_plot <- pheatmap::pheatmap(cubbt_ab_pos_true_agreement_ptau217,display_numbers=T,fontsize=12,
                                                                cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("AB PET + ~ p-tau217 Correct Agreement, N = ",nrow(cubbt_ab_pet_pos_plus_ab_gold_comparison_mat)))

## Early Tau PET ~ p-tau217 agreement charts
cubbt_tau_pos_agreement_matrix_ptau217 <- cubbt_tau_pos_agreement_matrix[which(rownames(cubbt_tau_pos_agreement_matrix) %in% c("c2n_ptau_217_tau_pet_pos","fuji_ptau_217_tau_pet_pos","AlzPath_ptau_217_tau_pet_pos","Janssen_ptau_217_tau_pet_pos")),
                                                                       which(rownames(cubbt_tau_pos_agreement_matrix) %in% c("c2n_ptau_217_tau_pet_pos","fuji_ptau_217_tau_pet_pos","AlzPath_ptau_217_tau_pet_pos","Janssen_ptau_217_tau_pet_pos"))]

colnames(cubbt_tau_pos_agreement_matrix_ptau217) <- c("C2N","Fuji","AlzPath","Janssen")
rownames(cubbt_tau_pos_agreement_matrix_ptau217) <- c("C2N","Fuji","AlzPath","Janssen")
cubbt_tau_pos_agreement_matrix_ptau217_plot <- pheatmap::pheatmap(cubbt_tau_pos_agreement_matrix_ptau217,display_numbers=T,fontsize=12,
                                                                 cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("Early Tau PET + ~ p-tau217 Agreement, N = ",nrow(cubbt_tau_pet_pos_comparison_mat)))

cubbt_tau_pos_agreement_matrix_plus_gold_standard_ptau217 <- cubbt_tau_pos_agreement_matrix_plus_gold_standard[which(rownames(cubbt_tau_pos_agreement_matrix_plus_gold_standard) %in% c("c2n_ptau_217_tau_pet_pos","fuji_ptau_217_tau_pet_pos","AlzPath_ptau_217_tau_pet_pos","Janssen_ptau_217_tau_pet_pos","early_tau_pet_pos")),
                                                                                                             which(rownames(cubbt_tau_pos_agreement_matrix_plus_gold_standard) %in% c("c2n_ptau_217_tau_pet_pos","fuji_ptau_217_tau_pet_pos","AlzPath_ptau_217_tau_pet_pos","Janssen_ptau_217_tau_pet_pos","early_tau_pet_pos"))]

colnames(cubbt_tau_pos_agreement_matrix_plus_gold_standard_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","Tau PET")
rownames(cubbt_tau_pos_agreement_matrix_plus_gold_standard_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","Tau PET")
cubbt_tau_pos_agreement_matrix_plus_gold_standard_ptau217_plot <- pheatmap::pheatmap(cubbt_tau_pos_agreement_matrix_plus_gold_standard_ptau217,display_numbers=T,fontsize=12,
                                                                                    cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("Early Tau PET + ~ p-tau217 Agreement, N = ",nrow(cubbt_tau_pet_pos_plus_tau_gold_comparison_mat)))

cubbt_tau_pos_true_agreement_ptau217 <- cubbt_tau_pos_true_agreement_matrix[which(rownames(cubbt_tau_pos_true_agreement_matrix) %in% c("c2n_ptau_217_tau_pet_pos","fuji_ptau_217_tau_pet_pos","AlzPath_ptau_217_tau_pet_pos","Janssen_ptau_217_tau_pet_pos","early_tau_pet_pos")),
                                                                          which(rownames(cubbt_tau_pos_true_agreement_matrix) %in% c("c2n_ptau_217_tau_pet_pos","fuji_ptau_217_tau_pet_pos","AlzPath_ptau_217_tau_pet_pos","Janssen_ptau_217_tau_pet_pos","early_tau_pet_pos"))]

colnames(cubbt_tau_pos_true_agreement_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","Tau PET")
rownames(cubbt_tau_pos_true_agreement_ptau217) <- c("C2N","Fuji","AlzPath","Janssen","Tau PET")
cubbt_tau_pos_true_agreement_ptau217_plot <- pheatmap::pheatmap(cubbt_tau_pos_true_agreement_ptau217,display_numbers=T,fontsize=12,
                                                               cluster_cols = FALSE,cluster_rows = FALSE,main = paste0("Early Tau PET + ~ p-tau217 Correct Agreement, N = ",nrow(cubbt_tau_pet_pos_plus_tau_gold_comparison_mat)))

## comparing variables within platforms
## AB PET
cubbt_pos_comparison_cs_plus_ab_gold_comparison <- cubbt_pos_comparison_cs_plus_ab_gold %>%
  dplyr::mutate(quanterix_ab_ratio_ab_pet_correct = quanterix_ab_ratio_ab_pet_pos == amyloid_pos_pet, 
                quanterix_ab_ratio_ab_pet_correct_text = case_when(
                  quanterix_ab_ratio_ab_pet_correct ~ "Correct",
                  !quanterix_ab_ratio_ab_pet_correct ~ "Incorrect"),
                roche_ab_ratio_ab_pet_correct = roche_ab_ratio_ab_pet_pos == amyloid_pos_pet,
                roche_ab_ratio_ab_pet_correct_text = case_when(
                  roche_ab_ratio_ab_pet_correct ~ "Correct",
                  !roche_ab_ratio_ab_pet_correct ~ "Incorrect"),
                fuji_ab_ratio_ab_pet_correct = fuji_ab_ratio_ab_pet_pos == amyloid_pos_pet,
                fuji_ab_ratio_ab_pet_correct_text = case_when(
                  fuji_ab_ratio_ab_pet_correct ~ "Correct",
                  !fuji_ab_ratio_ab_pet_correct ~ "Incorrect"),
                c2n_ab_ratio_ab_pet_correct = c2n_ab_ratio_ab_pet_pos == amyloid_pos_pet,
                c2n_ab_ratio_ab_pet_correct_text = case_when(
                  c2n_ab_ratio_ab_pet_correct ~ "Correct",
                  !c2n_ab_ratio_ab_pet_correct ~ "Incorrect"),
                ab_ratio_ab_pet_agreement_count = 
                  quanterix_ab_ratio_ab_pet_correct + roche_ab_ratio_ab_pet_correct + fuji_ab_ratio_ab_pet_correct + c2n_ab_ratio_ab_pet_correct,
                ab_ratio_ab_pet_agreement_cat = case_when(
                  ab_ratio_ab_pet_agreement_count == 4 ~ "All Correct",
                  ab_ratio_ab_pet_agreement_count == 0 ~ "None Correct",
                  ab_ratio_ab_pet_agreement_count > 0 ~ "Some Correct"),
                c2n_ptau217_ab_pet_correct = (c2n_ptau_217_ab_pet_pos == amyloid_pos_pet),
                c2n_ptau217_ab_pet_correct_text = case_when(
                  c2n_ptau217_ab_pet_correct ~ "Correct",
                  !c2n_ptau217_ab_pet_correct ~ "Incorrect"),
                fuji_ptau217_ab_pet_correct = (fuji_ptau_217_ab_pet_pos == amyloid_pos_pet),
                fuji_ptau217_ab_pet_correct_text = case_when(
                  fuji_ptau217_ab_pet_correct ~ "Correct",
                  !fuji_ptau217_ab_pet_correct ~ "Incorrect"),
                alzpath_ptau217_ab_pet_correct = (AlzPath_ptau_217_ab_pet_pos == amyloid_pos_pet),
                alzpath_ptau217_ab_pet_correct_text = case_when(
                  alzpath_ptau217_ab_pet_correct ~ "Correct",
                  !alzpath_ptau217_ab_pet_correct ~ "Incorrect"),
                janssen_ptau217_ab_pet_correct = (Janssen_ptau_217_ab_pet_pos == amyloid_pos_pet),
                janssen_ptau217_ab_pet_correct_text = case_when(
                  janssen_ptau217_ab_pet_correct ~ "Correct",
                  !janssen_ptau217_ab_pet_correct ~ "Incorrect"),
                ptau_217_ab_pet_agreement_count = 
                  c2n_ptau217_ab_pet_correct +
                  fuji_ptau217_ab_pet_correct +
                  alzpath_ptau217_ab_pet_correct +
                  janssen_ptau217_ab_pet_correct,
                ptau_217_ab_pet_agreement_cat = case_when(
                  ptau_217_ab_pet_agreement_count == 4 ~ "All Correct",
                  ptau_217_ab_pet_agreement_count == 0 ~ "None Correct",
                  ptau_217_ab_pet_agreement_count > 0 ~ "Some Correct"),
                c2n_ptau_ratio_ab_pet_correct = (c2n_ptau_ratio_ab_pet_pos == amyloid_pos_pet),
                c2n_ptau_ratio_ab_pet_correct_text = case_when(
                  c2n_ptau_ratio_ab_pet_correct ~ "Correct",
                  !c2n_ptau_ratio_ab_pet_correct ~ "Incorrect"))

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,adni_diagnoses,by = c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,cdr,by = c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,demographics %>%
                                                               dplyr::select(-VISCODE2),by="RID")
cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,apoeres,by="RID")

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,lab_data,by="RID")

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,egfr_df,by="RID")

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_ab_gold_comparison,hachinski %>%
                                                               dplyr::select(RID,hypertension_ind,stroke_ind),by="RID")

cubbt_pos_comparison_cs_plus_ab_gold_comparison$age <- round(as.numeric(lubridate::year(cubbt_pos_comparison_cs_plus_ab_gold_comparison$EXAMDATE_c2n)) -
                                                        as.numeric(cubbt_pos_comparison_cs_plus_ab_gold_comparison$birth_year) +
                                                        ((as.numeric(lubridate::month(cubbt_pos_comparison_cs_plus_ab_gold_comparison$EXAMDATE_c2n)) - 
                                                            as.numeric(cubbt_pos_comparison_cs_plus_ab_gold_comparison$birth_month)) / 12), digits=1)

cubbt_pos_comparison_cs_plus_ab_gold_comparison$age_band <- cut(cubbt_pos_comparison_cs_plus_ab_gold_comparison$age,breaks = seq(55,95,5))

cubbt_pos_comparison_cs_plus_ab_gold_comparison <- cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
  dplyr::mutate(
    white_ind = case_when(
      race == "White" ~ "White",
      race != "White" ~ "Non-White"
    ),
    cognitively_impaired = case_when(
      DX == "CU" ~ "Cognitively Unimpaired",
      DX != "CU" ~ "Cognitively Impaired"
  )
  )

count_ab_ratio_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                         dplyr::select(ab_ratio_ab_pet_agreement_count,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                       by=ab_ratio_ab_pet_agreement_count) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("No. of AB42/40 Platforms Correctly Predicting AB PET")

c2n_ab_ratio_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                         dplyr::select(c2n_ab_ratio_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                       by=c2n_ab_ratio_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ C2N AB42/40 Correct?")

fuji_ab_ratio_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                         dplyr::select(fuji_ab_ratio_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                       by=fuji_ab_ratio_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ Fuji AB42/40 Correct?")

quanterix_ab_ratio_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                             dplyr::select(quanterix_ab_ratio_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                           by=quanterix_ab_ratio_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ Quanterix AB42/40 Correct?")

roche_ab_ratio_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                             dplyr::select(roche_ab_ratio_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                           by=roche_ab_ratio_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ Roche AB42/40 Correct?")

count_ptau217_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                              dplyr::select(ptau_217_ab_pet_agreement_count,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                            by=ptau_217_ab_pet_agreement_count) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("No. of p-tau217 Platforms Correctly Predicting AB PET")

c2n_ptau217_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                            dplyr::select(c2n_ptau217_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                          by=c2n_ptau217_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ C2N p-tau217 Correct?")

fivenum(cubbt_pos_comparison_cs_plus_ab_gold_comparison$CENTILOIDS[which(cubbt_pos_comparison_cs_plus_ab_gold_comparison$c2n_ptau217_ab_pet_correct==TRUE)])

alzpath_ptau217_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                             dplyr::select(alzpath_ptau217_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                           by=alzpath_ptau217_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ AlzPath p-tau217 Correct?")

fuji_ptau217_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                               dplyr::select(fuji_ptau217_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                             by=fuji_ptau217_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ Fuji p-tau217 Correct?")

janssen_ptau217_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                               dplyr::select(janssen_ptau217_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                             by=janssen_ptau217_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ Janssen p-tau217 Correct?")

c2n_ptau_ratio_ab_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
                                                               dplyr::select(c2n_ptau_ratio_ab_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                             by=c2n_ptau_ratio_ab_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("AB PET ~ C2N p-tau217% Correct?")

## 
cubbt_pos_comparison_cs_plus_tau_gold_comparison <- cubbt_pos_comparison_cs_plus_tau_gold %>%
  dplyr::mutate(c2n_ptau217_tau_pet_correct = (c2n_ptau_217_tau_pet_pos == early_tau_pet_pos),
                c2n_ptau217_tau_pet_correct_text = case_when(
                  c2n_ptau217_tau_pet_correct ~ "Correct",
                  !c2n_ptau217_tau_pet_correct ~ "Incorrect"),
                fuji_ptau217_tau_pet_correct = (fuji_ptau_217_tau_pet_pos == early_tau_pet_pos),
                fuji_ptau217_tau_pet_correct_text = case_when(
                  fuji_ptau217_tau_pet_correct ~ "Correct",
                  !fuji_ptau217_tau_pet_correct ~ "Incorrect"),
                alzpath_ptau217_tau_pet_correct = (AlzPath_ptau_217_tau_pet_pos == early_tau_pet_pos),
                alzpath_ptau217_tau_pet_correct_text = case_when(
                  alzpath_ptau217_tau_pet_correct ~ "Correct",
                  !alzpath_ptau217_tau_pet_correct ~ "Incorrect"),
                janssen_ptau217_tau_pet_correct = (Janssen_ptau_217_tau_pet_pos == early_tau_pet_pos),
                janssen_ptau217_tau_pet_correct_text = case_when(
                  janssen_ptau217_tau_pet_correct ~ "Correct",
                  !janssen_ptau217_tau_pet_correct ~ "Incorrect"),
                ptau_217_tau_pet_agreement_count = 
                  c2n_ptau217_tau_pet_correct +
                  fuji_ptau217_tau_pet_correct +
                  alzpath_ptau217_tau_pet_correct +
                  janssen_ptau217_tau_pet_correct,
                ptau_217_tau_pet_agreement_cat = case_when(
                  ptau_217_tau_pet_agreement_count == 4 ~ "All Correct",
                  ptau_217_tau_pet_agreement_count == 0 ~ "None Correct",
                  ptau_217_tau_pet_agreement_count > 0 ~ "Some Correct"),
                c2n_ptau_ratio_tau_pet_correct = (c2n_ptau_ratio_ab_pet_pos == early_tau_pet_pos),
                c2n_ptau_ratio_tau_pet_correct_text = case_when(
                  c2n_ptau_ratio_tau_pet_correct ~ "Correct",
                  !c2n_ptau_ratio_tau_pet_correct ~ "Incorrect"))

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,adni_diagnoses,by = c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,cdr,by = c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,amyloid_pet,by = c("RID","VISCODE3"))

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,demographics %>%
                                                                      dplyr::select(-VISCODE2),by="RID")
cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,apoeres,by="RID")

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,lab_data,by="RID")

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,egfr_df,by="RID")

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- dplyr::left_join(cubbt_pos_comparison_cs_plus_tau_gold_comparison,hachinski %>%
                                                                      dplyr::select(RID,hypertension_ind,stroke_ind),by="RID")

cubbt_pos_comparison_cs_plus_tau_gold_comparison$age <- round(as.numeric(lubridate::year(cubbt_pos_comparison_cs_plus_tau_gold_comparison$EXAMDATE_c2n)) -
                                                               as.numeric(cubbt_pos_comparison_cs_plus_tau_gold_comparison$birth_year) +
                                                               ((as.numeric(lubridate::month(cubbt_pos_comparison_cs_plus_tau_gold_comparison$EXAMDATE_c2n)) - 
                                                                   as.numeric(cubbt_pos_comparison_cs_plus_tau_gold_comparison$birth_month)) / 12), digits=1)

cubbt_pos_comparison_cs_plus_tau_gold_comparison$age_band <- cut(cubbt_pos_comparison_cs_plus_tau_gold_comparison$age,breaks = seq(55,95,5))

cubbt_pos_comparison_cs_plus_tau_gold_comparison <- cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
  dplyr::mutate(
    white_ind = case_when(
      race == "White" ~ "White",
      race != "White" ~ "Non-White"
    ),
    cognitively_impaired = case_when(
      DX == "CU" ~ "Cognitively Unimpaired",
      DX != "CU" ~ "Cognitively Impaired"
    )
  )

count_ptau217_tau_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
                                                             dplyr::select(ptau_217_tau_pet_agreement_count,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                           by=ptau_217_tau_pet_agreement_count) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("No. of p-tau217 Platforms Correctly Predicting Early Tau PET")

c2n_ptau217_tau_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
                                                           dplyr::select(c2n_ptau217_tau_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                         by=c2n_ptau217_tau_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("Early Tau PET ~ C2N p-tau217 Correct?")

alzpath_ptau217_tau_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
                                                               dplyr::select(alzpath_ptau217_tau_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                             by=alzpath_ptau217_tau_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("Early Tau PET ~ AlzPath p-tau217 Correct?")

fuji_ptau217_tau_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
                                                            dplyr::select(fuji_ptau217_tau_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                          by=fuji_ptau217_tau_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("Early Tau PET ~ Fuji p-tau217 Correct?")

janssen_ptau217_tau_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
                                                               dplyr::select(janssen_ptau217_tau_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                             by=janssen_ptau217_tau_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("Early Tau PET ~ Janssen p-tau217 Correct?")

c2n_ptau_ratio_tau_pet_correct_tbl <- gtsummary::tbl_summary(cubbt_pos_comparison_cs_plus_tau_gold_comparison %>%
                                                              dplyr::select(c2n_ptau_ratio_tau_pet_correct_text,white_ind,age,sex,apoe4_status,cognitively_impaired,CDRSB,yrs_education,ckd_ind,diabetes_ind,hypertension_ind,stroke_ind,CENTILOIDS),
                                                            by=c2n_ptau_ratio_tau_pet_correct_text) %>%
  gtsummary::add_p() %>%
  gtsummary::add_q() %>%
  gtsummary::modify_caption("Early Tau PET ~ C2N p-tau217% Correct?")

## testing impact of age on p-tau217 variables
mediator_model <- lm(CENTILOIDS ~ age,
                     data = cubbt_pos_comparison_cs_plus_ab_gold_comparison)
exposure_model_c2n <- lm(ptau217_plasma_cubbt_c2n ~ CENTILOIDS + age,
                     data = cubbt_pos_comparison_cs_plus_ab_gold_comparison)

mediation_results_c2n <- mediation::mediate(mediator_model,exposure_model_c2n,
                                        treat = "age",mediator = "CENTILOIDS")

exposure_model_fuji <- lm(ptau217_plasma_cubbt_fuji ~ CENTILOIDS + age,
                         data = cubbt_pos_comparison_cs_plus_ab_gold_comparison)

mediation_results_fuji <- mediation::mediate(mediator_model,exposure_model_fuji,
                                            treat = "age",mediator = "CENTILOIDS")

exposure_model_alzpath <- lm(ptau217_plasma_cubbt_AlzPath ~ CENTILOIDS + age,
                          data = cubbt_pos_comparison_cs_plus_ab_gold_comparison)

mediation_results_alzpath <- mediation::mediate(mediator_model,exposure_model_alzpath,
                                             treat = "age",mediator = "CENTILOIDS")

exposure_model_janssen <- lm(ptau217_plasma_cubbt_Janssen ~ CENTILOIDS + age,
                             data = cubbt_pos_comparison_cs_plus_ab_gold_comparison)

mediation_results_janssen <- mediation::mediate(mediator_model,exposure_model_janssen,
                                                boot = TRUE, sims = 1000,
                                                treat = "age",mediator = "CENTILOIDS")

## graph comparing AB PET ~ p-tau217 and AB42/40
centiloid_ab_ratio_ptau_comparison_df <-cubbt_pos_comparison_cs_plus_ab_gold_comparison %>%
  dplyr::mutate(centiloid_band = as.factor(case_when(
    CENTILOIDS <= 0 ~ "0 or lower",
    CENTILOIDS <= 20 ~ "1 to 19",
    CENTILOIDS < 50 ~ "20 to 40",
    CENTILOIDS >=50 ~ "40 or higher"
  ))) %>%
  dplyr::group_by(centiloid_band) %>%
  dplyr::summarize(across(tidyr::contains("ab_pet_correct"),list(mean=mean,sd=sd),.names="{.col}{.fn}"),count=n()) %>%
  tidyr::pivot_longer(cols = tidyr::contains("ab_pet_correct"),names_to = c("biomarker",".value"),names_sep = "correct") %>%
  dplyr::mutate(se = sqrt((mean*(1-mean))/count),
                moe = se*1.96,
                upper = mean+moe,
                lower = mean-moe,
                biomarker_cat = case_when(
                  stringr::str_detect(biomarker,"ab_ratio") ~ "AB42/40",
                  stringr::str_detect(biomarker,"ptau217") ~ "p-tau217"),
                platform = case_when(
                  stringr::str_detect(biomarker,"quanterix") ~ "Quanterix",
                  stringr::str_detect(biomarker,"roche") ~ "Roche",
                  stringr::str_detect(biomarker,"janssen") ~ "Janssen",
                  stringr::str_detect(biomarker,"AlzPath") ~ "AlzPath",
                  stringr::str_detect(biomarker,"fuji") ~ "Fuji",
                  stringr::str_detect(biomarker,"c2n") ~ "C2N"),
                biomarker_name = paste(platform,biomarker_cat)) %>%
  tidyr::drop_na(centiloid_band,biomarker_cat,platform)

centiloid_ab_ratio_ptau_comparison_plot <- ggplot2::ggplot(aes(x=centiloid_band,y=mean,color=biomarker_cat,shape=platform),data=centiloid_ab_ratio_ptau_comparison_df) +
  theme_blank() +
  ggplot2::geom_point(size=3,position=position_dodge(width=.5)) +
  ggplot2::geom_errorbar(aes(ymin=lower,ymax=upper,width=0.2),position=position_dodge(width=.5)) +
  ggplot2::ylab("Accuracy") +
  ggplot2::xlab("Centiloid Range") +
  ggplot2::theme(text=element_text(size=15))