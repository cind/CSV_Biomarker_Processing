#this creates a column with both APOE alleles combined
apoeres <- read.csv("~/APOERES.csv") 
apoeres$apoe <- paste(paste("E", apoeres$APGEN1, sep = ""), paste("E", apoeres$APGEN2, sep = ""), sep = "/")

## creates five additional utility variables for APOE information
## APOE2, APOE3, APOE4: the total number of that allele each subject has
## apoe4_status: whether a subject has at least one APOE-e4 allele
## apoe4_status_text: a leveled factor for each subject's APOE-e4 genotype
apoeres <- apoeres %>% 
  dplyr::mutate(APOE2 = as.numeric(APGEN1==2)+as.numeric(APGEN2==2),
                APOE3 = as.numeric(APGEN1==3)+as.numeric(APGEN2==3),
                APOE4 = as.numeric(APGEN1==4)+as.numeric(APGEN2==4)) %>%
  dplyr::select(RID,APOE4) %>% dplyr::mutate(apoe4_status=(APOE4>0),
                                             apoe4_status_text = factor(case_when(
                                               APOE4==0 ~ "Non-Carrier",
                                               APOE4==1 ~ "Heterozygotes",
                                               APOE4==2 ~ "Homozygotes")
                                               ,levels = c("Non-Carrier","Heterozygotes","Homozygotes"))) %>%
  dplyr::distinct_at(vars(RID),.keep_all = TRUE)
