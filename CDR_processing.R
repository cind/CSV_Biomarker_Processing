## table name on LONI: 	Clinical Dementia Rating [ADNI1,GO,2,3,4]

library(dplyr) 

cdr <- read.csv("~/CDR_24Apr2024.csv")
cdr <- cdr[cdr$CDGLOBAL>=0,]
cdr$CDRSB <- rowSums(cdr[,c('CDMEMORY','CDORIENT','CDJUDGE','CDCOMMUN','CDHOME','CDCARE')])
cdr <- cdr %>%
  dplyr::rename(VISDATE = EXAMDATE)
