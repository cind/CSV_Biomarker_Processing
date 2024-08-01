## table name on LONI: 	Clinical Dementia Rating [ADNI1,GO,2,3,4]

library(tidyverse) 

cdr <- read.csv("~/Downloads/CDR_31Jul2024.csv")
cdr <- cdr[cdr$CDGLOBAL>=0,]
cdr$CDRSB <- rowSums(cdr[,c('CDMEMORY','CDORIENT','CDJUDGE','CDCOMMUN','CDHOME','CDCARE')])
