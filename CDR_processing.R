library(dplyr) 

cdr <- read.csv("~/CDR_24Apr2024.csv")
cdr <- cdr[cdr$CDGLOBAL>=0,]
cdr$CDRSB <- rowSums(cdr[,c('CDMEMORY','CDORIENT','CDJUDGE','CDCOMMUN','CDHOME','CDCARE')])
cdr <- cdr %>%
  dplyr::rename(VISDATE = EXAMDATE)
