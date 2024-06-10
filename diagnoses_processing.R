#adding diagnoses as CN, MCI, or AD
library(dplyr)

diagnoses <- read.csv("~/DXSUM_PDXCONV.csv")
diagnoses[diagnoses == ""] <- NA

diagnoses_all <- diagnoses %>%
  dplyr::mutate(diagnosis = dplyr::case_when(DIAGNOSIS==1 | DXCHANGE==1 | DXCHANGE==7 | DXCHANGE==9 | DXCURREN == 1 ~ 'CN',
                                         DIAGNOSIS==2 | DXCHANGE==2 | DXCHANGE==4 | DXCHANGE==8 | DXCURREN == 2 ~ 'MCI',
                                         DIAGNOSIS==3 | DXCHANGE==3 | DXCHANGE==5 | DXCHANGE==6 | DXCURREN == 3 ~ 'AD')) %>%
  dplyr::rename(DX.DATE = EXAMDATE)
