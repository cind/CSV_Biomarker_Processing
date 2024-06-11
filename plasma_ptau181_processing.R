## table name: University of Gothenburg Longitudinal Plasma P-tau181 [ADNI1,GO,2]

## Data is longitudinal
plasma_ptau <- 
  readr::read_delim("~/Projects/Amprion Project/Source Data/UGOTPTAU181_06_18_20.csv") %>% 
  dplyr::select(RID,VISCODE2,EXAMDATE,PLASMAPTAU181)

## Some subjects have two measurements at a given visit; took mean two measurements for those subjects
plasma_ptau <- plasma_ptau %>% 
  dplyr::group_by(RID,EXAMDATE) %>% 
  dplyr::mutate(PLASMAPTAU181=mean(PLASMAPTAU181)) %>% 
  dplyr::distinct_at(vars(RID,EXAMDATE),.keep_all = TRUE) 

## Add longitudinal utility variables
plasma_ptau <- plasma_ptau %>% 
  dplyr::arrange(RID,EXAMDATE) %>% 
  dplyr::group_by(RID) %>% 
  dplyr::mutate(plasma_ptau_change_from_bl=PLASMAPTAU181-PLASMAPTAU181[1L],
                baseline=PLASMAPTAU181[1L],
                time=EXAMDATE-EXAMDATE[1L]) %>% 
  dplyr::ungroup()

## recode incorrect VISCODE2
plasma_ptau$VISCODE2[which(plasma_ptau$RID=="1097"&plasma_ptau$EXAMDATE=="2013-01-31")]<-"m72"
