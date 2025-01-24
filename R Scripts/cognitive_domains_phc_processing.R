## table name: ADSP Phenotype Harmonization Consortium (PHC) - Composite Cognitive Scores, Version: 2023-12-18

## PreciseFilter is 1 when values are within an acceptable range of the standard error of the mean
## The visuospatial domain (PHC_VSP) was excluded because a very high number of values were not within the acceptable range
## PHC_VSP can be added by removing comments from relevant rows

cognitive_domains_phc <- readr::read_delim("~/ADSP_PHC_COGN_Dec2023.csv") %>% # edit to your source file name + location 
  dplyr::select(RID,VISCODE2,EXAMDATE,
                PHC_EXF,PHC_MEM,PHC_LAN,
                #PHC_VSP,PHC_VSP_PreciseFilter,
                PHC_EXF_PreciseFilter,PHC_MEM_PreciseFilter,PHC_LAN_PreciseFilter) %>%
  dplyr::mutate(PHC_EXF = case_when(
    PHC_EXF_PreciseFilter==1 ~ PHC_EXF,
    PHC_EXF_PreciseFilter==0 ~ NA
  ),
  PHC_MEM = case_when(
    PHC_MEM_PreciseFilter==1 ~ PHC_MEM,
    PHC_MEM_PreciseFilter==0 ~ NA
  ),
#  PHC_VSP = case_when(
#    PHC_VSP_PreciseFilter==1 ~ PHC_VSP,
#    PHC_VSP_PreciseFilter==0 ~ NA
#  ),
  PHC_LAN = case_when(
    PHC_LAN_PreciseFilter==1 ~ PHC_LAN,
    PHC_LAN_PreciseFilter==0 ~ NA
  ))
