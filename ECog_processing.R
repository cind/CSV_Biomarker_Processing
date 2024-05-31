############################################################################################################
#E-COG: Subject
############################################################################################################
d_s <- read.csv("`/ECOGPT_08Apr2024.csv")
d_s[d_s == "9"] <- NA
total_s <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8',
             'LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9',
             'VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
             'VISSPAT6', 'VISSPAT7', 'VISSPAT8',
             'PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5',
             'ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6',
             'DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
total.coef_s <- c(0.89,0.83,0.83,0.87,0.82,0.90,0.83,0.90,0.73,0.78,0.72,0.81,0.75,
                  0.63,0.55,0.68,0.65,0.82,0.84,0.91,0.94,0.89,0.90,0.89,0.95,0.86,
                  0.91,0.80,0.87,0.74,0.91,0.91,0.90,0.91,0.85,0.85,0.88,0.84,0.85)
mem_s <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8')
mem.coef_s <- c(0.26,0.45,0.47,0.21,0.47,0.23,0.41,0.17)
lang_s <-  c('LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9')
lang.coef_s <- c(0.55,0.45,0.55,0.45,0.44,0.63,0.71,0.62,0.67)
visspat_s <- c('VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
               'VISSPAT6', 'VISSPAT7', 'VISSPAT8')
visspat.coef_s <- c(0.56,0.51,0.25,0.23,0.35,0.37,0.33)
plan_s <- c('PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5')
plan.coef_s <- c(0.20,0.27,0.29,0.58,0.42)
organ_s <- c('ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6')
organ.coef_s <- c(0.46,0.35,0.38,0.36,0.23,0.38)
divatt_s <- c('DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
divatt.coef_s <- c(0.42,0.37,0.42,0.44)

d_s$EcogGlobal <- rowMeans(d_s[,total_s]*total.coef_s,na.rm=T)

############################################################################################################
#E-COG: Study Partner
############################################################################################################
d_p <- read.csv("~/ECOGSP_08Apr2024.csv")
d_p[d_p == "9"] <- NA
total_p <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8',
             'LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9',
             'VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
             'VISSPAT6', 'VISSPAT7', 'VISSPAT8',
             'PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5',
             'ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6',
             'DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
total.coef_p <- c(0.89,0.83,0.83,0.87,0.82,0.90,0.83,0.90,0.73,0.78,0.72,0.81,0.75,
                  0.63,0.55,0.68,0.65,0.82,0.84,0.91,0.94,0.89,0.90,0.89,0.95,0.86,
                  0.91,0.80,0.87,0.74,0.91,0.91,0.90,0.91,0.85,0.85,0.88,0.84,0.85)
mem_p <- c('MEMORY1', 'MEMORY2', 'MEMORY3', 'MEMORY4', 'MEMORY5', 'MEMORY6', 'MEMORY7', 'MEMORY8')
mem.coef_p <- c(0.26,0.45,0.47,0.21,0.47,0.23,0.41,0.17)
lang_p <-  c('LANG1', 'LANG2', 'LANG3', 'LANG4', 'LANG5', 'LANG6', 'LANG7', 'LANG8', 'LANG9')
lang.coef_p <- c(0.55,0.45,0.55,0.45,0.44,0.63,0.71,0.62,0.67)
visspat_p <- c('VISSPAT1', 'VISSPAT2', 'VISSPAT3', 'VISSPAT4',
               'VISSPAT6', 'VISSPAT7', 'VISSPAT8')
visspat.coef_p <- c(0.56,0.51,0.25,0.23,0.35,0.37,0.33)
plan_p <- c('PLAN1', 'PLAN2', 'PLAN3', 'PLAN4', 'PLAN5')
plan.coef_p <- c(0.20,0.27,0.29,0.58,0.42)
organ_p <- c('ORGAN1', 'ORGAN2', 'ORGAN3', 'ORGAN4', 'ORGAN5', 'ORGAN6')
organ.coef_p <- c(0.46,0.35,0.38,0.36,0.23,0.38)
divatt_p <- c('DIVATT1', 'DIVATT2', 'DIVATT3', 'DIVATT4')
divatt.coef_p <- c(0.42,0.37,0.42,0.44)

d_p$EcogGlobal <- rowMeans(d_p[,total_p]*total.coef_p,na.rm=T)
