source("/Users/zackhausle/apoe_sex_project_header.R")

# WMH ----
# # Changed this to WMH with ICV (called wmlICV here) 
# # For analysis WMH will be used, 
# # and for other measures WML=WMH/median_wmlICV_by_subject will be used as covariate

wml=getWML()
colnames(wml)[colnames(wml)=='wml.DATE']='EXAMDATE'

wml <- wml %>%
  dplyr::mutate(
    log_WMH = log(WMH)
  )

dataset_builder <- function(df,response_var){
  
apoe=getAPOE()
df <- dplyr::left_join(df,apoe,by="RID")

dem=getDemog()
df <- dplyr::left_join(df,dem,by="RID")
df$dt=abs(as.Date(df$dem.DATE,format="%Y-%m-%d")-as.Date(df$EXAMDATE,format="%Y-%m-%d")) 
df <- df[order(df$RID,df$EXAMDATE,df$dt),]
df <- df[!duplicated(df[c("RID","EXAMDATE")]),] #Keep the demographic info collected closest in date to the EXAMDATE
df$Age=as.numeric(as.Date(df$EXAMDATE,format="%Y-%m-%d")-df$refDOB)/365 #Age at the time of df scan
rm(dem)
df=subset(df, select = -c(refDOB,dem.DATE,dt) )

dxsm=getDXSUM()
df <- dplyr::left_join(df,dxsm,by="RID")
df$dt=abs(as.Date(df$DX.DATE,format="%Y-%m-%d")-as.Date(df$EXAMDATE,format="%Y-%m-%d"))
df<-df[df$dt<=180,]#limiting to observations which have the other covariates collected within 180 days
df <- df[order(df$RID,df$EXAMDATE,df$dt),]
df <- df[!duplicated(df[c("RID","EXAMDATE")]),] #Keep the DXSM exam closest in date to the df EXAMDATE
df=subset(df, select = -c(DX.DATE) )
df=df[!is.na(df$DXSUM),]
rm(dxsm)

df$Age = as.numeric(df$Age)
df <- df %>%
  dplyr::filter(Age<85)

df$Sex<-relevel(factor(df$Sex),ref = 'M')
df$APOE<-relevel(factor(df$APOE),ref = 'E3/E3')
df$int_sex_apoe <- interaction(df$Sex, df$APOE) #interaction variable that allows teasing out details of individual groups
df = df[df$DXSUM %in% dx,]
df$DXSUM=factor(df$DXSUM)

df$Agediff=df$Age-min(df$Age)
# df$Edudiff=df$PTEDUCAT-16 #Since we have been predicting for 16 years of education

df <- df %>% dplyr::filter(APOE!="E2/E2")

df <- df %>% dplyr::mutate(
  reg = {{ response_var }}
)

df <- df %>%
  tidyr::drop_na(reg,Agediff,APOE,RID,PTEDUCAT,APOE,Sex)

df$NumOfScans=c(1)
df<-df[order(df$RID,df$EXAMDATE),]
df$AvgYrsBtnVisits=c(0)
df$TotalFollowupYrs=c(0)
for (id in unique(df$RID)){
  df$NumOfScans[df$RID==id]=length(df$RID[df$RID==id])
  if (length(df$RID[df$RID==id])>1){
    df$AvgYrsBtnVisits[df$RID==id]=mean(diff(as.numeric(as.Date(df$EXAMDATE[df$RID==id])))/365)
    df$TotalFollowupYrs[df$RID==id]=round(sum(diff(as.numeric(as.Date(df$EXAMDATE[df$RID==id])))/365),2)
  }
}

return(df %>%
         dplyr::select(-dt))}

test <- dataset_builder(wml,log_WMH)

tmpdd <- test

demog_info=rbind(demog_info,demog_distribution(tmpdd,"WMH"))

Participant_Records=merge(Participant_Records,data.frame(RID=tmpdd$RID[!duplicated(tmpdd$RID)],WMH=1),by=c("RID"),all = TRUE)

plotAPOEnSexinAgeBins(tmpdd,"WMH")

plt_labs <- labs(y =  "WMH",
                 x = 'Age in Years',
                 colour = 'APOE')

gamm_modeling <- function(df,
                          measure_name,
                          additional_vars=NULL,
                          additional_predicted_values = NULL, ## must be a list of desired fixed values in covariates for prediction if supplied
                          knots=NULL){
  
if(!is.null(knots)) {
formula_m1 <- paste0("reg ~ s(Agediff, bs='cs', k = ",knots,", by= APOE) +
    s(RID,Agediff, bs='re')  + s(RID,bs='re')  +
    PTEDUCAT + APOE + Sex")
  
formula_m2 <- paste0("reg ~ s(Agediff, bs = 'cs', k =", knots, ", by = interaction(APOE,Sex)) +
                      s(RID, Agediff, bs = 're') + s(RID,bs='re')  + 
                      PTEDUCAT + APOE + Sex + APOE:Sex")
}
  
if(is.null(knots)){
formula_m1 <- "reg ~ s(Agediff, bs='cs', by= APOE) +
    s(RID,Agediff, bs='re')  + s(RID,bs='re')  +
    PTEDUCAT + APOE + Sex"

formula_m2 <- "reg ~ s(Agediff, bs = 'cs', by = interaction(APOE,Sex)) +
                      s(RID, Agediff, bs = 're') + s(RID,bs='re')  + 
                      PTEDUCAT + APOE + Sex + APOE:Sex"
}
  
if(!is.null(additional_vars)){

formula_m1 <- paste(formula_m1,additional_vars,sep = " + ")

formula_m2 <- paste(formula_m2,additional_vars,sep = " + ")
}
  
formula_m1 <- as.formula(formula_m1)
  
formula_m2 <- as.formula(formula_m2)

m1 <- bam(formula = formula_m1, data = df,method = "REML") 

m1$model$Age=m1$model$Agediff+min(df$Age)

m2 <- bam(formula_m2,
          data = df,
          method = "REML")

m2$model$Age=m2$model$Agediff+min(df$Age)

## NEED TO FIX HARDCODED WMLICV 
m1_predict_values <- list(PTEDUCAT=16,RID=NULL,Sex='F')
if(!is.null(additional_predicted_values)){
m1_predict_values <- append(m1_predict_values,additional_predicted_values)
}

m2_predict_values <- list(PTEDUCAT=16,RID=NULL)
if(!is.null(additional_predicted_values)){
  m2_predict_values <- append(m2_predict_values,additional_predicted_values)
}

m1_predict<-predict_gam(m1,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=m1_predict_values) 
m1_predict$APOE=factor(m1_predict$APOE)
m1_predict$APOE = relevel(m1_predict$APOE,ref='E3/E3')
m1_predict$Age=m1_predict$Agediff+min(df$Age)#Getting back from age difference to actual age

m2_predict<-predict_gam(m2,exclude_terms = c('s(RID)','s(RID,Agediff)'),length_out = 100,
                        values=m2_predict_values)

m2_predict$APOE=factor(m2_predict$APOE)
m2_predict$APOE = relevel(m2_predict$APOE,ref='E3/E3')
m2_predict$Sex=factor(m2_predict$Sex)
m2_predict$Sex = relevel(m2_predict$Sex,'M')
m2_predict$Age=m2_predict$Agediff+min(df$Age)#Getting back from age difference to actual age

if(!is.null(additional_vars)){
covars_for_sig_or_not <- paste("PTEDUCAT",additional_vars,sep = ",")
} else {
  covars_for_sig_or_not <- "PTEDUCAT"
}

#Printing progression with significance info 
SigOrNot=CompareSplineSigSlope(measure_name,"PTEDUCATCovar",df,m1,m2,covars_for_sig_or_not,'F')

m1_predict$z=0
m1_predict$z= (m1_predict$fit-mean(df$reg[df$APOE=='E3/E3' & df$DXSUM=='CN']) )/stats::sd(df$reg[df$APOE=='E3/E3' & df$DXSUM=='CN'])
m1_predict$Significance=0
m1_predict$Significance[m1_predict$APOE=='E3/E3']=SigOrNot$deriv_e3$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E2']=SigOrNot$deriv_e2$Significance
m1_predict$Significance[m1_predict$APOE=='E3/E4']=SigOrNot$deriv_e4$Significance
m1_predict$Significance[m1_predict$APOE=='E4/E4']=SigOrNot$deriv_e4_homozygous$Significance
m1_predict$Measure=measure_name
m1_predict$SigAPOE=paste0(m1_predict$Significance,m1_predict$APOE)

m2_predict$z=0
m2_predict$z= (m2_predict$fit-mean(df$reg[df$APOE=='E3/E3' & df$DXSUM=='CN']) )/stats::sd(df$reg[df$APOE=='E3/E3' & df$DXSUM=='CN'])
m2_predict$Significance=0
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='M']=SigOrNot$deriv_me3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E3' & m2_predict$Sex=='F']=SigOrNot$deriv_fe3$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='M']=SigOrNot$deriv_me2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E2' & m2_predict$Sex=='F']=SigOrNot$deriv_fe2$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4$Significance
m2_predict$Significance[m2_predict$APOE=='E3/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4$Significance
m2_predict$Significance[m2_predict$APOE=='E4/E4' & m2_predict$Sex=='M']=SigOrNot$deriv_me4_homozygous$Significance
m2_predict$Significance[m2_predict$APOE=='E4/E4' & m2_predict$Sex=='F']=SigOrNot$deriv_fe4_homozygous$Significance
m2_predict$Measure= measure_name
m2_predict$SigAPOE=paste0(m2_predict$Significance,m2_predict$APOE)

prediction_plots_jpg_name <- paste0(measure_name,"PTEDUCATCovarModelPredictionPlots.jpg")

jpeg(prediction_plots_jpg_name,width=6000,height = 1500)
gg=plot_m1andm2predicts(m1_predict,m2_predict)
print(gg)
dev.off()

prediction_plots <<- plot_m1andm2predicts(m1_predict,m2_predict) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank())}

wml_icv_value <-mean(test$wmlICV,na.rm=T)

test_gamm <- gamm_modeling(df = test,
                           measure_name = "WMH",
                           additional_vars=c("wmlICV"),
                           additional_predicted_values = list(wmlICV=wml_icv_value),
                           knots=3)

spectrum_plot_across_sex_df_wmh <- spectrum_plot_across_sex_df
spectrum_plot_across_apoe_df_wmh <- spectrum_plot_across_apoe_df

wmh_between_apoe_diff_ggplot <-tmp_between_apoe_diff_ggplot
wmh_between_sex_diff_ggplot<- tmp_between_sex_diff_ggplot

wmh_prediction_plots <- prediction_plots

sink("ModelSummaryoutput.txt",append = TRUE)
print("WMH")
print("m1 APOE effect model summary")
print(summary(m1))
print("--------")
print("m2 APOExSex effect model summary")
print(summary(m2))
print("------------------------------")
sink()