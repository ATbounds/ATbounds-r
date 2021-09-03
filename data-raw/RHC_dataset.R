## code to prepare Right Heart Catheterization (RHC) Study dataset

rm(list = ls())

rhc <- read.csv("https://hbiostat.org/data/repo/rhc.csv")

# Demographics (3 variables)
reg_demog <- c("age","sex","edu")
# Categorical variables with multiple levels (20 variables)
reg_multi <- c("cat1","ca","ninsclas","race","income")
# Categories of admission diagnosis (10 variables)
reg_diagn <- c("resp","card","neuro","gastr","renal","meta","hema","seps","trauma","ortho")
# Categories of comorbidities illness (12 variables)
reg_comor <- c("cardiohx","chfhx","dementhx","psychhx","chrpulhx",
               "renalhx","liverhx","gibledhx","malighx","immunhx",
               "transhx","amihx")
# Other regressors (20 variables)
reg_other <- c("das2d3pc","dnr1","surv2md1","aps1","scoma1",
               "wtkilo1","temp1","meanbp1","resp1","hrt1",
               "pafi1","paco21","ph1","wblc1","hema1",
               "sod1","pot1","crea1","bili1","alb1")
# Outcome variable: 1 if a patient survives after 30 days of admission; 0 if a patient dies within 30 days
y_rhc <- 2-as.integer(rhc[,"death"])
# Treatment variable: 1 if RHC; 0 if no RHC 
t_rhc <- as.integer(rhc[,"swang1"])-1

# construction of covariates
all_reg_names <- c(reg_demog,reg_multi,reg_diagn,reg_comor,reg_other)
rhc_covariates <- rhc[,all_reg_names]

x_rhc <- fastDummies::dummy_cols(rhc_covariates,
                                 select_columns=c(reg_multi,reg_diagn,"dnr1","sex"),
                                 remove_most_frequent_dummy=TRUE,
                                 remove_selected_columns=TRUE)

# construction of wt0 (1 if weight == 0, and 0 otherwise)
wt0 <- as.integer(rhc[,"wtkilo1"]==0)

# construction of cat2 dummies
rhc_cat2 <- rhc[,"cat2"]
rhc_cat2_dummies <- fastDummies::dummy_cols(rhc_cat2)
# Remove the original cat2 variable and the dummy for NA
rhc_cat2_dummies <- rhc_cat2_dummies[,2:7]   
# Replace NA with 0
rhc_cat2_dummies[is.na(rhc_cat2_dummies)]<-0

# 72 covariates are constructed following Hirano and Imbens (2001)

x_rhc <- cbind(x_rhc,rhc_cat2_dummies,wt0)

reg_rhc <- {}
for (j in 1:ncol(x_rhc)){
  reg_rhc <- cbind(reg_rhc, as.numeric(x_rhc[,j]))
}

colnames(reg_rhc) <- colnames(x_rhc)

RHC <- cbind(y_rhc,t_rhc,reg_rhc)

colnames(RHC)[1]<-"survival"
colnames(RHC)[2]<-"RHC"
colnames(RHC)[38]<-"cat1_Colon_Cancer"
colnames(RHC)[41]<-"cat1_Lung_Cancer"
colnames(RHC)[42]<-"cat1_MOSF_Malignancy"
colnames(RHC)[43]<-"cat1_MOSF_Sepsis"
colnames(RHC)[48]<-"ninsclas_Medicare_and_Medicaid"
colnames(RHC)[49]<-"ninsclas_No_insurance"
colnames(RHC)[50]<-"ninsclas_Private_and_Medicare"
colnames(RHC)[53]<-"income3"
colnames(RHC)[54]<-"income1"
colnames(RHC)[55]<-"income2"
colnames(RHC)[68]<-"cat2_Cirrhosis"
colnames(RHC)[69]<-"cat2_Colon_Cancer"
colnames(RHC)[70]<-"cat2_Coma"
colnames(RHC)[71]<-"cat2_Lung_Cancer"
colnames(RHC)[72]<-"cat2_MOSF_Malignancy"
colnames(RHC)[73]<-"cat2_MOSF_Sepsis"

RHC <- as.data.frame(RHC)

usethis::use_data(RHC, overwrite = TRUE)
