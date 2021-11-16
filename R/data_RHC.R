#' RHC
#'
#' The right heart catheterization (RHC) dataset is publicly available on the Vanderbilt Biostatistics website.
#' RHC is a diagnostic procedure for directly measuring cardiac function in critically ill patients. 
#' The dependent variable is 1 if a patient survived after 30 days of admission, 0 if a patient died within 30 days.
#' The treatment variable is 1 if RHC was applied within 24 hours of admission, and 0 otherwise.
#' The sample size was n = 5735, and 2184 patients were treated with RHC. 
#' Connors et al. (1996) used a propensity score matching approach to study the efficacy of RHC,
#' using data from the observational study called SUPPORT (Murphy and Cluff, 1990). 
#' Many authors used this dataset subsequently.
#' The 72 covariates are constructed, following Hirano and Imbens (2001).
#' 
#' @references Connors, A.F., Speroff, T., Dawson, N.V., Thomas, C., Harrell, F.E., Wagner, D., Desbiens, N., Goldman, L., Wu, A.W., Califf, R.M. and Fulkerson, W.J., 1996. 
#' The effectiveness of right heart catheterization in the initial care of critically III patients. JAMA, 276(11), pp.889-897.
#' \doi{10.1001/jama.1996.03540110043030}
#'
#' @references Hirano, K., Imbens, G.W. Estimation of Causal Effects using Propensity Score Weighting: An Application to Data on Right Heart Catheterization, 2001.
#' Health Services & Outcomes Research Methodology 2, pp.259–278. 
#' \doi{10.1023/A:1020371312283}
#' 
#' @references D. J. Murphy, L. E. Cluff, SUPPORT: Study to understand prognoses and preferences for outcomes and risks of treatments—study design, 1990.
#' Journal of Clinical Epidemiology, 43, pp. 1S–123S
#' \url{https://www.jclinepi.com/issue/S0895-4356(00)X0189-8}
#' .
#' @format A data frame with 5735 rows and 74 variables:
#' \describe{
#' \item{survival}{Outcome: 1 if a patient survived after 30 days of admission, and 0 if a patient died within 30 days}
#' \item{RHC}{Treatment: 1 if RHC was applied within 24 hours of admission, and 0 otherwise.}
#' \item{age}{Age in years}
#' \item{edu}{Years of education}
#' \item{cardiohx}{Cardiovascular symptoms}
#' \item{chfhx}{Congestive Heart Failure}
#' \item{dementhx}{Dementia, stroke or cerebral infarct, Parkinson’s disease}
#' \item{psychhx}{Psychiatric history, active psychosis or severe depression}
#' \item{chrpulhx}{Chronic pulmonary disease, severe pulmonary disease}
#' \item{renalhx}{Chronic renal disease, chronic hemodialysis or peritoneal dialysis}
#' \item{liverhx}{Cirrhosis, hepatic failure}
#' \item{gibledhx}{Upper GI bleeding}
#' \item{malighx}{Solid tumor, metastatic disease, chronic leukemia/myeloma, acute leukemia, lymphoma}
#' \item{immunhx}{Immunosuppression, organ transplant, HIV, Diabetes Mellitus, Connective Tissue Disease}
#' \item{transhx}{transfer (> 24 hours) from another hospital}
#' \item{amihx}{Definite myocardial infarction}
#' \item{das2d3pc}{DASI - Duke Activity Status Index}
#' \item{surv2md1}{Estimate of prob. of surviving 2 months}
#' \item{aps1}{APACHE score}
#' \item{scoma1}{Glasgow coma score}
#' \item{wtkilo1}{Weight}
#' \item{temp1}{Temperature}
#' \item{meanbp1}{Mean Blood Pressure}
#' \item{resp1}{Respiratory Rate}
#' \item{hrt1}{Heart Rate}
#' \item{pafi1}{PaO2/FI02 ratio}
#' \item{paco21}{PaCO2}
#' \item{ph1}{PH}
#' \item{wblc1}{WBC}
#' \item{hema1}{Hematocrit}
#' \item{sod1}{Sodium}
#' \item{pot1}{Potassium}
#' \item{crea1}{Creatinine}
#' \item{bili1}{Bilirubin}
#' \item{alb1}{Albumin}
#' \item{cat1_CHF}{1 if the primary disease category is CHF, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_Cirrhosis}{1 if the primary disease category is Cirrhosis, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_Colon_Cancer}{1 if the primary disease category is Colon Cancer, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_Coma}{1 if the primary disease category is Coma, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_COPD}{1 if the primary disease category is COPD, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_Lung_Cancer}{1 if the primary disease category is Lung Cancer, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_MOSF_Malignancy}{1 if the primary disease category is MOSF w/Malignancy, and 0 otherwise (Omitted category = ARF).}
#' \item{cat1_MOSF_Sepsis}{1 if the primary disease category is MOSF w/Sepsis, and 0 otherwise (Omitted category = ARF).}
#' \item{ca_Metastatic}{1 if cancer is metastatic, and 0 otherwise (Omitted category = no cancer).}
#' \item{ca_Yes}{1 if cancer is localized, and 0 otherwise (Omitted category = no cancer).}
#' \item{ninsclas_Medicaid}{1 if medical insurance category is Medicaid, and 0 otherwise (Omitted category = Private).}
#' \item{ninsclas_Medicare}{1 if medical insurance category is Medicare, and 0 otherwise (Omitted category = Private).}
#' \item{ninsclas_Medicare_and_Medicaid}{1 if medical insurance category is Medicare & Medicaid, and 0 otherwise (Omitted category = Private).}
#' \item{ninsclas_No_insurance}{1 if medical insurance category is No Insurance, and 0 otherwise (Omitted category = Private).}
#' \item{ninsclas_Private_and_Medicare}{1 if medical insurance category is Private & Medicare, and 0 otherwise (Omitted category = Private).}
#' \item{race_black}{1 if Black, and 0 otherwise (Omitted category = White).}
#' \item{race_other}{1 if Other, and 0 otherwise (Omitted category = White).}
#' \item{income3}{1 if Income >$50k, and 0 otherwise (Omitted category = under $11k).}
#' \item{income1}{1 if Income $11–$25k, and 0 otherwise (Omitted category = under $11k).}
#' \item{income2}{1 if Income $25–$50k, and 0 otherwise (Omitted category = under $11k).}
#' \item{resp_Yes}{Respiratory diagnosis}
#' \item{card_Yes}{Cardiovascular diagnosis}
#' \item{neuro_Yes}{Neurological diagnosis}
#' \item{gastr_Yes}{Gastrointestinal diagnosis}
#' \item{renal_Yes}{Renal diagnosis}
#' \item{meta_Yes}{Metabolic diagnosis}
#' \item{hema_Yes}{Hematological diagnosis}
#' \item{seps_Yes}{Sepsis diagnosis}
#' \item{trauma_Yes}{Trauma diagnosis}
#' \item{ortho_Yes}{Orthopedic diagnosis}
#' \item{dnr1_Yes}{Do Not Resuscitate status on day 1}
#' \item{sex_Female}{Female}
#' \item{cat2_Cirrhosis}{1 if the secondary disease category is Cirrhosis, and 0 otherwise (Omitted category = NA).}
#' \item{cat2_Colon_Cancer}{1 if secondary disease category is Colon Cancer, and 0 otherwise (Omitted category = NA).}
#' \item{cat2_Coma}{1 if the secondary disease category is Coma, and 0 otherwise (Omitted category = NA).}
#' \item{cat2_Lung_Cancer}{1 if the secondary disease category is Lung Cancer, and 0 otherwise (Omitted category = NA).}
#' \item{cat2_MOSF_Malignancy}{1 if the secondary disease category is MOSF w/Malignancy, and 0 otherwise (Omitted category = NA).}
#' \item{cat2_MOSF_Sepsis}{1 if the secondary disease category is MOSF w/Sepsis, and 0 otherwise (Omitted category = NA).}
#' \item{wt0}{weight = 0 (missing)}
#' }
#' 
#' @source The dataset is publicly available on the Vanderbilt Biostatistics website at 
#' \url{https://hbiostat.org/data/}.
"RHC"
