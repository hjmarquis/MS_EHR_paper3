
library(tidyverse)
library(zoo)
library(openxlsx)
library(lubridate)

# Load --------------------------------------------------------------------
wkpath  = "D:/BoxUPitt/Box Sync/Boston/MS CLIMB Data/"
wkpath2 = "D:/BoxUPitt/Box Sync/Boston/MS CLIMB Paper2 Code&Result/Step3_BuildModel/IntermediateDataResult/"


MS_map  = read.xlsx(paste0(wkpath,"CLIMB Cohort/Spec95_i2b2_Mapping2017.xlsx"), sheet = 1)
colnames(MS_map) = c("PatientNum","PatientID")
MS_cohort = read.csv(paste0(wkpath2,"Cleaned_MS_cohort.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c(rep("character",4),rep("Date",4),rep("integer",6),
                                    rep("numeric",3),"integer",rep("Date",4),
                                    rep("numeric",3),"integer"))
MS_attack = read.csv(paste0(wkpath2,"Cleaned_MS_attack.csv"),
                     stringsAsFactors = FALSE,
                     colClasses = c("character",rep("integer",6),"Date",rep("numeric",3),
                                    rep("character",12),rep("Date",2),"integer",
                                    rep("Date",5),rep("integer",2)))
MS_attack$PatientNum <- MS_map$PatientNum[match(MS_attack$PatientID,MS_map$PatientID)]
MS_trt    = read.csv(paste0(wkpath2,"Cleaned_MS_treatment.csv"), stringsAsFactors = FALSE)

# New patients from Chart Review (of 200 unique patients, May 2020)

## RELAPSE DATA ##
group1 <- read.xlsx(paste0(wkpath, "Chart Review Annotation/CHANL_Group1_verified.xlsx"),sheet = 2, detectDates = TRUE)
group2 <- read.xlsx(paste0(wkpath, "Chart Review Annotation/CHANL_Group2_Verified.xlsx"),sheet = 3, detectDates = TRUE)
select2 <- read.xlsx(paste0(wkpath, "Chart Review Annotation/2patients.xlsx"),sheet = 3, detectDates = TRUE)
group1 <- group1 %>% dplyr::select(patient_number, relapse1_type,	relapse1_date)
group2 <- group2 %>% dplyr::select(patient_number, relapse1_type,	relapse1_date)
select2 <- select2 %>% dplyr::select(patient_number, relapse1_type,	relapse1_date)
select2$relapse1_date <- as.character(select2$relapse1_date)
relapse_new <- rbind(group1, group2, select2); rm(group1, group2, select2)
colnames(relapse_new)[1] <- "PatientNum"
relapse_new$PatientID <- MS_map$PatientID[match(relapse_new$PatientNum,MS_map$PatientNum)]
relapse_new <- relapse_new %>% dplyr::select(PatientID, PatientNum, everything())
# Remove "null" relapse entries
relapse_new <- relapse_new[!(relapse_new$relapse1_type=="null"),]
relapse_new$relapse1_date <- as.Date(relapse_new$relapse1_date)
# Fix 1 manual error for PatientNum 139534
relapse_new$relapse1_date[is.na(relapse_new$relapse1_date)] <- as.Date("2009-01-15")
# Define types of relapse: (according to Chart Review Annotation/CHANL_Data_Key)
relapse_new$clinical <- ifelse(relapse_new$relapse1_type %in% c(1,3), 1, 0)
relapse_new$radiographic <- ifelse(relapse_new$relapse1_type %in% c(2,3), 1, 0)
relapse_new$relapse1_type <- NULL
# All new patients belong to EHR
relapse_new$EHR <- 1

# Merge with older relapse data
MS_attack_new <- MS_attack %>% dplyr::select(PatientID, EHR, onset, clinical, radiographic)
MS_attack_new$PatientNum <- MS_map$PatientNum[match(MS_attack_new$PatientID, MS_map$PatientID)]
MS_attack_new <- MS_attack_new %>% dplyr::select(PatientID, PatientNum, EHR, onset, everything())
relapse_new <- relapse_new %>% dplyr::select(PatientID, PatientNum, EHR, everything())
names(MS_attack_new)[4] = names(relapse_new)[4] <- "relapse_date"

MS_attack_new <- rbind(MS_attack_new, relapse_new); rm(relapse_new)

## DMT DATA ##
group1 <- read.xlsx(paste0(wkpath, "Chart Review Annotation/CHANL_Group1_verified.xlsx"),sheet = 3, detectDates = TRUE)
group2 <- read.xlsx(paste0(wkpath, "Chart Review Annotation/CHANL_Group2_Verified.xlsx"),sheet = 2, detectDates = TRUE)
select2 <- select2 <- read.xlsx(paste0(wkpath, "Chart Review Annotation/2patients.xlsx"),sheet = 2, detectDates = TRUE)
group1 <- group1 %>% dplyr::select(patient_number, dmt1_type,	dmt1_type_cmt, dmt1_start_date, dmt1_end_date)
group2 <- group2 %>% dplyr::select(patient_number, dmt1_type,	dmt1_type_cmt, dmt1_start_date, dmt1_end_date)
select2 <- select2 %>% dplyr::select(patient_number, dmt1_type,	dmt1_type_cmt, dmt1_start_date, dmt1_end_date)
select2$dmt1_start_date <- as.character(select2$dmt1_start_date)
dmt_new <- rbind(group1, group2, select2); rm(group1, group2, select2)
colnames(dmt_new)[1] <- "PatientNum"
dmt_new$PatientID <- MS_map$PatientID[match(dmt_new$PatientNum,MS_map$PatientNum)]
dmt_new <- dmt_new %>% dplyr::select(PatientID, PatientNum, everything())
# Remove "null" DMT entries
dmt_new <- dmt_new %>% filter(dmt1_type != "null")
# Apply number to drug mapping
dmt_new$medication_desc <- sapply(dmt_new$dmt1_type, function(type) {
  if (type == 1) {
    return("Interferon-beta")
  } else if (type == 2) {
    return("Glatiramer acetate")
  } else if (type == 3) {
    return("Natalizumab")
  } else if (type == 4) {
    return("Dimethyl fumarate")
  } else if (type == 5) {
    return("Fingolimod")
  } else if (type == 6) {
    return("B-cell depleting")
  } else if (type == 7) {
    return("Teriflunomide")
  } else if (type == 99) {
    # Other
    return(NA)
  }
})
# Replace "other" DMT name with description
dmt_new$medication_desc[is.na(dmt_new$medication_desc)] <- dmt_new$dmt1_type_cmt[is.na(dmt_new$medication_desc)]
# Fix 1 manual error for PatientNum 119962
dmt_new$medication_desc[dmt_new$medication_desc == "null"] <- "Cytoxan"
# Standardize names to capitalized
dmt_new$medication_desc <- toupper(dmt_new$medication_desc)
dmt_new <- dmt_new %>% dplyr::select(PatientID, PatientNum, medication_desc, dmt1_start_date, dmt1_end_date)
colnames(dmt_new)[c(4,5)] <- c("start_date", "stop_date")
# Change stop date "null" to NA
dmt_new$stop_date[dmt_new$stop_date == "null"] <- NA
# Convert to date
dmt_new$start_date <- as.Date(dmt_new$start_date)
dmt_new$stop_date <- as.Date(dmt_new$stop_date)
# Consider all dates "validated" due to manual chart review
dmt_new$val_start <- dmt_new$start_date
dmt_new$val_stop <- dmt_new$stop_date


MS_trt$PatientNum <- MS_map$PatientNum[match(MS_trt$PatientID,MS_map$PatientID)]
MS_trt$start_date <- with(MS_trt, as.Date(paste(start_year, start_month, start_day, sep = "-")))
MS_trt$stop_date <- with(MS_trt, as.Date(paste(stop_year, stop_month, stop_day, sep = "-")))

MS_trt_new <- MS_trt %>% 
  dplyr::select(PatientID, PatientNum, medication_desc, start_date, stop_date, val_start, val_stop)
# Merge
MS_trt_new <- rbind(MS_trt_new, dmt_new); rm(dmt_new)


############################################
############### EHR Data ###################
############################################
MS_EHR    = read.xlsx(paste0(wkpath,"EHR/MS_first_last_date_whole_list.xlsx"),sheet = 1)
colnames(MS_EHR)[1] = "PatientNum"
MS_EHR$First_date <- convertToDate(MS_EHR$First_date)
MS_EHR$Last_date <- convertToDate(MS_EHR$Last_date)
MS_EHR$PatientID = MS_map$PatientID[match(MS_EHR$PatientNum,MS_map$PatientNum)] # Add CLIMB PatientID
# Demographics data
MS_demographics <- read.xlsx(paste0(wkpath, "EHR/ms_demographics.xlsx"))
MS_demographics <- MS_demographics %>% dplyr::select(PATIENT_NUM, BIRTH_DATE, SEX_CD, RACE_CD)
MS_demographics$BIRTH_DATE <- convertToDate(MS_demographics$BIRTH_DATE)
# NLP Note data:each row each record; each patient multiple records; contains CUI, used for mentions of MS, relapses,.. 
MS_NLP    = read.csv(paste0(wkpath,"EHR/NLP_note_level_using_code_from_zq_raw_DocTimeRel_excluded_20180604.csv"),
                     header=FALSE,skip=1,sep = "|",stringsAsFactors = FALSE)
MS_NLP    = MS_NLP[,c(1:3,6)]
colnames(MS_NLP)  = c("PatientNum","EncounterNum","CUI","StartDate")
MS_NLP$StartDate  = as.Date(MS_NLP$StartDate,"%Y-%m-%d %H:%M:%S")
# CUI Dictionary
CUIdictAll = read.xlsx(paste0(wkpath,"EHR/AllCUI_Database.xlsx"), sheet = 1)
colnames(CUIdictAll) = c("ConceptCd","Desc")

## ICD codes
ICDPheCode <- read.csv(paste0(wkpath, "EHR/MS_EHR_ICDCodes_05142020.csv"), header = FALSE,
                       col.names = c("patient_num", "encounter_num","start_date","concept_cd","phecode","phecode_description"),
                       stringsAsFactors = FALSE)
ICDPheCode$start_date <- as.Date(ICDPheCode$start_date)
ICDPheCode$phecode[ICDPheCode$concept_cd == "LPA268"] = "335_"
ICDPheCode$phecode_description[ICDPheCode$concept_cd == "LPA268"] = "Multiple sclerosis"
tmp = ICDPheCode[ICDPheCode$concept_cd%in%c("V700","V762"),-6]
ICDPheCode  = ICDPheCode[!ICDPheCode$concept_cd%in%c("V700","V762"),]
colnames(ICDPheCode)  = c("PatientNum","EncounterNum","StartDate","ConceptCd","PheCode","Desc")
## CPT codes
CPTGrouped  = read.csv(paste0(wkpath,"EHR/MS_EHR_CPTCodes_05122020.csv"),stringsAsFactors = FALSE,
                       col.names = c("patient_num","encounter_num", "start_date","concept_cd", "cpt_group"))
CPTGrouped$start_date <- as.Date(CPTGrouped$start_date)
colnames(CPTGrouped)  = c("PatientNum","EncounterNum","StartDate","ConceptCd","CPTGroup")
colnames(tmp) = colnames(CPTGrouped)
CPTGrouped  = rbind(tmp, CPTGrouped)
CPTGrouped$ConceptCd = toupper(CPTGrouped$ConceptCd)
# CPTGrouped$CPTGroup  = substr(CPTGrouped$CPTGroup,2,nchar(CPTGrouped$CPTGroup))
CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("V700","CG0463")] = "outpatient clinic visit"
CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("V762","CG0101","CG0202","CQ0091")] ="cancer screening"
CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C90765","C90766","C90780","C90781","CC8950","CC8951")] = "hydration, therapeutic, prophylactic, diagnostic injections and infusions, and chemotherapy and other highly complex drug or highly complex biologic agent administration"
CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C70540","C70542","C70543")] = "MRI orbit"
CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C70551","C70552","C70553")] = "MRI brain"
CPTGrouped$CPTGroup[CPTGrouped$ConceptCd%in%c("C72141","C72142","C72146","C72147","C72156","C72157")] = "MRI spine (cervical, thoracic)"

## CUI codes
CUISelected  = read.csv(paste0(wkpath,"EHR/MS_AllEncounters_CUIs_08072019.csv"), stringsAsFactors = FALSE, header = FALSE,col.names = c('patient_num', 'encounter_num', 'start_date','concept_cd'))
CUISelected$start_date <- format(as.POSIXct(CUISelected$start_date, format = '%Y-%m-%d %H:%M:%S'), format = '%Y-%m-%d')
CUISelected$start_date <- as.Date(CUISelected$start_date)
CUISelected <- CUISelected %>% filter(start_date >= as.Date("1980-01-01"))
colnames(CUISelected) = c("PatientNum","EncounterNum","StartDate","ConceptCd")

ICDPheCode$ConceptCd  = paste0("PheCode.",ICDPheCode$PheCode)
CPTGrouped$ConceptCd  = paste0("CPTGroup.",CPTGrouped$CPTGroup)
CUISelected$ConceptCd = paste0("CUI.",CUISelected$ConceptCd)

# Combine
ICD_CPT_CUI_Comb = rbind(ICDPheCode[,c("PatientNum","EncounterNum","StartDate","ConceptCd")],
                         CPTGrouped[,c("PatientNum","EncounterNum","StartDate","ConceptCd")],
                         CUISelected[,c("PatientNum","EncounterNum","StartDate","ConceptCd")])
ICD_CPT_CUI_Comb = plyr::count(ICD_CPT_CUI_Comb)[,c("PatientNum", "StartDate", "ConceptCd")]
ICD_CPT_CUI_Comb$PatientID = MS_map$PatientID[match(ICD_CPT_CUI_Comb$PatientNum,MS_map$PatientNum)] 
ICD_CPT_CUI_Comb <- ICD_CPT_CUI_Comb %>% dplyr::select(PatientID, PatientNum, everything())
ICD_CPT_CUI_Comb$ConceptCd = as.factor(ICD_CPT_CUI_Comb$ConceptCd)


## RXNORM data
# Target drugs: RTX, NTZ, DMF, FGL
RXNORM_comb <- read.xlsx(paste0(wkpath,"MS DMT/MS_Medication_updatedData_05052020.xlsx"), sheet = 2)
RXNORM_comb$start_date <- convertToDate(RXNORM_comb$start_date)
RXNORM_comb$RX_Generic_Name <- sapply(RXNORM_comb$RX_Norm_code, function(code) {
  if (code == 121191) {return("Rituximab")}
  if (code == 1373478) {return("Dimethyl Fumarate")}
  if (code == 354770) {return("Natalizumab")}
  if (code == 1012892) {return("Fingolimod")}
})
RXNORM_comb$CLIMBID <- MS_map$PatientID[match(RXNORM_comb$patient_num,MS_map$PatientNum)]
RXNORM_comb <- RXNORM_comb %>% dplyr::select(patient_num, CLIMBID, everything())
# Interferon 
RXNORM_interferon <- read.xlsx(paste0(wkpath,"MS DMT/MS_AdditionalRXNorm_data_afterUpdate_05072020.xlsx"), sheet = 7)
RXNORM_interferon$start_date <- convertToDate(RXNORM_interferon$start_date)
# Additional DMTs
RXNORM_addl_drugs <- read.xlsx(paste0(wkpath,"MS DMT/MS_AdditionalRXNorm_data_afterUpdate_05072020.xlsx"), sheet = 5)
RXNORM_addl_drugs$start_date <- convertToDate(RXNORM_addl_drugs$start_date)
RXNORM_addl_drugs <- RXNORM_addl_drugs %>% filter(RX_Generic_Name != "Mitoxantrone") # Avoid duplication, as mitoxantrone is included in chemotherapy file below
# Chemotherapy drugs
RXNORM_chemo <- read.xlsx(paste0(wkpath, "MS DMT/MS_EHR_RXNORM_06042020.xlsx"), sheet = 3)
RXNORM_chemo$start_date <- convertToDate(RXNORM_chemo$start_date)
# Combine
RXNORM_comb <- rbind(RXNORM_comb, RXNORM_interferon, RXNORM_addl_drugs, RXNORM_chemo)
rm(RXNORM_interferon, RXNORM_addl_drugs, RXNORM_chemo)
## Steroid, hospitalization, and MRI CPT codes
CPT_list  = read.csv(paste0(wkpath,"EHR/CPT_codelist_of_interest.csv"),header=FALSE,skip=1,sep = "|",stringsAsFactors = FALSE)
colnames(CPT_list) = c("CPT_code","Category","Description")
CPT_list$Category = sapply(CPT_list$Category,function(x){
  if(x=="ED visit"){
    "ED"
  } else if(x=="Hospital Discharge"){
    "HD"
  } else if(x=="Solumedrol (methylprednisolone)"){
    "Solu"
  } else{
    "MRI"
  }
})
CPT_list <- CPT_list[,-3]

MS_CPT <- read.csv(paste0(wkpath, "EHR/MS_selectedCPT_05202020.csv"), header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("PatientNum", "EncounterNum", "StartDate", "Concept_cd", "Desc"),
                   colClasses = c("integer", "integer", "Date", "character", "character"),
                   fileEncoding="UTF-8-BOM")
MS_CPT$Concept_cd     = toupper(MS_CPT$Concept_cd)


######################################### 
######## Subset post-2006 data ##########
#########################################
ICD_CPT_CUI_Comb <- ICD_CPT_CUI_Comb %>% filter(StartDate >= as.Date("2006-01-01"))
RXNORM_comb <- RXNORM_comb %>% filter(start_date >= as.Date("2006-01-01"))
MS_attack_new <- MS_attack_new %>% filter(relapse_date >= as.Date("2006-01-01"))
MS_NLP <- MS_NLP %>% filter(StartDate >= as.Date("2006-01-01"))
MS_CPT <- MS_CPT %>% filter(StartDate >= as.Date("2006-01-01"))

# Subset codes data into code type
ICD_subset <- ICD_CPT_CUI_Comb[str_detect(ICD_CPT_CUI_Comb$ConceptCd, "PheCode"),]
CPT_subset <- ICD_CPT_CUI_Comb[str_detect(ICD_CPT_CUI_Comb$ConceptCd, "CPT"),]
CUI_subset <- ICD_CPT_CUI_Comb[str_detect(ICD_CPT_CUI_Comb$ConceptCd, "CUI"),]
# Only subset patients with EHR data, post-2006 data
MS_trt_new_EHR <- MS_trt_new %>% filter(!is.na(PatientNum) & start_date >= as.Date("2006-01-01")) 

save(list = objects(), file = "MS CLIME analysis/data-processed/common.rda")