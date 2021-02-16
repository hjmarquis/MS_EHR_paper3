
trt_subset <- trt_subset %>% filter(start_date <= as.Date("2016-12-31"))
##########################################################################################
######################                  Add covariates               #####################
##########################################################################################
# Sex
trt_subset$FEMALE <- ifelse(trt_subset$PatientNum %in% MS_demographics$PATIENT_NUM[MS_demographics$SEX_CD == 'F'], 1, 0)

# Race
trt_subset$RACE <- ifelse(trt_subset$PatientNum %in% MS_demographics$PATIENT_NUM[MS_demographics$RACE_CD == 'WHITE'], 1, 0)

# Age at first MS code (MS ICD: 'PheCode.335_') (in years)
trt_subset$FIRSTMSICD_DATE <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICD_subset %>% filter(PatientNum == pnum & ConceptCd == 'PheCode.335_')
  tmp <- tmp[order(tmp$StartDate), ]
  return(as.character(tmp$StartDate[1]))
}) %>% as.Date()
trt_subset$BIRTHDAY <- sapply(trt_subset$PatientNum, function(p_num) {
  MS_demographics$BIRTH_DATE[MS_demographics$PATIENT_NUM == p_num] %>% as.character()
}) %>% as.Date()
trt_subset$AGE_AT_FIRSTMSICD <- with(trt_subset, difftime(FIRSTMSICD_DATE,BIRTHDAY,units='weeks')/52.25) %>% as.numeric()

# Follow-up duration: (in years): time period between date of occurrence of first of any ICD code in EHR and time of treatment initiation
trt_subset$FIRSTICD_DATE <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- ICD_subset %>% filter(PatientNum == pnum)
  tmp <- tmp[order(tmp$StartDate),]
  return(tmp$StartDate[1])
}) %>% as.Date()
trt_subset$FOLLOWUP_DURA <- with(trt_subset, as.numeric(difftime(start_date, FIRSTICD_DATE, units='weeks')/52.25))

# Disease duration proxy: Age at the time of DMT start - age at first MS ICD code (in years)
trt_subset$DISEASE_DURA <- 
  with(trt_subset, as.numeric(difftime(start_date, BIRTHDAY)/365)) - trt_subset$AGE_AT_FIRSTMSICD

# Reorder columns
trt_subset <- trt_subset %>% dplyr::select(PatientNum, PatientID, medication_desc, start_date, FIRSTMSICD_DATE,
                                           FIRSTICD_DATE, BIRTHDAY, everything())

# Health utilization (total and 3 months) prior to treatment initiation: number of all ICD codes + notes 
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  # Overall 
  tmp <- ICD_subset %>% filter(PatientNum == pnum & StartDate < tx_start_date)
  tmp2 <- MS_NLP %>% filter(PatientNum == pnum)
  # 3 months prior to DMT start
  tmp3 <- ICD_subset %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(3))
  tmp4 <- MS_NLP %>% 
    filter(PatientNum == pnum & StartDate >= tx_start_date - months(3))
  
  
  return(c(nrow(tmp) + length(unique(tmp2$StartDate)),
           nrow(tmp3) + length(unique(tmp4$StartDate))))
})
trt_subset[,c("HUTIL_OVERALL", "HUTIL_3MONS")] <- t(tmp)
# Fix NAs created by dividing by 0
trt_subset$HUTIL_OVERALL[is.na(trt_subset$HUTIL_OVERALL)] <- 0
trt_subset$HUTIL_3MONS[is.na(trt_subset$HUTIL_3MONS)] <- 0

# Adjusted MS frequency of ICD code for MS (total and 3 months) prior to treatment initiation: Number of ICD codes for MS divided by total number of ICD codes (of any kind)   
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  # Overall
  tmp <- ICD_subset %>% filter(PatientNum == pnum & StartDate < tx_start_date)
  tmp2 <- tmp %>% filter(ConceptCd == "PheCode.335_")
  # 3 months prior to DMT start
  tmp3 <- ICD_subset %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(3))
  tmp4 <- tmp3 %>% filter(ConceptCd == "PheCode.335_")
  
  return(c(nrow(tmp2)/nrow(tmp),
           nrow(tmp4)/nrow(tmp3)))
})
trt_subset[,c("ADJ_MSICD_OVERALL", "ADJ_MSICD_3MONS")] <- t(tmp)
trt_subset$ADJ_MSICD_OVERALL[is.na(trt_subset$ADJ_MSICD_OVERALL)] <- 0
trt_subset$ADJ_MSICD_3MONS[is.na(trt_subset$ADJ_MSICD_3MONS)] <- 0

# Adjusted frequency of MS CUI code (3 months) prior to treatment initiation: Number of MS CUIs (C0026769 & C0751967), divided by health utilization
trt_subset$ADJ_MSCUI_3MONS <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  health_util <- trt_subset$HUTIL_3MONS[trt_subset$PatientNum == pnum]
  tmp <- CUI_subset %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(3)) %>% 
    filter(ConceptCd == "CUI.C0026769" | ConceptCd == "CUI.C0751967")
  return(nrow(tmp)/health_util)
})
trt_subset$ADJ_MSCUI_3MONS[is.na(trt_subset$ADJ_MSCUI_3MONS)] <- 0

# Adjusted steroid (total and 6 months), prior to treatment = Number of steroid prescriptions divided by health utilization (based on CPT codes)
steroid_codes <- CPT_list$CPT_code[CPT_list$Category == "Solu"]
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  health_util <- trt_subset$HUTIL_OVERALL[trt_subset$PatientNum == pnum]
  # Overall
  tmp <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date) %>% 
    filter(Concept_cd %in% steroid_codes)
  # 6 months prior to DMT start
  tmp2 <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(6)) %>% 
    filter(Concept_cd %in% steroid_codes)
  # Calculate 6 month health utilization metric
  tmp3 <- ICD_subset %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(6))
  tmp4 <- MS_NLP %>% 
    filter(PatientNum == pnum & StartDate >= tx_start_date - months(6))
  tmp5 <- nrow(tmp3) + length(unique(tmp4$StartDate)) 
  
  return(c(nrow(tmp)/health_util,
           nrow(tmp2)/tmp5))
})
trt_subset[,c("ADJ_STEROID_OVERALL", "ADJ_STEROID_6MONS")] <- t(tmp)
trt_subset$ADJ_STEROID_OVERALL[is.na(trt_subset$ADJ_STEROID_OVERALL)] <- 0
trt_subset$ADJ_STEROID_6MONS[is.na(trt_subset$ADJ_STEROID_6MONS)] <- 0

# Adjusted MRI usage (total and 6 months), prior to treatment = Number of MRI brain / spine (based on CPT codes) divided by health utilization
mri_codes <- CPT_list$CPT_code[CPT_list$Category == "MRI"]
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  health_util <- trt_subset$HUTIL_OVERALL[trt_subset$PatientNum == pnum]
  # Overall
  tmp <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date) %>% 
    filter(Concept_cd %in% mri_codes)
  # 6 months prior to DMT start
  tmp2 <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(6)) %>% 
    filter(Concept_cd %in% mri_codes)
  # Calculate 6 month health utilization metric
  tmp3 <- ICD_subset %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(6))
  tmp4 <- MS_NLP %>% 
    filter(PatientNum == pnum & StartDate >= tx_start_date - months(6))
  tmp5 <- nrow(tmp3) + length(unique(tmp4$StartDate)) 
  
  return(c(nrow(tmp)/health_util,
           nrow(tmp2)/tmp5))
})
trt_subset[,c("ADJ_MRI_OVERALL", "ADJ_MRI_6MONS")] <- t(tmp)
trt_subset$ADJ_MRI_OVERALL[is.na(trt_subset$ADJ_MRI_OVERALL)] <- 0
trt_subset$ADJ_MRI_6MONS[is.na(trt_subset$ADJ_MRI_6MONS)] <- 0

# Adjusted hospitalization (total and 3 months) prior to treatment initiation = Number of hospitalizations (based on CPT codes) divided by health utilization 
hosp_codes <- CPT_list$CPT_code[CPT_list$Category == "HD"]
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  health_util <- trt_subset$HUTIL_OVERALL[trt_subset$PatientNum == pnum]
  health_util3mons <- trt_subset$HUTIL_3MONS[trt_subset$PatientNum == pnum]
  # Overall
  tmp <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date) %>% 
    filter(Concept_cd %in% hosp_codes)
  # 3 months prior to DMT start
  tmp2 <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(3)) %>% 
    filter(Concept_cd %in% hosp_codes)
  
  return(c(nrow(tmp)/health_util,
           nrow(tmp2)/health_util3mons))
})
trt_subset[,c("ADJ_HOSP_OVERALL", "ADJ_HOSP_3MONS")] <- t(tmp)
trt_subset$ADJ_HOSP_OVERALL[is.na(trt_subset$ADJ_HOSP_OVERALL)] <- 0
trt_subset$ADJ_HOSP_3MONS[is.na(trt_subset$ADJ_HOSP_3MONS)] <- 0

# Adjusted emergency room visits (total and 3 months) prior to treatment initiation = Number of emergency room visits (based on CPT codes) divided by health utilization
ED_codes <- CPT_list$CPT_code[CPT_list$Category == "ED"]
tmp <- sapply(trt_subset$PatientNum, function(pnum) {
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  health_util <- trt_subset$HUTIL_OVERALL[trt_subset$PatientNum == pnum]
  health_util3mons <- trt_subset$HUTIL_3MONS[trt_subset$PatientNum == pnum]
  # Overall
  tmp <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date) %>% 
    filter(Concept_cd %in% ED_codes)
  # 3 months prior to DMT start
  tmp2 <- MS_CPT %>% 
    filter(PatientNum == pnum & StartDate < tx_start_date & StartDate >= tx_start_date - months(3)) %>% 
    filter(Concept_cd %in% ED_codes)
  
  return(c(nrow(tmp)/health_util,
           nrow(tmp2)/health_util3mons))
})
trt_subset[,c("ADJ_ED_OVERALL", "ADJ_ED_3MONS")] <- t(tmp)
trt_subset$ADJ_ED_OVERALL[is.na(trt_subset$ADJ_ED_OVERALL)] <- 0
trt_subset$ADJ_ED_3MONS[is.na(trt_subset$ADJ_ED_3MONS)] <- 0

# Duration of prior DMT treatment (in months)
trt_subset$PRIORDMT_DURA <- sapply(trt_subset$PatientNum, function(pnum) {
  tmp <- MS_trt_new_EHR %>% 
    filter(PatientNum == pnum & start_date < with(trt_subset, start_date[PatientNum == pnum])) 
  tmp <- tmp[order(tmp$start_date),]
  if (nrow(tmp) == 0) {return(0)} else {
    tmp2 <- difftime(tmp$stop_date, tmp$start_date, units = "weeks")/4
    return(sum(as.numeric(tmp2), na.rm = TRUE))
  }
})

# Number of relapses in 12 months prior to treatment initiation
trt_subset$PRIOR_RELAPSE_12MONS <- sapply(1:nrow(trt_subset), function(i) {
  aa = MS_attack_new[MS_attack_new$PatientNum==trt_subset$PatientNum[i],]
  bb = difftime(aa$relapse_date,trt_subset$start_date[i],units = "days")/30 # in months
  CC = abs(bb) <= 12 & bb <=0
  # Obtain relapse within designated time period
  valid_date = aa[CC,]
  # Return unique relapse dates
  return(length(with(valid_date, unique(relapse_date[clinical >= 1 | radiographic >= 1]))))
})

# Number of relapses in 24 months prior to treatment initiation
trt_subset$PRIOR_RELAPSE_24MONS <- sapply(1:nrow(trt_subset), function(i) {
  aa = MS_attack_new[MS_attack_new$PatientNum==trt_subset$PatientNum[i],]
  bb = difftime(aa$relapse_date,trt_subset$start_date[i],units = "days")/30 # in months
  CC = abs(bb) <= 24 & bb <=0
  # Obtain relapse within designated time period
  valid_date = aa[CC,]
  # Return unique relapse dates
  return(length(with(valid_date, unique(relapse_date[clinical >= 1 | radiographic >= 1]))))
})

# Response ----------------------------------------------------------------
# Censoring time
trt_subset$LAST_DATE = trt_subset$start_date
for(i in 1:nrow(trt_subset))
{
  pid = trt_subset$PatientID[i]
  tmp <- ICD_CPT_CUI_Comb %>% 
    filter(PatientID == pid) 
  trt_subset$LAST_DATE[i] = max(tmp$StartDate)
}
# MS_cohort = MS_EHR[!is.na(MS_EHR$PatientID), ]
# trt_subset$LAST_DATE <- as.Date(sapply(trt_subset$PatientID, function(pid) {
#   out = MS_EHR$Last_date[MS_EHR$PatientID==pid]
#   if(length(out) == 0)
#     out=NA
#   return(out)
#   }))
any(trt_subset$start_date > trt_subset$LAST_DATE)

trt_subset = trt_subset[trt_subset$start_date <= trt_subset$LAST_DATE,]

# Event information
trt_subset$RELAPSE = 0
trt_subset$RELAPSE_DATE = trt_subset$LAST_DATE
for(i in 1:nrow(trt_subset))
{
  pid = trt_subset$PatientID[i]
  tmp <- MS_attack_new %>% 
    filter(PatientID == pid) %>% 
    filter(relapse_date > trt_subset$start_date[i] &
             relapse_date <= trt_subset$LAST_DATE[i])%>% 
    filter(clinical > 0 | radiographic > 0)
  if (nrow(tmp) != 0) {
    trt_subset$RELAPSE[i] = 1
    trt_subset$RELAPSE_DATE[i] = min(tmp$relapse_date)
  } 
}

# Save --------------------------------------------------------------------
## Move all non-modeling variables to first few columns of dataframe
trt_subset <- trt_subset %>% select(PatientID, PatientNum, medication_desc, start_date, 
                                    BIRTHDAY, FIRSTICD_DATE, FIRSTMSICD_DATE, 
                                    LAST_DATE, RELAPSE_DATE, RELAPSE,
                                    everything())

saveRDS(trt_subset, file=paste('MS CLIME analysis/data-processed/low-d-data-',data.name,'.rds',sep=''))

# ICD, CPT, and CUI codes (total and 3 months) prior to treatment initiation
trt_subset_highD <- trt_subset

codes_list <- unique(ICD_CPT_CUI_Comb$ConceptCd) %>% as.character()

## Remove two-digit phecode if one-digit exists
tmp = codes_list[grep("PheCode.*_[0-9][0-9]",codes_list)]
tmp = tmp[substr(tmp,1,nchar(tmp)-1)%in%codes_list]
codes_list = codes_list[!codes_list%in%tmp]
## remove zero-digit phecode if one-digit exists
tmp = codes_list[grep("PheCode.*_[0-9]",codes_list)]
tmp = tmp[substr(tmp,1,nchar(tmp)-1)%in%codes_list]
tmp = substr(tmp,1,nchar(tmp)-1)
codes_list = codes_list[!codes_list%in%tmp]

# Count codes
tmp <- sapply(1:nrow(trt_subset_highD), function(i) {
  if (i %% 20 == 0) {print(i)}
  pnum <- trt_subset$PatientNum[i]
  tx_start_date <- trt_subset$start_date[trt_subset$PatientNum == pnum]
  tmp <- ICD_CPT_CUI_Comb %>% filter(PatientNum == pnum & StartDate < tx_start_date)
  # Overall
  tmp2 <- sapply(codes_list, function(code) {
    return(tmp %>% filter(ConceptCd == code) %>% nrow)
  })
  # 3 months 
  tmp3 <- sapply(codes_list, function(code) {
    return(tmp %>% filter(StartDate >= tx_start_date - months(3) & ConceptCd == code) %>% nrow)
  })
  names(tmp2) = names(tmp3) = NULL
  return(c(tmp2, tmp3))
})
tmp <- t(tmp)
colnames(tmp) <- paste0(rep(codes_list, 2), rep(c("_OVERALL", "_3MONS"), each = length(codes_list)))

trt_subset_comb <- cbind(trt_subset_highD, tmp)
trt_subset_highD <- cbind(trt_subset_highD[,c(1:10)], tmp)
saveRDS(trt_subset_highD, 
        file=paste('MS CLIME analysis/data-processed/high-d-data-',data.name,'.rds',sep=''))
saveRDS(trt_subset_comb, 
        file=paste('MS CLIME analysis/data-processed/comb-data-',data.name,'.rds',sep=''))