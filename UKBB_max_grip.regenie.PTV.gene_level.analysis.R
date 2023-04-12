library(data.table)
library(dplyr)

pheno <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_maxHGS_WBLM.pheno")
cov <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_maxHGS_WBLM.cov")
dat <- merge(pheno,cov,by=c("FID","IID"))

loco <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_maxHGS_WBLM_step1_1.loco")
loco_t <- as.data.frame(t(loco[,-1]))
colnames(loco_t) <- paste0("chr",1:23)
loco_t$FID <- row.names(loco_t)
getID <- function(x){return(strsplit(x,"_")[[1]][1])}
loco_t$FID <- sapply(loco_t$FID,getID)
loco_t$FID <- as.numeric(loco_t$FID)

bridge <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/WES_26041_450k_bridge.csv")
dat <- merge(dat,dplyr::select(bridge,FID=EID_26041,eid_sample=sample),by="FID")

# PTV-burden associations

ptv_files <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/ptv_files.list",header = F)

results <- as.data.frame(matrix(NA,1,5))
colnames(results) <- c("gene","N","beta","se","P")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  for(j in 141:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:140,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res,data = data)
    if(nrow(summary(res)$coeff)==2){
      asso <- data.frame("gene"=gene,"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4])
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.rda")


# cognitive function adjustment

cog_pheno <- fread("/camhpc/dept/human_genetics/UKBB_phenotypes/ukbb500k.cogphenoses.v3.txt")

dat <- merge(dat,select(cog_pheno,FID=f_eid,edu_baseline,rt_baseline,fiall_baseline),by="FID",all.x=T)

ptv_files <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/ptv_files.list",header = F)

results <- as.data.frame(matrix(NA,1,5))
colnames(results) <- c("gene","N","beta","se","P")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  for(j in 144:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:143,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res+edu_baseline,data = data)
    if(nrow(summary(res)$coeff)==3){
      asso <- data.frame("gene"=gene,"N"=nrow(filter(data,!is.na(edu_baseline))),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4])
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.edu_baseline_adj.rda")

results <- as.data.frame(matrix(NA,1,5))
colnames(results) <- c("gene","N","beta","se","P")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  for(j in 144:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:143,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res+rt_baseline,data = data)
    if(nrow(summary(res)$coeff)==3){
      asso <- data.frame("gene"=gene,"N"=nrow(filter(data,!is.na(rt_baseline))),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4])
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.rt_baseline_adj.rda")


# Assessment center adjustment

ses <- fread("/home/cchen11/ukbb_cogexome/data/pheno_v3/ukbb500k.cogphenoses.v1.txt",header=TRUE,stringsAsFactors=FALSE) %>% dplyr::select(f_eid,contains("assess_center"))
dat <- dat %>% left_join(dplyr::select(ses,FID=f_eid,contains("assess_center1")),by="FID")

ptv_files <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/ptv_files.list",header = F)

results <- as.data.frame(matrix(NA,1,6))
colnames(results) <- c("gene","N","beta","se","P","N_carrier")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

covars <- c(colnames(cov)[3:31],colnames(dplyr::select(dat,contains("assess"))))

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  for(j in 163:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:162,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(covars,collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(covars,collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res,data = data)
    if(nrow(summary(res)$coeff)==2){
      asso <- data.frame("gene"=gene,"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4],"N_carrier"=sum(data$gene==1,na.rm = T))
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.assess_center_adj.rda")


# Sex-specific analysis

ptv_files <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/ptv_files.list",header = F)

results_male <- as.data.frame(matrix(NA,1,6))
results_female <- as.data.frame(matrix(NA,1,6))

colnames(results_male) <- c("gene","N","beta","se","P","N_carrier")
colnames(results_female) <- c("gene","N","beta","se","P","N_carrier")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

covars <- c("age","age2","height","height2",paste0("PC",c(1:20)))

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  for(j in 141:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:140,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"

    data_male <- dplyr::filter(data,sex==1)
    data_female <- dplyr::filter(data,sex==0)

    data_male$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(covars,collapse = "+")),data = data_male))
    data_male$gene_res <- residuals(lm(formula = paste0("gene~",paste(covars,collapse = "+")),data = data_male))
    data_male$max_grip_res <- data_male$max_grip_res - data_male$loco
    res_male <- lm(max_grip_res~gene_res,data = data_male)
    if(nrow(summary(res_male)$coeff)==2){
      asso_male <- data.frame("gene"=gene,"N"=nrow(data_male),"beta"=summary(res_male)$coeff[2,1],"se"=summary(res_male)$coeff[2,2],"P"=summary(res_male)$coeff[2,4],"N_carrier"=sum(data_male$gene==1,na.rm = T))
      results_male <- rbind(results_male,asso_male)
    }
    data_female$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(covars,collapse = "+")),data = data_female))
    data_female$gene_res <- residuals(lm(formula = paste0("gene~",paste(covars,collapse = "+")),data = data_female))
    data_female$max_grip_res <- data_female$max_grip_res - data_female$loco
    res_female <- lm(max_grip_res~gene_res,data = data_female)
    if(nrow(summary(res_female)$coeff)==2){
      asso_female <- data.frame("gene"=gene,"N"=nrow(data_female),"beta"=summary(res_female)$coeff[2,1],"se"=summary(res_female)$coeff[2,2],"P"=summary(res_female)$coeff[2,4],"N_carrier"=sum(data_female$gene==1,na.rm = T))
      results_female <- rbind(results_female,asso_female)
    }
  }
  print(i)
}
results_male <- results_male[-1,]
results_female <- results_female[-1,]
save(results_male,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.male.rda")
save(results_female,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.female.rda")

# Disease sensitivity analysis

load("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/disease_sensitivity/UKBB_disease_extract_11122022.rda")
sensitivity1 <- filter(ukbb_icd,Osteoarthritis==1 | RA==1 | Osteoporosis==1 | Dupuytren==1 | Osteoarthritis_self==1 | RA_self==1 | Osteoporosis_self==1 | Dupuytren_self==1)$f_eid
sensitivity2 <- filter(ukbb_icd,Osteoarthritis==1 | RA==1 | Osteoporosis==1 | Dupuytren==1 | Osteoarthritis_self==1 | RA_self==1 | Osteoporosis_self==1 | Dupuytren_self==1 | cancer==1 | cancer_self==1)$f_eid

ptv_files <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/ptv_files.list",header = F)

results <- as.data.frame(matrix(NA,1,6))
colnames(results) <- c("gene","N","beta","se","P","Ncarriers")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  ptv_dat <- filter(ptv_dat,!FID %in% sensitivity1)
  for(j in 141:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:140,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res,data = data)
    if(nrow(summary(res)$coeff)==2){
      asso <- data.frame("gene"=gene,"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4],"Ncarriers"=sum(data$gene>0,na.rm = T))
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/disease_sensitivity/UKBB_max_grip.regenie.PTV.results.sensitivity1.rda")

results <- as.data.frame(matrix(NA,1,6))
colnames(results) <- c("gene","N","beta","se","P","Ncarriers")

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  ptv_dat <- filter(ptv_dat,!FID %in% sensitivity2)
  for(j in 141:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:140,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res,data = data)
    if(nrow(summary(res)$coeff)==2){
      asso <- data.frame("gene"=gene,"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4],"Ncarriers"=sum(data$gene>0,na.rm = T))
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/disease_sensitivity/UKBB_max_grip.regenie.PTV.results.sensitivity2.rda")


# Muscle mass adjustment

ptv_files <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/ptv_files.list",header = F)

results <- as.data.frame(matrix(NA,1,6))
colnames(results) <- c("gene","N","beta","se","P","Ncarriers")

getchr <- function(x){return(strsplit(strsplit(x,"chr")[[1]][2],"\\.")[[1]][1])}

for(i in 1:nrow(ptv_files)){
  chr <- getchr(ptv_files$V1[i])
  if(chr=="X"){next}
  ptv <- fread(paste("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/",ptv_files$V1[i],sep=""))
  ptv_dat <- merge(dat,ptv,by="eid_sample")
  for(j in 141:ncol(ptv_dat)){
    data <- dplyr::select(ptv_dat,c(1:140,all_of(j)))
    gene <- colnames(data)[ncol(data)]
    colnames(data)[ncol(data)] <- "gene"
    data <- merge(data,loco_t[,c("FID",paste0("chr",chr))],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$gene_res <- residuals(lm(formula = paste0("gene~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~gene_res+whole_body_lean,data = data)
    if(nrow(summary(res)$coeff)==3){
      asso <- data.frame("gene"=gene,"N"=nrow(filter(data,!is.na(whole_body_lean))),"beta"=summary(res)$coeff[2,1],
                         "se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4],"Ncarriers"=sum(filter(data,!is.na(whole_body_lean))$gene>0,na.rm = T))
      results <- rbind(results,asso)
    }
  }
  print(i)
}
results <- results[-1,]
save(results,file="/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/disease_sensitivity/UKBB_max_grip.regenie.PTV.results.updated.muscle_mass_adj.rda")




