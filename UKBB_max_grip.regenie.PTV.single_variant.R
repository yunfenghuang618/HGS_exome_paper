library(data.table)
library(dplyr)

ptv_raw <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/UKB_Freeze_456k_HC_LOF_canonical_bychr/chr2_lof_carriers_v3.txt")) %>% dplyr::filter(gene=="TTN")
ptv_raw <- ptv_raw[which(ptv_raw$alt_af<0.001),]

pheno <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM.pheno")
cov <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM.cov")
dat <- merge(pheno,cov,by=c("FID","IID"))

loco <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM_step1_1.loco")
loco_t <- as.data.frame(t(loco[,-1]))
colnames(loco_t) <- paste0("chr",1:23)
loco_t$FID <- row.names(loco_t)
getID <- function(x){return(strsplit(x,"_")[[1]][1])}
loco_t$FID <- sapply(loco_t$FID,getID)
loco_t$FID <- as.numeric(loco_t$FID)

bridge <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/WES_26041_450k_bridge.csv")
dat <- merge(dat,dplyr::select(bridge,FID=EID_26041,eid_sample=sample),by="FID")

results <- as.data.frame(matrix(NA,1,6))
colnames(results) <- c("Variant","N","beta","se","P","Ncarriers")

getid <- function(x){return(strsplit(x,":")[[1]][1])}

for (i in 1:length(ptv_raw$gene)){
  data <- dat
  carriers <- as.character(sapply(strsplit(ptv_raw$carriers[i],",")[[1]],getid))
  data$var <- ifelse(data$eid_sample %in% carriers,1,0)
  if(!is.na(ptv_raw$missing[i])){
    miss_ids <- as.character(strsplit(ptv_raw$missing[i],",")[[1]])
    data$var <- ifelse(data$eid_sample %in% miss_ids,NA,data$var)
  }
  if(sum(data$var==1,na.rm = T)>0){
    data <- filter(data,!is.na(var))
    data <- merge(data,loco_t[,c("FID","chr2")],by="FID")
    colnames(data)[ncol(data)] <- "loco"

    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$var_res <- residuals(lm(formula = paste0("var~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~var_res,data = data)
    if(nrow(summary(res)$coeff)==2){
      asso <- data.frame("Variant"=ptv_raw$variant_id[i],"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4],"Ncarriers"=sum(data$var==1,na.rm = T))
      results <- rbind(results,asso)
    }
  }
  if(ceiling(i/100)==(i/100)){
    print(i)
  }
}
results <- results[-1,]
save(results,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.single_variant.TTN.rda")


ptv_raw <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/UKB_Freeze_456k_HC_LOF_canonical_bychr/chr1_lof_carriers_v3.txt")) %>% dplyr::filter(gene=="OBSCN")
ptv_raw <- ptv_raw[which(ptv_raw$alt_af<0.001),]

pheno <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM.pheno")
cov <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM.cov")
dat <- merge(pheno,cov,by=c("FID","IID"))

loco <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM_step1_1.loco")
loco_t <- as.data.frame(t(loco[,-1]))
colnames(loco_t) <- paste0("chr",1:23)
loco_t$FID <- row.names(loco_t)
getID <- function(x){return(strsplit(x,"_")[[1]][1])}
loco_t$FID <- sapply(loco_t$FID,getID)
loco_t$FID <- as.numeric(loco_t$FID)

bridge <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/WES_26041_450k_bridge.csv")
dat <- merge(dat,dplyr::select(bridge,FID=EID_26041,eid_sample=sample),by="FID")

results <- as.data.frame(matrix(NA,1,6))
colnames(results) <- c("Variant","N","beta","se","P","Ncarriers")

getid <- function(x){return(strsplit(x,":")[[1]][1])}

for (i in 1:length(ptv_raw$gene)){
  data <- dat
  carriers <- as.character(sapply(strsplit(ptv_raw$carriers[i],",")[[1]],getid))
  data$var <- ifelse(data$eid_sample %in% carriers,1,0)
  if(!is.na(ptv_raw$missing[i])){
    miss_ids <- as.character(strsplit(ptv_raw$missing[i],",")[[1]])
    data$var <- ifelse(data$eid_sample %in% miss_ids,NA,data$var)
  }
  if(sum(data$var==1,na.rm = T)>0){
    data <- filter(data,!is.na(var))
    data <- merge(data,loco_t[,c("FID","chr1")],by="FID")
    colnames(data)[ncol(data)] <- "loco"
    
    data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$var_res <- residuals(lm(formula = paste0("var~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
    data$max_grip_res <- data$max_grip_res - data$loco
    res <- lm(max_grip_res~var_res,data = data)
    if(nrow(summary(res)$coeff)==2){
      asso <- data.frame("Variant"=ptv_raw$variant_id[i],"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4],"Ncarriers"=sum(data$var==1,na.rm = T))
      results <- rbind(results,asso)
    }
  }
  if(ceiling(i/100)==(i/100)){
    print(i)
  }
}
results <- results[-1,]
save(results,file="~/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.single_variant.OBSCN.rda")

