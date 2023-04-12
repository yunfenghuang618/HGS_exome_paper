library(data.table)

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

# pathway PTV-burden

args <- commandArgs(trailingOnly = T)
pathway <- args[1]

ptvsum <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/pathwaysum_456k_v3/ptvsum/",pathway,".maf001_ptvsum.txt"))

results <- as.data.frame(matrix(NA,1,5))
colnames(results) <- c("pathway","N","beta","se","P")

for(i in 2:ncol(ptvsum)){
  data <- merge(dat,dplyr::select(ptvsum,eid_sample,all_of(i)),by="eid_sample")
  pathway_name <- colnames(data)[ncol(data)]
  colnames(data)[ncol(data)] <- "pathway"
  data <- merge(data,loco_t[,c("FID","chr23")],by="FID")
  colnames(data)[ncol(data)] <- "loco"
  data$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
  data$pathway_res <- residuals(lm(formula = paste0("pathway~",paste(colnames(cov)[3:31],collapse = "+")),data = data))
  data$max_grip_res <- data$max_grip_res - data$loco
  res <- lm(max_grip_res~pathway_res,data = data)
  if(nrow(summary(res)$coeff)==2){
    asso <- data.frame("pathway"=pathway_name,"N"=nrow(data),"beta"=summary(res)$coeff[2,1],"se"=summary(res)$coeff[2,2],"P"=summary(res)$coeff[2,4])
    results <- rbind(results,asso)
  }
}

results <- results[-1,]
save(results,file=paste0("~/muscle_targets/regenie/ukbb_exome_scripts/UKBB_max_grip.regenie.PTV.",pathway,".rda"))

