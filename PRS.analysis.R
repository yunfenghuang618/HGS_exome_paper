library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(forestplot)

pheno <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_maxHGS_WBLM.pheno")
cov <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_maxHGS_WBLM.cov")
dat <- merge(pheno,cov,by=c("FID","IID"))

loco <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_maxHGS_WBLM_step1_1.loco")
loco_t <- as.data.frame(t(loco[,-1]))
colnames(loco_t) <- paste0("chr",1:23)
loco_t$FID <- row.names(loco_t)
getID <- function(x){return(strsplit(x,"_")[[1]][1])}
loco_t$FID <- sapply(loco_t$FID,getID)
loco_t$FID <- as.numeric(loco_t$FID)

bridge <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/WES_26041_450k_bridge.csv")
dat <- merge(dat,dplyr::select(bridge,FID=EID_26041,eid_sample=sample),by="FID")

# PRS analysis

# Exclude PTV carriers and same number of randomly selected non-carriers from the UKBB to derive PRS
ptv_invitae_carriers <- get(load("~/muscle_targets/regenie/ukbb_exome_scripts/invitae_PTV_carriers_Feb15_2022.rda"))
carrier <- unique(strsplit(paste(paste0(ptv_invitae_carriers$homo,collapse = ";"),paste0(ptv_invitae_carriers$hetero,collapse = ";"),sep=";"),";")[[1]])
ptv_carrier <- dat$FID[dat$eid_sample %in% carrier]
ptv_noncarrier_match <- sample(dat$FID[!dat$eid_sample %in% carrier],size = length(ptv_carrier),replace = F)

gwas_keep <- fread("~/muscle_targets/regenie/UKBB_qc_pass.id",header = F)
gwas_keep <- filter(gwas_keep, !(V1 %in% c(ptv_carrier,ptv_noncarrier_match)))
write.table(gwas_keep,file="~/muscle_targets/regenie/UKBB_qc_pass_invitae_ptv_exclude.id",row.names = F,quote = F,col.names = F,sep="\t")
gwas_exclude <- rbind(cbind(ptv_carrier,ptv_carrier),cbind(ptv_noncarrier_match,ptv_noncarrier_match))
write.table(gwas_exclude,file="~/muscle_targets/regenie/UKBB_maxHGS_WBLM_invitae_ptv_exclude_step2/PRScs/UKBB_maxHGS_invitae_ptv_exclude_02222022.id",row.names = F,quote = F,col.names = F,sep="\t")

# Interaction test between PTV-burden and PRS

ptv_invitae_carriers <- get(load("~/muscle_targets/regenie/ukbb_exome_scripts/invitae_PTV_carriers_Feb15_2022.rda"))
carrier <- unique(strsplit(paste(paste0(ptv_invitae_carriers$homo,collapse = ";"),paste0(ptv_invitae_carriers$hetero,collapse = ";"),sep=";"),";")[[1]])
ptv_invitae_carriers.homo <- strsplit(paste(ptv_invitae_carriers$homo,collapse = ";"),";")[[1]]
ptv_invitae_carriers.hetero <- strsplit(paste(ptv_invitae_carriers$hetero,collapse = ";"),";")[[1]]

invitae_genes <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/invitae_PTV_carriers.genelist.tsv")
ptv_invitae_carriers.hetero.dom <- strsplit(paste(ptv_invitae_carriers$hetero[ptv_invitae_carriers$gene %in% invitae_genes$V1[grepl("AD",invitae_genes$V3)]],collapse = ";"),";")[[1]]
ptv_invitae_carriers.hetero.rec <- strsplit(paste(ptv_invitae_carriers$hetero[ptv_invitae_carriers$gene %in% invitae_genes$V1[invitae_genes$V3=="AR"]],collapse = ";"),";")[[1]]

dat$carrier <- ifelse(dat$eid_sample %in% ptv_invitae_carriers.homo,"homo",
                      ifelse(dat$eid_sample %in% ptv_invitae_carriers.hetero.dom,"Heterozygous, autosomal dominant genes",
                             ifelse(dat$eid_sample %in% ptv_invitae_carriers.hetero.rec,"Heterozygous, autosomal recessive genes",
                                    ifelse(dat$eid_sample %in% carrier,"hetero_other","Non-carrier"))))
dat <- dplyr::filter(dat,carrier != "homo" & carrier != "hetero_other")

prs <- fread("~/muscle_targets/regenie/UKBB_maxHGS_WBLM_invitae_ptv_exclude_step2/UKBB_maxHGS_invitae_ptv_exclude.PRScs_pst_eff_a1_b0.5_phiauto_chr1.sscore")
prs <- dplyr::select(prs,eid=`#IID`,prs_sum=SCORE1_SUM)

for(i in 2:22){
  prs_temp <- fread(paste0("~/muscle_targets/regenie/UKBB_maxHGS_WBLM_invitae_ptv_exclude_step2/UKBB_maxHGS_invitae_ptv_exclude.PRScs_pst_eff_a1_b0.5_phiauto_chr",i,".sscore"))
  prs_temp <- dplyr::select(prs_temp,eid=`#IID`,prs_temp=SCORE1_SUM)
  prs <- merge(prs,prs_temp,by="eid")
  prs$prs_sum <- prs$prs_sum + prs$prs_temp
  prs <- dplyr::select(prs,eid,prs_sum)
  #print(i)
}

dat <- merge(dat,prs,by.x="IID",by.y="eid")
dat$prs_sum_std <- (dat$prs_sum - mean(dat$prs_sum)) / sd(dat$prs_sum)

unrel <- fread("~/lipids_heiko/LoF_UKB/2020Update/ukbb_derived_phenotypes_20200219.csv") %>% dplyr::select(IND,UNREL)
dat <- dplyr::filter(dat,IID %in% dplyr::filter(unrel,UNREL==1)$IND)
dat$ptv_carrier <- ifelse(dat$carrier!="Non-carrier",1,0)

dat$carrier <- factor(dat$carrier,levels = c("Non-carrier","Heterozygous, autosomal dominant genes","Heterozygous, autosomal recessive genes"))
res_prs_ptv <- lm(formula = paste0("max_grip~prs_sum_std+carrier+prs_sum_std:carrier+",paste(colnames(cov)[3:31],collapse = "+")),data=dat)
res_prs_ptv2 <- lm(formula = paste0("max_grip~prs_sum_std+carrier+",paste(colnames(cov)[3:31],collapse = "+")),data=dat)

dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data=dat))
cut_prs_100 <- quantile(dat$prs_sum_std, probs=seq(0,1,0.01))
dat$prs_sum_std_pct <- cut(dat$prs_sum_std, cut_prs_100, include.lowest = TRUE)

max_grip_res_mean_noLoF <- with(dat[which(dat$carrier=="Non-carrier"),],
                                tapply(max_grip_res, prs_sum_std_pct, mean))
max_grip_res_mean_LOFdom <- with(dat[which(dat$carrier=="Heterozygous, autosomal dominant genes"),],
                                 tapply(max_grip_res, prs_sum_std_pct, mean))
max_grip_res_mean_LOFrec <- with(dat[which(dat$carrier=="Heterozygous, autosomal recessive genes"),],
                                 tapply(max_grip_res, prs_sum_std_pct, mean))

output <- cbind.data.frame(c(seq(1,100,1),seq(1,100,1),seq(1,100,1)),
                           c(rep("Non-carrier",100),rep("Heterozygous, autosomal dominant genes",100),rep("Heterozygous, autosomal recessive genes",100)),
                           c(max_grip_res_mean_noLoF,max_grip_res_mean_LOFdom,max_grip_res_mean_LOFrec))
colnames(output) <- c("PRS_percentile", "LoF", "Mean_residualized_hand_grip_strength")
output$PRS_percentile <- as.numeric(output$PRS_percentile)
output$Mean_residualized_hand_grip_strength <- as.numeric(output$Mean_residualized_hand_grip_strength)

pdf(file="~/muscle_targets/regenie/ukbb_exome_scripts/UKBB_max_grip.PRS.PTV_invitae.v2.02242022.pdf",width = 14,height = 12,bg = "white")
p10 <- ggplot(output, aes(x=PRS_percentile, y=Mean_residualized_hand_grip_strength, colour = LoF))
p10 <- p10 + geom_point() + geom_smooth(method = "lm", formula = y ~ x, se = TRUE)
p10 <- p10 + xlab("PRS by percentile") + ylab("Mean residualized hand grip strength") + labs(color="PTV in NMD genes")
p10 <- p10 + theme_pubr() + theme(panel.border = element_blank(), axis.line = element_line(colour = "black"))
p10 <- p10 + theme(axis.text = element_text(size=16),axis.title = element_text(size=18),legend.text = element_text(size=14),legend.title = element_text(size=16))
print(p10)
dev.off()
