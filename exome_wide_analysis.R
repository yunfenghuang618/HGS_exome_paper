################################
## Exome-wide burden analysis ##
## Author: Yunfeng Huang      ##
## Date:   Sep-20-2021        ##
################################

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

# PTV-burden by PLI>=0.9 & PLI<0.9 and all genes

ptv_tol <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/081721_ptvsum/step2/ukbb_456k_pLI_tol_maf001_ptvsum.txt")
ptv_tol_dat <- merge(dat,ptv_tol,by="eid_sample")

ptv_tol_dat <- merge(ptv_tol_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(ptv_tol_dat)[ncol(ptv_tol_dat)] <- "loco"
ptv_tol_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = ptv_tol_dat))
ptv_tol_dat$gene_res <- residuals(lm(formula = paste0("ptvsum~",paste(colnames(cov)[3:31],collapse = "+")),data = ptv_tol_dat))
ptv_tol_dat$max_grip_res <- ptv_tol_dat$max_grip_res - ptv_tol_dat$loco
ptv_res_tol <- lm(max_grip_res~gene_res,data = ptv_tol_dat)


ptv_intol <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/081721_ptvsum/step2/ukbb_456k_pLI_intol_maf001_ptvsum.txt")
ptv_intol_dat <- merge(dat,ptv_intol,by="eid_sample")

ptv_intol_dat <- merge(ptv_intol_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(ptv_intol_dat)[ncol(ptv_intol_dat)] <- "loco"
ptv_intol_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = ptv_intol_dat))
ptv_intol_dat$gene_res <- residuals(lm(formula = paste0("ptvsum~",paste(colnames(cov)[3:31],collapse = "+")),data = ptv_intol_dat))
ptv_intol_dat$max_grip_res <- ptv_intol_dat$max_grip_res - ptv_intol_dat$loco
ptv_res_intol <- lm(max_grip_res~gene_res,data = ptv_intol_dat)


ptv_all <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/081721_ptvsum/step2/ukbb_456k_gw_maf001_ptvsum.txt")
ptv_all_dat <- merge(dat,ptv_all,by="eid_sample")

ptv_all_dat <- merge(ptv_all_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(ptv_all_dat)[ncol(ptv_all_dat)] <- "loco"
ptv_all_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = ptv_all_dat))
ptv_all_dat$gene_res <- residuals(lm(formula = paste0("ptvsum~",paste(colnames(cov)[3:31],collapse = "+")),data = ptv_all_dat))
ptv_all_dat$max_grip_res <- ptv_all_dat$max_grip_res - ptv_all_dat$loco
ptv_res_all <- lm(max_grip_res~gene_res,data = ptv_all_dat)


# missense-burden by PLI>0.9 & PLI<=0.9 and all genes

ms_tol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_pLI_tol_maf001_CADD2_0_missensesum.txt")
ms_tol_dat <- merge(dat,dplyr::select(ms_tol,eid_sample,CADD2_0_sum=burdensum),by="eid_sample")
ms_tol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_pLI_tol_maf001_CADD2_20_missensesum.txt")
ms_tol_dat <- merge(ms_tol_dat,dplyr::select(ms_tol,eid_sample,CADD2_20_sum=burdensum),by="eid_sample")
ms_tol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_pLI_tol_maf001_CADD2_30_missensesum.txt")
ms_tol_dat <- merge(ms_tol_dat,dplyr::select(ms_tol,eid_sample,CADD2_30_sum=burdensum),by="eid_sample")

ms_tol_dat <- merge(ms_tol_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(ms_tol_dat)[ncol(ms_tol_dat)] <- "loco"

ms_tol_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_tol_dat))
ms_tol_dat$CADD2_0_res <- residuals(lm(formula = paste0("CADD2_0_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_tol_dat))
ms_tol_dat$CADD2_20_res <- residuals(lm(formula = paste0("CADD2_20_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_tol_dat))
ms_tol_dat$CADD2_30_res <- residuals(lm(formula = paste0("CADD2_30_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_tol_dat))
ms_tol_dat$max_grip_res <- ms_tol_dat$max_grip_res - ms_tol_dat$loco
ms_res_tol_CADD2_0 <- lm(max_grip_res~CADD2_0_res,data = ms_tol_dat)
ms_res_tol_CADD2_20 <- lm(max_grip_res~CADD2_20_res,data = ms_tol_dat)
ms_res_tol_CADD2_30 <- lm(max_grip_res~CADD2_30_res,data = ms_tol_dat)


ms_intol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_pLI_intol_maf001_CADD2_0_missensesum.txt")
ms_intol_dat <- merge(dat,dplyr::select(ms_intol,eid_sample,CADD2_0_sum=burdensum),by="eid_sample")
ms_intol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_pLI_intol_maf001_CADD2_20_missensesum.txt")
ms_intol_dat <- merge(ms_intol_dat,dplyr::select(ms_intol,eid_sample,CADD2_20_sum=burdensum),by="eid_sample")
ms_intol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_pLI_intol_maf001_CADD2_30_missensesum.txt")
ms_intol_dat <- merge(ms_intol_dat,dplyr::select(ms_intol,eid_sample,CADD2_30_sum=burdensum),by="eid_sample")

ms_intol_dat <- merge(ms_intol_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(ms_intol_dat)[ncol(ms_intol_dat)] <- "loco"

ms_intol_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_intol_dat))
ms_intol_dat$CADD2_0_res <- residuals(lm(formula = paste0("CADD2_0_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_intol_dat))
ms_intol_dat$CADD2_20_res <- residuals(lm(formula = paste0("CADD2_20_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_intol_dat))
ms_intol_dat$CADD2_30_res <- residuals(lm(formula = paste0("CADD2_30_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_intol_dat))
ms_intol_dat$max_grip_res <- ms_intol_dat$max_grip_res - ms_intol_dat$loco
ms_res_intol_CADD2_0 <- lm(max_grip_res~CADD2_0_res,data = ms_intol_dat)
ms_res_intol_CADD2_20 <- lm(max_grip_res~CADD2_20_res,data = ms_intol_dat)
ms_res_intol_CADD2_30 <- lm(max_grip_res~CADD2_30_res,data = ms_intol_dat)


ms_all <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_gw_maf001_CADD2_0_missensesum.txt")
ms_all_dat <- merge(dat,dplyr::select(ms_all,eid_sample,CADD2_0_sum=burdensum),by="eid_sample")
ms_all <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_gw_maf001_CADD2_20_missensesum.txt")
ms_all_dat <- merge(ms_all_dat,dplyr::select(ms_all,eid_sample,CADD2_20_sum=burdensum),by="eid_sample")
ms_all <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missensesum_456k_v2/ukbb_456k_gw_maf001_CADD2_30_missensesum.txt")
ms_all_dat <- merge(ms_all_dat,dplyr::select(ms_all,eid_sample,CADD2_30_sum=burdensum),by="eid_sample")

ms_all_dat <- merge(ms_all_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(ms_all_dat)[ncol(ms_all_dat)] <- "loco"

ms_all_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_all_dat))
ms_all_dat$CADD2_0_res <- residuals(lm(formula = paste0("CADD2_0_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_all_dat))
ms_all_dat$CADD2_20_res <- residuals(lm(formula = paste0("CADD2_20_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_all_dat))
ms_all_dat$CADD2_30_res <- residuals(lm(formula = paste0("CADD2_30_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = ms_all_dat))
ms_all_dat$max_grip_res <- ms_all_dat$max_grip_res - ms_all_dat$loco
ms_res_all_CADD2_0 <- lm(max_grip_res~CADD2_0_res,data = ms_all_dat)
ms_res_all_CADD2_20 <- lm(max_grip_res~CADD2_20_res,data = ms_all_dat)
ms_res_all_CADD2_30 <- lm(max_grip_res~CADD2_30_res,data = ms_all_dat)


# synonymous-burden by PLI>=0.9 & PLI<0.9

syn_tol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/synonymoussum_456k_v2/ukbb_456k_pLI_tol_maf001_synonymoussum.txt")
syn_tol_dat <- merge(dat,dplyr::select(syn_tol,eid_sample,syn_sum=burdensum),by="eid_sample")

syn_tol_dat <- merge(syn_tol_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(syn_tol_dat)[ncol(syn_tol_dat)] <- "loco"

syn_tol_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = syn_tol_dat))
syn_tol_dat$syn_sum_res <- residuals(lm(formula = paste0("syn_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = syn_tol_dat))
syn_tol_dat$max_grip_res <- syn_tol_dat$max_grip_res - syn_tol_dat$loco
syn_res_tol <- lm(max_grip_res~syn_sum_res,data = syn_tol_dat)

syn_intol <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/synonymoussum_456k_v2/ukbb_456k_pLI_intol_maf001_synonymoussum.txt")
syn_intol_dat <- merge(dat,dplyr::select(syn_intol,eid_sample,syn_sum=burdensum),by="eid_sample")

syn_intol_dat <- merge(syn_intol_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(syn_intol_dat)[ncol(syn_intol_dat)] <- "loco"

syn_intol_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = syn_intol_dat))
syn_intol_dat$syn_sum_res <- residuals(lm(formula = paste0("syn_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = syn_intol_dat))
syn_intol_dat$max_grip_res <- syn_intol_dat$max_grip_res - syn_intol_dat$loco
syn_res_intol <- lm(max_grip_res~syn_sum_res,data = syn_intol_dat)


syn_all <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/synonymoussum_456k_v2/ukbb_456k_gw_maf001_synonymoussum.txt")
syn_all_dat <- merge(dat,dplyr::select(syn_all,eid_sample,syn_sum=burdensum),by="eid_sample")

syn_all_dat <- merge(syn_all_dat,loco_t[,c("FID","chr23")],by="FID")
colnames(syn_all_dat)[ncol(syn_all_dat)] <- "loco"

syn_all_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = syn_all_dat))
syn_all_dat$syn_sum_res <- residuals(lm(formula = paste0("syn_sum~",paste(colnames(cov)[3:31],collapse = "+")),data = syn_all_dat))
syn_all_dat$max_grip_res <- syn_all_dat$max_grip_res - syn_all_dat$loco
syn_res_all <- lm(max_grip_res~syn_sum_res,data = syn_all_dat)

# count # of genes & # of variants

ptv_ngenes <- fread("~/muscle_targets/regenie/ukbb_exome_scripts/081721_ptvsum/step2/ukbb_456k_n.genes_maf001.txt")

pli <- fread("/camhpc/home/rtian/sleep/useful_files/pLI_cleaned_gnomad.v2.1.1.txt")
tol <- pli$gene[pli$pLI<0.9]
intol <- pli$gene[pli$pLI>=0.9]

CADD2_0_ngenes <- 0
CADD2_0_ngenes_tol <- 0
CADD2_0_ngenes_intol <- 0
CADD2_20_ngenes <- 0
CADD2_20_ngenes_tol <- 0
CADD2_20_ngenes_intol <- 0
CADD2_30_ngenes <- 0
CADD2_30_ngenes_tol <- 0
CADD2_30_ngenes_intol <- 0

for(i in 1:22){
  CADD2_0_genes <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missense_456k_matrix_v2/ukbb_456k_v2_chr",i,".maf001.CADD2_0.info.tsv"))
  CADD2_0_ngenes <- CADD2_0_ngenes + length(CADD2_0_genes$gene[CADD2_0_genes$n_variants>0])
  CADD2_0_ngenes_tol <- CADD2_0_ngenes_tol + length(CADD2_0_genes$gene[CADD2_0_genes$n_variants>0 & CADD2_0_genes$gene %in% tol])
  CADD2_0_ngenes_intol <- CADD2_0_ngenes_intol + length(CADD2_0_genes$gene[CADD2_0_genes$n_variants>0 & CADD2_0_genes$gene %in% intol])
  
  CADD2_20_genes <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missense_456k_matrix_v2/ukbb_456k_v2_chr",i,".maf001.CADD2_20.info.tsv"))
  CADD2_20_ngenes <- CADD2_20_ngenes + length(CADD2_20_genes$gene[CADD2_20_genes$n_variants>0])
  CADD2_20_ngenes_tol <- CADD2_20_ngenes_tol + length(CADD2_20_genes$gene[CADD2_20_genes$n_variants>0 & CADD2_20_genes$gene %in% tol])
  CADD2_20_ngenes_intol <- CADD2_20_ngenes_intol + length(CADD2_20_genes$gene[CADD2_20_genes$n_variants>0 & CADD2_20_genes$gene %in% intol])
  
  CADD2_30_genes <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missense_456k_matrix_v2/ukbb_456k_v2_chr",i,".maf001.CADD2_30.info.tsv"))
  CADD2_30_ngenes <- CADD2_30_ngenes + length(CADD2_30_genes$gene[CADD2_30_genes$n_variants>0])
  CADD2_30_ngenes_tol <- CADD2_30_ngenes_tol + length(CADD2_30_genes$gene[CADD2_30_genes$n_variants>0 & CADD2_30_genes$gene %in% tol])
  CADD2_30_ngenes_intol <- CADD2_30_ngenes_intol + length(CADD2_30_genes$gene[CADD2_30_genes$n_variants>0 & CADD2_30_genes$gene %in% intol])
  
  print(i)
} 

syn_ngenes <- 0
syn_ngenes_tol <- 0
syn_ngenes_intol <- 0

for(i in 1:22){
  syn_genes <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/synonymous_456k_matrix_v3/ukbb_456k_v3_synonymous_chr",i,".maf001.info.tsv"))
  syn_ngenes <- syn_ngenes + length(syn_genes$gene[syn_genes$n_variants>0])
  syn_ngenes_tol <- syn_ngenes_tol + length(syn_genes$gene[syn_genes$n_variants>0 & syn_genes$gene %in% tol])
  syn_ngenes_intol <- syn_ngenes_intol + length(syn_genes$gene[syn_genes$n_variants>0 & syn_genes$gene %in% intol])
}


ptv_nvars <- 0
ptv_nvars_tol <- 0
ptv_nvars_intol <- 0

for(i in 1:22){
  ptv_vars <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/ptv_456k_v4/ukbb_456k_v4_chr",i,".maf001.ptvburden.info.tsv"))
  ptv_nvars <- ptv_nvars + sum(ptv_vars$n_variants[ptv_vars$n_variants>0])
  ptv_nvars_tol <- ptv_nvars_tol + sum(ptv_vars$n_variants[ptv_vars$n_variants>0 & ptv_vars$gene %in% tol])
  ptv_nvars_intol <- ptv_nvars_intol + sum(ptv_vars$n_variants[ptv_vars$n_variants>0 & ptv_vars$gene %in% intol])
}

CADD2_0_nvars <- 0
CADD2_0_nvars_tol <- 0
CADD2_0_nvars_intol <- 0
CADD2_20_nvars <- 0
CADD2_20_nvars_tol <- 0
CADD2_20_nvars_intol <- 0
CADD2_30_nvars <- 0
CADD2_30_nvars_tol <- 0
CADD2_30_nvars_intol <- 0

for(i in 1:22){
  CADD2_0_vars <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missense_456k_matrix_v2/ukbb_456k_v2_chr",i,".maf001.CADD2_0.info.tsv"))
  CADD2_0_nvars <- CADD2_0_nvars + sum(CADD2_0_vars$n_variants[CADD2_0_vars$n_variants>0])
  CADD2_0_nvars_tol <- CADD2_0_nvars_tol + sum(CADD2_0_vars$n_variants[CADD2_0_vars$n_variants>0 & CADD2_0_vars$gene %in% tol])
  CADD2_0_nvars_intol <- CADD2_0_nvars_intol + sum(CADD2_0_vars$n_variants[CADD2_0_vars$n_variants>0 & CADD2_0_vars$gene %in% intol])
  
  CADD2_20_vars <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missense_456k_matrix_v2/ukbb_456k_v2_chr",i,".maf001.CADD2_20.info.tsv"))
  CADD2_20_nvars <- CADD2_20_nvars + sum(CADD2_20_vars$n_variants[CADD2_20_vars$n_variants>0])
  CADD2_20_nvars_tol <- CADD2_20_nvars_tol + sum(CADD2_20_vars$n_variants[CADD2_20_vars$n_variants>0 & CADD2_20_vars$gene %in% tol])
  CADD2_20_nvars_intol <- CADD2_20_nvars_intol + sum(CADD2_20_vars$n_variants[CADD2_20_vars$n_variants>0 & CADD2_20_vars$gene %in% intol])
  
  CADD2_30_vars <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/missense_456k_matrix_v2/ukbb_456k_v2_chr",i,".maf001.CADD2_30.info.tsv"))
  CADD2_30_nvars <- CADD2_30_nvars + sum(CADD2_30_vars$n_variants[CADD2_30_vars$n_variants>0])
  CADD2_30_nvars_tol <- CADD2_30_nvars_tol + sum(CADD2_30_vars$n_variants[CADD2_30_vars$n_variants>0 & CADD2_30_vars$gene %in% tol])
  CADD2_30_nvars_intol <- CADD2_30_nvars_intol + sum(CADD2_30_vars$n_variants[CADD2_30_vars$n_variants>0 & CADD2_30_vars$gene %in% intol])
  
  print(i)
} 

syn_nvars <- 0
syn_nvars_tol <- 0
syn_nvars_intol <- 0

for(i in 1:22){
  syn_vars <- fread(paste0("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/synonymous_456k_matrix_v3/ukbb_456k_v3_synonymous_chr",i,".maf001.info.tsv"))
  syn_nvars <- syn_nvars + sum(syn_vars$n_variants[syn_vars$n_variants>1])
  syn_nvars_tol <- syn_nvars_tol + sum(syn_vars$n_variants[syn_vars$n_variants>1 & syn_vars$gene %in% tol])
  syn_nvars_intol <- syn_nvars_intol + sum(syn_vars$n_variants[syn_vars$n_variants>1 & syn_vars$gene %in% intol])
}

## plot forestplot - figure 1A

data1 <- data.frame("Gene"=rep(c("pLI>=0.9","pLI<0.9","all"),each=5),"variant"=rep(c("PTV","CADD2_0","CADD2_20","CADD2_30","syn"),3),"Beta"=rep(NA,15),"SE"=rep(NA,15),"P"=rep(NA,15))

data1[1,3:5] <- summary(ptv_res_intol)$coefficients[2,c(1,2,4)]
data1[2,3:5] <- summary(ms_res_intol_CADD2_0)$coefficients[2,c(1,2,4)]
data1[3,3:5] <- summary(ms_res_intol_CADD2_20)$coefficients[2,c(1,2,4)]
data1[4,3:5] <- summary(ms_res_intol_CADD2_30)$coefficients[2,c(1,2,4)]
data1[5,3:5] <- summary(syn_res_intol)$coefficients[2,c(1,2,4)]

data1[6,3:5] <- summary(ptv_res_tol)$coefficients[2,c(1,2,4)]
data1[7,3:5] <- summary(ms_res_tol_CADD2_0)$coefficients[2,c(1,2,4)]
data1[8,3:5] <- summary(ms_res_tol_CADD2_20)$coefficients[2,c(1,2,4)]
data1[9,3:5] <- summary(ms_res_tol_CADD2_30)$coefficients[2,c(1,2,4)]
data1[10,3:5] <- summary(syn_res_tol)$coefficients[2,c(1,2,4)]

data1[11,3:5] <- summary(ptv_res_all)$coefficients[2,c(1,2,4)]
data1[12,3:5] <- summary(ms_res_all_CADD2_0)$coefficients[2,c(1,2,4)]
data1[13,3:5] <- summary(ms_res_all_CADD2_20)$coefficients[2,c(1,2,4)]
data1[14,3:5] <- summary(ms_res_all_CADD2_30)$coefficients[2,c(1,2,4)]
data1[15,3:5] <- summary(syn_res_all)$coefficients[2,c(1,2,4)]

data1$P.bonf <- p.adjust(data1$P,"bonferroni")

data1$lower <- data1$Beta-1.96*data1$SE
data1$upper <- data1$Beta+1.96*data1$SE
data1$Beta_round <- round(data1$Beta,2)

fg1_dat <- cbind(c(NA,NA,data1$Beta[11],data1$Beta[1],data1$Beta[6],NA,data1$Beta[14],data1$Beta[4],data1$Beta[9],
                   NA,data1$Beta[13],data1$Beta[3],data1$Beta[8],NA,data1$Beta[12],data1$Beta[2],data1$Beta[7],
                   NA,data1$Beta[15],data1$Beta[5],data1$Beta[10]),
                 c(NA,NA,data1$lower[11],data1$lower[1],data1$lower[6],NA,data1$lower[14],data1$lower[4],data1$lower[9],
                   NA,data1$lower[13],data1$lower[3],data1$lower[8],NA,data1$lower[12],data1$lower[2],data1$lower[7],
                   NA,data1$lower[15],data1$lower[5],data1$lower[10]),
                 c(NA,NA,data1$upper[11],data1$upper[1],data1$upper[6],NA,data1$upper[14],data1$upper[4],data1$upper[9],
                   NA,data1$upper[13],data1$upper[3],data1$upper[8],NA,data1$upper[12],data1$upper[2],data1$upper[7],
                   NA,data1$upper[15],data1$upper[5],data1$upper[10]))

tabletext <- cbind(c("Test","PTV-burden","All genes","pLI >= 0.9","pLI < 0.9","Missense (CADD>30)","All genes","pLI >= 0.9","pLI < 0.9",
                     "Missense (CADD 20-30)","All genes","pLI >= 0.9","pLI < 0.9","Missense (CADD 0-20)","All genes","pLI >= 0.9","pLI < 0.9",
                     "Synonymous","All genes","pLI >= 0.9","pLI < 0.9"),
                   c("No. of genes",NA,as.character(ptv_ngenes$n),as.character(ptv_ngenes$n2),as.character(ptv_ngenes$n1),
                     NA,as.character(CADD2_30_ngenes),as.character(CADD2_30_ngenes_intol),as.character(CADD2_30_ngenes_tol),
                     NA,as.character(CADD2_20_ngenes),as.character(CADD2_20_ngenes_intol),as.character(CADD2_20_ngenes_tol),
                     NA,as.character(CADD2_0_ngenes),as.character(CADD2_0_ngenes_intol),as.character(CADD2_0_ngenes_tol),
                     NA,as.character(syn_ngenes),as.character(syn_ngenes_intol),as.character(syn_ngenes_tol)),
                   c("P.adj",NA,as.character(format(data1$P.bonf[c(11,1,6)],scientific = T,digits = 2)),
                     NA,as.character(format(data1$P.bonf[c(14,4,9)],scientific = T,digits = 2)),
                     NA,as.character(format(data1$P.bonf[c(13,3,8)],scientific = T,digits = 2)),
                     NA,as.character(format(data1$P.bonf[c(12,2,7)],scientific = T,digits = 2)),
                     NA,as.character(format(data1$P.bonf[c(15,5,10)],scientific = T,digits = 2))))

fg1a <- grid.grabExpr(print(fg1_dat %>% forestplot(labeltext = tabletext, is.summary = c(rep(TRUE, 2), rep(FALSE, 3), TRUE, rep(FALSE, 3), TRUE, 
                                                                                         rep(FALSE, 3),TRUE, rep(FALSE, 3), TRUE, rep(FALSE, 3)),clip = c(-0.3, 0.05), xlog = FALSE, 
                                                   col = fpColors(box = c("royalblue","red"),line = "darkblue",summary = "royalblue"),
                                                   txt_gp = fpTxtGp(ticks = gpar(cex = 0.75),xlab  = gpar(cex = 0.83),label = gpar(cex=0.75)), xlab = "Beta for grip strength (kg)",
                                                   xticks = c(-0.3,-0.2,-0.1,0,0.05),zero = 0,
                                                   hrzl_lines = list("2" = gpar(lwd = 1.2)),mar = unit(rep(10, times = 4), "mm"),graph.pos = 3)))

## Mahattan and QQ plots - figure 1b,c

# gene location file
gene_pos <- read.table("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/UKB_Freeze_456k_siteonly/UKB_Freeze_456k_gene_pos.txt",header=TRUE,stringsAsFactors=FALSE)
# recode gene names in ptv gene location file
gene_pos$GENE <- gsub(",", "_", gene_pos$GENE)
gene_pos$GENE <- gsub("\\.", "_", gene_pos$GENE)
gene_pos$GENE <- gsub("-", "_", gene_pos$GENE)
gene_pos$CHR<-gsub("X", "23", gene_pos$CHR)

plt<-function(data,type){
  #compute cumulative BP
  data$CHR<-as.numeric(data$CHR)
  data <- data[with(data, order(CHR, BP)),]
  nCHR <- length(unique(data$CHR))
  data$BPcum <- NA
  s <- 0
  nbp <- c()
  for (i in unique(data$CHR)){
    nbp[i] <- max(data[data$CHR == i,]$BP)
    data[data$CHR == i,"BPcum"] <- data[data$CHR == i,"BP"] + s
    s <- s + nbp[i]
  }
  
  
  axis.set <- data %>% 
    group_by(CHR) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(data$P)))) + 2 
  sig <- 0.05/dim(data)[1]
  
  
  p1<-ggplot(data, aes(x=BPcum, y=-log10(P))) +
    # Add lines
    geom_hline(yintercept=-1*log10(0.05/dim(data)[1]), linetype="dashed", color = "#d73027") + 
    # Show all points
    geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
    scale_color_manual(values = rep(c("#357EBDB2", "#46B8DAB2"), 22 )) +
    # Annotate point Bon sig
    geom_point(data=data[data$P<sig,], color="#d73027",size=2) +
    #geom_text_repel(seed=3, data=data[data$P<sig,], aes(label=GENE), size=2, point.padding=0.1, segment.color=NA, box.padding = 0.1) +
    geom_label_repel(seed=3, data=data[data$P<sig,], aes(label=GENE), size=8/.pt, point.padding=0.1, segment.color=NA, box.padding = 0.1) +
    # Annotate point FDR<0.05
    geom_point(data=data[data$P>sig&data$FDR<0.05,], color="#fdae61",size=1.5) +
    #geom_text_repel(seed=3, data=data[data$P>sig&data$FDR<0.05,], aes(label=GENE), size=1.5, point.padding=0.1,nudge_y = 0, segment.color=NA, box.padding = 0.1) +
    geom_label_repel(seed=3, data=data[data$P>sig&data$FDR<0.05,], aes(label=GENE), size=6/.pt, point.padding=0.1,nudge_y = 0, segment.color=NA, box.padding = 0.1) +
    # Custom X axis:
    scale_x_continuous(name="Chromosome", label = c(seq(1,22,2),"X"), breaks = axis.set$center[seq(1,23,2)]) +
    scale_y_continuous(name=expression(-log[10](italic(P))), limits=c(0,max(c(-log10(data$P))*1.03,-1.03*log10(0.05/dim(data)[1]))), expand = c(0, 0) ) +  # remove space between plot area and x axis
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",
      axis.text.y=element_text(size = 9,color = "black"),
      axis.text.x=element_text(size=9,color = "black"),
      axis.title=element_text(size=10,color = "black"),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"))
  
  # QQ plot
  N_GENE<-dim(data)[1]
  data <- data[order(data$P),]
  data$P_exp <- ppoints(length(data$P))
  data$clower <- -log10(qbeta(p = (1 - 0.95) / 2, shape1 = 1:N_GENE, shape2 = N_GENE:1))
  data$cupper <- -log10(qbeta(p = (1 + 0.95) / 2, shape1 = 1:N_GENE, shape2 = N_GENE:1))
  obs_chisq <- format(round(median(qchisq(data$P,1,lower.tail=FALSE)),3), trim=TRUE, digis=3, nsmall=3)
  lambdagc <- format(round(median(qchisq(data$P,1,lower.tail=FALSE))/0.455,3), trim=TRUE, digis=3, nsmall=3)
  obs_chisq
  # [1] "0.476"
  lambdagc
  # [1] "1.046"
  p1_qq <- ggplot(data, aes(x = -1*log10(P_exp), y=-1*log10(P))) + 
    # Add points
    geom_point(size=1) +  
    # Annotate point Bon sig
    geom_point(data=data[data$P<sig,], color="#d73027",size=2) +
    #geom_text_repel(seed=3, data=data[data$P<sig,], aes(label=GENE), size=2, point.padding=0.1, segment.color=NA, box.padding = 0.1) +
    # Annotate point FDR<0.05
    geom_point(data=data[data$P>sig&data$FDR<0.05,], color="#fdae61",size=1.5) +
    #geom_text_repel(seed=3, data=data[data$P>sig&data$FDR<0.05,], aes(label=GENE), size=1.5, point.padding=0.1,nudge_y = 0, segment.color=NA, box.padding = 0.1) +
    # Add diagonal line
    geom_abline(intercept=0, slope=1, colour="gray20") + 
    # Add 95% CI
    geom_ribbon(mapping = aes(x = -1*log10(P_exp), ymin = clower, ymax = cupper), alpha = 0.1) +
    # Add significant line
    geom_abline(intercept= -1*log10(0.05/dim(data)[1]), slope=0, linetype = "dashed", colour="#d73027") + 
    # Custom the theme:
    theme_bw() +
    theme( 
      legend.position="none",  
      axis.text.y=element_text(size = 9,color = "black"),
      axis.text.x=element_text(size=9,color = "black"),
      axis.ticks = element_line(size = 0.5),
      axis.title=element_text(size=10,color = "black"),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.grid = element_blank()
    ) +
    # Custom axis:
    scale_x_continuous(name=expression(Expected ~ ~-log[10](italic(P))), limits=c(0,max(-log10(data$P_exp))*1.03), expand = c(0, 0) ) +
    scale_y_continuous(name=expression(Observed ~ ~-log[10](italic(P))), limits=c(0,max(c(-log10(data$P))*1.03,-1.03*log10(0.05/dim(data)[1]))), expand = c(0, 0) ) 
  
  return(list(p1,p1_qq))
}

setwd("/camhpc/home/yhuang5/muscle_targets/regenie/")

data_ptv <- get(load("/camhpc/home/yhuang5/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.rda"))
ptv_carriers <- get(load("/camhpc/home/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/0823_ptv_carrier.rda"))
data_ptv <- merge(data_ptv,ptv_carriers,by="gene") %>% filter(N_carriers>=10)
ptv_final <- merge(gene_pos,data_ptv,by.x="GENE",by.y="gene")
ptv_final$Bonf <- p.adjust(ptv_final$P,"bonferroni")
ptv_final$FDR <- p.adjust(ptv_final$P,"fdr")

man_qq <- plt(data = ptv_final,type = "PTV")
#pdf(file = "/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.exomewide_Man_QQ.final.pdf",width = 7,height =8.2,family = "ArialMT")
fg1_all <- ggarrange(fg1a,ggarrange(man_qq[[1]],man_qq[[2]],labels = c("b","c"),widths = c(2,1),ncol = 2,nrow = 1),nrow = 2,labels = "a",heights = c(1,1.1))
#print(fg1_all)
#dev.off()
fg1_all
ggsave(filename = "UKBB_max_grip.regenie.PTV.exomewide_Man_QQ.final.pdf",
       device = "pdf",path = "/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/",width = 180,height = 210,units = "mm",dpi = 300)

