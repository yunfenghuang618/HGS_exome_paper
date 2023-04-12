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

# HPA tissue-based PTV-burden

protatlas <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/ptvsum_tissue/tissue_category_rna_any_Tissue.tsv")
tissues <- strsplit(paste(protatlas$`RNA tissue specific nTPM`,collapse = ";"),";")[[1]]
gettissue <- function(x){return(strsplit(x,": ")[[1]][1])}
tissues <- unique(sapply(tissues,gettissue))

tissue_res <- as.data.frame(matrix(NA,length(tissues),4))
colnames(tissue_res) <- c("tissue","beta","se","P")

for(i in 1:length(tissues)){
  tissue_res$tissue[i] <- tissues[i]
  tissue_ptv <- get(load(paste("/camhpc/home/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/ptvsum_tissue/120321_ptvsum_",tissues[i],".rda",sep = "")))
  tissue_dat <- merge(dat,tissue_ptv,by="eid_sample")
  
  tissue_dat <- merge(tissue_dat,loco_t[,c("FID","chr23")],by="FID")
  colnames(tissue_dat)[ncol(tissue_dat)] <- "loco"
  tissue_dat$max_grip_res <- residuals(lm(formula = paste0("max_grip~",paste(colnames(cov)[3:31],collapse = "+")),data = tissue_dat))
  tissue_dat$gene_res <- residuals(lm(formula = paste0("ptvsum~",paste(colnames(cov)[3:31],collapse = "+")),data = tissue_dat))
  tissue_dat$max_grip_res <- tissue_dat$max_grip_res - tissue_dat$loco
  ptv_res_tissue <- lm(max_grip_res~gene_res,data = tissue_dat)
  tissue_res$beta[i] <- summary(ptv_res_tissue)$coeff[2,1]
  tissue_res$se[i] <- summary(ptv_res_tissue)$coeff[2,2]
  tissue_res$P[i] <- summary(ptv_res_tissue)$coeff[2,4]
  print(i)
}
tissue_res$FDR_q <- p.adjust(tissue_res$P,"fdr")

tissue_ngenes <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/ptvsum_tissue/ptvsum_elevatedgenes.ngenes",header = F)
getngenes <- function(x){return(strsplit(x," ")[[1]][1])}
tissue_ngenes$ngenes <- sapply(tissue_ngenes$V2,getngenes)
tissue_res$ngenes <- as.numeric(tissue_ngenes$ngenes)

# forestplot

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

pdf(file="~/muscle_targets/regenie/ukbb_exome_scripts/UKBB_max_grip.regenie.exomewide.v6.pdf",width = 20,height = 18,bg = "white")
fg1_dat %>% forestplot(labeltext = tabletext, is.summary = c(rep(TRUE, 2), rep(FALSE, 3), TRUE, rep(FALSE, 3), TRUE, rep(FALSE, 3),
                                                             TRUE, rep(FALSE, 3), TRUE, rep(FALSE, 3)),clip = c(-0.3, 0.05), xlog = FALSE, 
                       col = fpColors(box = c("royalblue","red"),
                                      line = "darkblue",
                                      summary = "royalblue"),
                       txt_gp = fpTxtGp(ticks = gpar(cex = 2),xlab  = gpar(cex = 2.5),label = gpar(cex=2.5)), xlab = "Beta for grip strength (Kg)",
                       xticks = c(-0.3,-0.2,-0.1,0,0.05),zero = 0,
                       hrzl_lines = list("2" = gpar(lwd = 1.2)),mar = unit(rep(10, times = 4), "mm"),
                       graph.pos = 3)
dev.off()

tissue_res <- arrange(tissue_res,P)
tissue_res$lower <- tissue_res$beta-1.96*tissue_res$se
tissue_res$upper <- tissue_res$beta+1.96*tissue_res$se
tissue_res$tissue[tissue_res$tissue=="endometrium 1"] <- "endometrium"
tissue_res$tissue[tissue_res$tissue=="stomach 1"] <- "stomach"
tissue_res$tissue[tissue_res$tissue=="skin 1"] <- "skin"

fg5_dat <- cbind(c(NA,tissue_res$beta),c(NA,tissue_res$lower),c(NA,tissue_res$upper))
tabletext <- cbind(c("Tissue",tissue_res$tissue),c("No. of genes",tissue_res$ngenes),c("FDR",format(tissue_res$FDR_q,scientific = T,digits = 2)))

pdf(file="~/muscle_targets/regenie/ukbb_exome_scripts/UKBB_max_grip.regenie.tissue_elevation.v2.pdf",width = 20,height = 18,bg = "white")
fg5_dat %>% forestplot(labeltext = tabletext, is.summary = c(TRUE, rep(FALSE, 36)),clip = c(-0.3, 0.2), xlog = FALSE, 
                       col = fpColors(box = "royalblue",
                                      line = "darkblue",
                                      summary = "royalblue"),
                       txt_gp = fpTxtGp(ticks = gpar(cex = 2),xlab  = gpar(cex = 2.5),label = gpar(cex=2.5)), xlab = "Beta for grip strength (Kg)",
                       xticks = c(-0.3,-0.2,-0.1,0,0.1,0.2),zero = 0,
                       hrzl_lines = list("2" = gpar(lwd = 1.2)),mar = unit(rep(10, times = 4), "mm"),graph.pos = 3)
dev.off()
