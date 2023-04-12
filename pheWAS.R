library(logistf)
library(data.table)
library(dplyr)

irnt <- function(x, k=3/8, randomties=FALSE, seed=NULL){
  if (randomties == FALSE){
    set.seed(seed)
    return(qnorm((rank(x,na.last="keep")-k)/(sum(!is.na(x))+1-(2*k))))
  }else if (randomties == TRUE) {
    # Modified from Pain et al. Eur J Hum Genet 2018
    set.seed(seed)
    return(qnorm((rank(x,na.last="keep",ties.method='random')-k)/(sum(!is.na(x))+1-(2*k))))
  }
}

args <- commandArgs(trailingOnly = TRUE)
gene <- as.character(args[1])
phenofile <- as.character(args[2])
ptvfile <- as.character(args[3])

print("Readling phenotype file")

all_pheno <- fread(phenofile, header = TRUE)
covar <- fread("/home/cchen11/ukbb_gwas/covar_sampleQCed.ukb37427_sep2019.tsv",header=TRUE,stringsAsFactors=FALSE)[,c("IID", "sex", "age_mean0", "age_mean0_2", "sex_age_mean0","sex_age_mean0_2")]
colnames(covar)[1] <- c("IND")
all_phenocovar <- left_join(all_pheno, covar, by="IND")

idmap <- fread("/camhpc/dept/human_genetics/UKBB_exome_variant_450k/WES_26041_450k_bridge.csv")[,c(1,2,5)]
colnames(idmap) <- c("IID", "exomeID", "IND" )
all_phenocovarid <- inner_join(idmap, all_phenocovar, by="IND")

derive_pheno <- fread("~/lipids_heiko/LoF_UKB/2020Update/ukbb_derived_phenotypes_20200219.csv")
phenos <- as.data.frame(inner_join(all_phenocovarid, filter(derive_pheno,BRITISH==1 & UNREL==1) %>% select(IND,all_of(paste0("PC",1:20))), by="IND"))

bl_phenos <- c("IID", "exomeID", "IND", "sex", "age_mean0", "age_mean0_2", "sex_age_mean0", "sex_age_mean0_2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")

pheno_list <- setdiff(names(phenos), bl_phenos)

exclude_i <- c()
for(i in pheno_list){
  if(length(levels(factor(phenos[,i])))==2){
    if(sum(phenos[,i],na.rm=T) < 100){
      exclude_i <- append(exclude_i,i)
    }
  }
}

if (is.null(exclude_i)==FALSE){
  phenos <- phenos[,-which(names(phenos) %in% exclude_i) ]  
}

pheno_list <- setdiff(names(phenos), bl_phenos)

print("Reading carrier tables")

lof_table <- fread(ptvfile) %>% select(eid_sample,all_of(gene))
phenos <- inner_join(phenos,select(lof_table,IID=eid_sample,all_of(gene)),by="IID")
colnames(phenos)[ncol(phenos)] <- "lof"

x <- phenos[,"lof"]
covar_mat <- as.matrix(phenos[,c("sex", "age_mean0", "age_mean0_2", "sex_age_mean0", "sex_age_mean0_2", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20")])

phewas_res <- data.frame()
phewas_key <- read.csv("/home/cchen11/ukbb_cogexome/data/pheno_phewas_v1/phewas_pheno_key_456.csv", header = TRUE, stringsAsFactors=FALSE)
names(phewas_key) <- c("pheno","long_pheno")

for(i in 1:length(pheno_list)){
  #if(i==10){break}
  pheno_name <- pheno_list[i]
  print(pheno_name)
  if ((pheno_name %in% phewas_key$pheno)==TRUE){
    long_pheno <- phewas_key$long_pheno[which(phewas_key$pheno == pheno_name)]  
  } else{
    long_pheno <- pheno_name
  }
  y <- phenos[,pheno_name]
  if(length(levels(factor(y)))==2){
    pheno_type <- "binary"
    n_cases <- length(which(y == 1))
    n_controls <- length(which(y == 0))
    if(n_cases < 100 | n_controls < 100){
      next()
    }
    # both sexes
    test_type <- "logistic_regression"
    lreg <- glm(y~x+covar_mat,family = "binomial")
    if(!"x" %in% rownames(summary(lreg)$coefficients)){
      next()
    }
    beta <- summary(lreg)$coefficients["x",1]
    p <- summary(lreg)$coefficients["x",4]
    if(p < 0.01){
      test_type <- "firth_logistic_regression"
      lreg <- logistf(y~x+covar_mat)
      beta <- lreg$coefficients["x"]
      p <- lreg$prob["x"]
      if(p == 0){
        lreg <- glm(y~x+covar_mat,family = "binomial")
        p <- summary(lreg)$coefficients["x",4]
      }
    }
  }
  if(length(levels(factor(y)))>2){
    # not binary. do linear regression
    pheno_type <- "continuous"
    test_type <- "linear_regression"    
    # remove outliers - 5 SDs
    sdy <- sd(y,na.rm=T)
    meany <- mean(y,na.rm=T)
    outliers_i <- union(which(y > meany + 5*sdy),which(y < meany - 5*sdy))
    if(length(outliers_i) > 0){
      y[outliers_i] <- NA
    }
    n_cases <- length(which(!is.na(y)))
    if(n_cases < 100){
      next()
    }
    y_int <- irnt(y) # qnorm((rank(y,na.last="keep")-0.5)/sum(!is.na(y)))
    n_controls <- NA
    # also check distribution of phenotype
    # if number of unique values <=12, skip since it might be categorical
    n_values <- length(unique(y))
    if(n_values <= 12){
      next()
    }
    # both sexes
    beta <- NA
    p <- NA
    lreg <- lm(y~x+covar_mat)
    lreg_int <- lm(y_int~x+covar_mat)
    if(nrow(summary(lreg)$coefficients) == (ncol(covar_mat) + 2)){
      beta <- summary(lreg)$coefficients["x",1]
      p <- summary(lreg)$coefficients["x",4]
      beta_int <- summary(lreg_int)$coefficients["x",1]
      p_int <- summary(lreg_int)$coefficients["x",4]
    }
  }
  res <- data.frame(gene,long_pheno,pheno_name,beta,p,test_type)
  names(res) <- c("gene","long_pheno","pheno","beta","p","test_type")
  phewas_res <- rbind(phewas_res,res)
}

fwrite(phewas_res,paste0("/camhpc/home/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/phewas.maf0.001.",strsplit(strsplit(phenofile,"/")[[1]][7],"[.]")[[1]][1],".",gene,".tsv"),sep="\t")

