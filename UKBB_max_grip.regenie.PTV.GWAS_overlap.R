library(data.table)
library(dplyr)

gencode <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/ALS_progression/GTAC_ChiaYen/gencode.v39.annotation.gtf.gz",skip = 1,header = F)
gwas_catalog <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/efotraits_EFO_0006941-associations-2022-10-23.csv")
gwas_catalog <- filter(gwas_catalog,as.numeric(RAF)>=0.01)
getchr <- function(x){return(strsplit(x,":")[[1]][1])}
getpos <- function(x){return(strsplit(x,":")[[1]][2])}
gwas_catalog$chr <- sapply(gwas_catalog$Location,getchr)
gwas_catalog$pos <- sapply(gwas_catalog$Location,getpos)
gwas_catalog$pos <- as.numeric(gwas_catalog$pos)
gwas_catalog$start <- gwas_catalog$pos-500000
gwas_catalog$end <- gwas_catalog$pos+500000

gwas_catalog$genes_within500kb <- NA
for(i in 1:nrow(gwas_catalog)){
  sub <- filter(gencode,V1==paste0("chr",gwas_catalog$chr[i]) & !(V4>gwas_catalog$end[i] | V5<gwas_catalog$start[i]) & V3=="gene" & grepl("protein_coding",V9))
  if(nrow(sub)==0){next}
  genes <- c()
  for(j in 1:nrow(sub)){
    terms <- strsplit(sub$V9[j],"; ")[[1]]
    genes <- c(genes,strsplit(terms[grepl("gene_name",terms)],'\\"')[[1]][2])
  }
  gwas_catalog$genes_within500kb[i] <- paste(genes,collapse = ";")
  if(ceiling(i/20)==(i/20)){print(i)}
}

pchange <- function(x){return(paste0(strsplit(x," ")[[1]][1],"e-",strsplit(strsplit(x," ")[[1]][3],"-")[[1]][2]))}
gwas_catalog$`P-value` <- as.numeric(sapply(gwas_catalog$`P-value`,pchange))

ptv_burden <- fread("/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.updated.final.tsv")
ptv_burden_gwas_catalog <- filter(ptv_burden,GENE %in% unique(strsplit(paste(gwas_catalog$genes_within500kb,collapse = ";"),";")[[1]]))

fwrite(gwas_catalog,file="/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/GWAS_catalog.HGS.loci.maf0.01.tsv",sep="\t")
fwrite(ptv_burden_gwas_catalog,file="/edgehpc/dept/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/UKBB_max_grip.regenie.PTV.GWAS_catalog.HGS.loci.maf0.01.500kb.tsv",sep="\t")
