library(data.table)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(RNOmni)
library(ggpubr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(cowplot)

load("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/gencode.v38.annotation.exon.rda")

plot.gene.name <- "TTN"
chr <- unique(gene.anno$V1[gene.anno$gene.name==plot.gene.name])
exon.pos <- dplyr::filter(gene.anno,gene.name==plot.gene.name & grepl("Ensembl_canonical",V9))
exon.pos <- dplyr::filter(gene.anno,gene.name==plot.gene.name & grepl("OBSCN-209",V9))
getpid <- function(x){return(strsplit(strsplit(strsplit(x,"; ")[[1]][10],' "')[[1]][2],'\\.')[[1]][1])}
plot.pid <- getpid(exon.pos$V9[1])
gettid <- function(x){return(strsplit(strsplit(strsplit(x,"; ")[[1]][2],' "')[[1]][2],'\\.')[[1]][1])}
plot.tid <- gettid(exon.pos$V9[1])
getEN <- function(x){return(strsplit(strsplit(x,"; ")[[1]][7]," ")[[1]][2])}
exon.pos$EN <- sapply(exon.pos$V9,getEN)
exon.pos <- exon.pos[,c(10,11,4,5)]
colnames(exon.pos) <- c("gene","EN","start","end")
exon.pos$median <- (exon.pos$start + exon.pos$end)/2
exon.pos$EN <- as.numeric(exon.pos$EN)

sv_res <- get(load(paste0("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_max_grip.regenie.PTV.results.single_variant.",plot.gene.name,".rda")))
getchr <- function(x){return(strsplit(x,":")[[1]][1])}
getpos <- function(x){return(strsplit(x,":")[[1]][2])}
getalleles <- function(x){return(strsplit(x,":")[[1]][3])}
sv_res$chr <- sapply(sv_res$Variant,getchr)
sv_res$pos <- sapply(sv_res$Variant,getpos)
sv_res$pos <- as.numeric(sv_res$pos)
sv_res$alleles <- sapply(sv_res$Variant,getalleles)
sv_res$varid <- paste(sv_res$chr,sv_res$pos,sv_res$alleles,sep="_")

clinvar <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/lookup_for_Ben/clinvar.vcf.GRC38.reformatted.withDiseaseID.txt") %>% dplyr::filter(grepl("pathogenic",ignore.case = T,CLNSIG) & grepl("Conflicting",CLNSIG)==F)
clinvar$ID1 <- paste(clinvar$CHR,clinvar$POS,clinvar$REF,clinvar$ALT,sep="_")
clinvar$ID2 <- paste(clinvar$CHR,clinvar$POS,clinvar$ALT,clinvar$REF,sep="_")
clinvar$ID <- ifelse(clinvar$ID1 %in% sv_res$varid,clinvar$ID1,ifelse(clinvar$ID2 %in% sv_res$varid,clinvar$ID2,NA))

sv_res$ClinVar <- ifelse(sv_res$varid %in% as.character(na.omit(clinvar$ID)),"ClinVar Pathogenic/Likely pathogenic","n/a")
sv_res$ClinVar <- factor(sv_res$ClinVar,levels = c("n/a","ClinVar Pathogenic/Likely pathogenic"))
sv_res <- left_join(sv_res,dplyr::select(clinvar,varid=ID,CLNDN),by="varid")

sv_res <- mutate(sv_res,clinvar_pheno=case_when(
  grepl("cardiomyopathy",CLNDN,ignore.case = T) & grepl("dystrophy",CLNDN,ignore.case = T) ~ "Cardiomyopathy and muscular dystrophy",
  grepl("cardiomyopathy",CLNDN,ignore.case = T) & grepl("dystrophy",CLNDN,ignore.case = T)==F ~ "Cardiomyopathy",
  grepl("cardiomyopathy",CLNDN,ignore.case = T)==F & grepl("dystrophy",CLNDN,ignore.case = T) ~ "Muscular dystrophy",
  !is.na(CLNDN) & grepl("cardiomyopathy",CLNDN,ignore.case = T)==F & grepl("dystrophy",CLNDN,ignore.case = T)==F ~ "Other"))
sv_res$clinvar_pheno <- factor(sv_res$clinvar_pheno,levels = c("Cardiomyopathy and muscular dystrophy","Cardiomyopathy","Muscular dystrophy","Other"))

chr_num <- strsplit(chr,"r")[[1]][2]
colnames(sv_res)[13] <- "ClinVar(Pathogenic/Likely Pathogenic)"
dis <- (range(sv_res$pos)[2]-range(sv_res$pos)[1])/50

p1 <- ggplot(sv_res) + 
  geom_point(aes(x=pos,y=beta),shape=20,color="#619CFF") +
  geom_point(data = sv_res[sv_res$ClinVar=="ClinVar Pathogenic/Likely pathogenic",],aes(x=pos,y=beta,shape=`ClinVar(Pathogenic/Likely Pathogenic)`),color="#F8766D") +
  scale_shape_manual(name="ClinVar Pathogenic/Likely pathogenic", values=c(17,18,20,3), labels=c("Cardiomyopathy and muscular dystrophy","Cardiomyopathy","Muscular dystrophy","Other")) +
  geom_line(aes(x=pos,y=beta), stat = 'smooth',method='loess', formula = 'y ~ x',alpha=1,color="red") +
  geom_ribbon(aes(x=pos, y=beta),stat='smooth', method = "loess", se=TRUE, alpha=0.1, formula = 'y ~ x',color=NA) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.6) +
  xlim(178500000,178900000) +
  theme_pubr() +
  labs(x=paste0("Chromosome ",chr_num," position"), y="Effect size on hand grip strength (kg)") +
  scale_y_continuous(limits=c(round(min(sv_res$beta))-40,ceiling(max(sv_res$beta))+15)) +
  theme(legend.justification = c(0,0),legend.direction = "vertical",legend.position = c(0.01,0.8),legend.text = element_text(size=7),legend.title = element_text(size=7)) + 
  guides(shape=guide_legend(nrow=1,byrow=TRUE))
lg1 <- get_legend(p1)

if(nrow(exon.pos)<=10){
  step=1
} else if(nrow(exon.pos)<=50){
  step=5
} else if(nrow(exon.pos)<=200){
  step=10
} else if(nrow(exon.pos)<=500){
  step=50
}

for (i in 1:nrow(exon.pos)){
  p1 <- p1 +	
    geom_segment(x = exon.pos$median[exon.pos$EN==i], y = floor(min(sv_res$beta))-9, xend = exon.pos$median[exon.pos$EN==i], yend = floor(min(sv_res$beta))-12, color="black")
  if(i==1 | i==nrow(exon.pos)){
    p1 <- p1 + annotate(geom="text", x=exon.pos$median[exon.pos$EN==i], y=floor(min(sv_res$beta))-4, label=as.character(i), color="black",size=4)
  } else if(i%%step==0){
    p1 <- p1 + annotate(geom="text", x=exon.pos$median[exon.pos$EN==i], y=floor(min(sv_res$beta))-6.5, label=as.character(i), color="black",size=4)
  }
}
p1 <- p1 + geom_segment(x = exon.pos$median[exon.pos$EN==1], y = floor(min(sv_res$beta))-10.5, xend = exon.pos$median[exon.pos$EN==nrow(exon.pos)], yend = floor(min(sv_res$beta))-10.5, color="black")

edb <- EnsDb.Hsapiens.v86
domain_info <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/protein_domains/InterPro.entry.list")
domains <- as.data.frame(proteins(edb, filter = GeneNameFilter(plot.gene.name), columns = c("tx_id", listColumns(edb, "protein_domain")))) %>% dplyr::filter(protein_id==plot.pid & protein_domain_source=="pfam")
domains <- merge(domains,domain_info,by.x="interpro_accession",by.y="ENTRY_AC",all.x=T)
prngs <- IRanges(start = domains$prot_dom_start, end = domains$prot_dom_end)
names(prngs) <- rep(plot.pid,length(prngs))
res <- proteinToGenome(prngs, edb)

for(i in 1:length(res)){
  domains$exon_s[i] <- min(res[[i]]$exon_rank)
  domains$exon_e[i] <- max(res[[i]]$exon_rank)
  domains$gstart[i] <- exon.pos$median[exon.pos$EN==min(res[[i]]$exon_rank)]
  domains$gend[i] <- exon.pos$median[exon.pos$EN==max(res[[i]]$exon_rank)]
}
domains$`Protein domain` <- factor(domains$ENTRY_NAME,levels = unique(domains$ENTRY_NAME))

p2 <- ggplot(sv_res) + 
  geom_point(aes(x=pos,y=beta),shape=20,color="grey") +
  geom_point(data = sv_res[sv_res$ClinVar=="ClinVar Pathogenic/Likely pathogenic",],aes(x=pos,y=beta),color="midnightblue",size=3) +
  geom_line(aes(x=pos,y=beta), stat = 'smooth',method='loess', formula = 'y ~ x',alpha=1,color="red") +
  geom_ribbon(aes(x=pos, y=beta),stat='smooth', method = "loess", se=TRUE, alpha=0.1, formula = 'y ~ x',color=NA) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.6) +
  theme_pubr() +
  labs(x=paste0("Chromosome ",chr_num," position"), y="Effect size on hand grip strength (kg)") +
  scale_y_continuous(limits=c(round(min(sv_res$beta))-40,ceiling(max(sv_res$beta))+20))
p2 <- p2 + geom_rect(data=domains,aes(xmin=gstart,xmax=gend,fill=`Protein domain`),ymin=floor(min(sv_res$beta))-18,ymax=floor(min(sv_res$beta))-15)
p2 <- p2 + geom_segment(data=domains,aes(x=gstart,xend=gstart,color=`Protein domain`),y=floor(min(sv_res$beta))-18,yend=floor(min(sv_res$beta))-15)
p2 <- p2 + theme(axis.text.y=element_text(size=14),
                 axis.text.x=element_text(size=14),
                 axis.title.y=element_text(size=14),
                 axis.title.x=element_text(size=14),
                 legend.position = c(0.01,0.01),
                 legend.justification = c(0,0),
                 legend.direction = "vertical",legend.text = element_text(size=5),legend.title = element_text(size=7)) + guides(fill=guide_legend(nrow=1,byrow=TRUE))
lg2 <- get_legend(p2)

p_nolg <- ggplot(sv_res) + 
  geom_point(aes(x=pos,y=beta),shape=20,color="#619CFF") +
  geom_point(data = sv_res[sv_res$ClinVar=="ClinVar Pathogenic/Likely pathogenic",],aes(x=pos,y=beta,shape=`ClinVar(Pathogenic/Likely Pathogenic)`),color="#F8766D") +
  scale_shape_manual(name="ClinVar Pathogenic/Likely pathogenic", values=c(17,18,20,3), labels=c("Cardiomyopathy and muscular dystrophy","Cardiomyopathy","Muscular dystrophy","Other")) +
  geom_line(aes(x=pos,y=beta), stat = 'smooth',method='loess', formula = 'y ~ x',alpha=1,color="midnightblue") +
  geom_ribbon(aes(x=pos, y=beta),stat='smooth', method = "loess", se=TRUE, alpha=0.1, formula = 'y ~ x',color=NA) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.6) +
  xlim(178500000,178900000) +
  theme_pubr() +
  labs(x=paste0("Chromosome ",chr_num," position"), y="Effect size on hand grip strength (kg)") +
  scale_y_continuous(limits=c(round(min(sv_res$beta))-40,ceiling(max(sv_res$beta))+20)) + 
  theme(legend.position = "none")
for (i in 1:nrow(exon.pos)){
  p_nolg <- p_nolg +	
    geom_segment(x = exon.pos$median[exon.pos$EN==i], y = floor(min(sv_res$beta))-9, xend = exon.pos$median[exon.pos$EN==i], yend = floor(min(sv_res$beta))-12, color="black")
  if(i==1 | i==nrow(exon.pos)){
    p_nolg <- p_nolg + annotate(geom="text", x=exon.pos$median[exon.pos$EN==i], y=floor(min(sv_res$beta))-4, label=as.character(i), color="black",size=7/.pt)
  } else if(i%%step==0){
    p_nolg <- p_nolg + annotate(geom="text", x=exon.pos$median[exon.pos$EN==i], y=floor(min(sv_res$beta))-6.5, label=as.character(i), color="black",size=7/.pt)
  }
}
p_nolg <- p_nolg + geom_segment(x = exon.pos$median[exon.pos$EN==1], y = floor(min(sv_res$beta))-10.5, xend = exon.pos$median[exon.pos$EN==nrow(exon.pos)], yend = floor(min(sv_res$beta))-10.5, color="black")
p_nolg <- p_nolg + geom_text(x = exon.pos$median[exon.pos$EN==1] + 15000, y = floor(min(sv_res$beta))-10.5, label="TTN", fontface = "italic", hjust="left", size=9/.pt)
p_nolg <- p_nolg + geom_rect(data=domains,aes(xmin=gstart,xmax=gend,fill=`Protein domain`),ymin=floor(min(sv_res$beta))-18,ymax=floor(min(sv_res$beta))-15)
p_nolg <- p_nolg + geom_segment(data=domains,aes(x=gstart,xend=gstart,color=`Protein domain`),y=floor(min(sv_res$beta))-18,yend=floor(min(sv_res$beta))-15)
p_nolg <- p_nolg + theme(axis.text.y=element_text(size=9),
                         axis.text.x=element_text(size=9),
                         axis.title.y=element_text(size=11),
                         axis.title.x=element_text(size=11),
                         legend.position = "none",
                         plot.margin=margin(10,10,10,10))
p_addlg <- ggdraw(p_nolg) + draw_grob(lg1,0.1,0.7,0.15,0.15) + draw_plot(lg2,0.1,0.15,0.15,0.15)

# TTN GWAS locus plot

library(gassocplot2)
library(ieugwasr)
library(susieR)
library(gridExtra)

gwas <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_max_grip.regenie.sumstats") %>% dplyr::filter(A1FREQ>=0.01 & A1FREQ<=0.99 & INFO>=0.8)
ttn_hg38 <- fread("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/UKBB_max_grip.regenie.GWAS.chr2.hg38") %>% dplyr::filter(V2>=178500000 & V2<=178900000)

ttn_hg38 <- inner_join(ttn_hg38,dplyr::select(gwas,V4=ID,ALLELE1,ALLELE0,BETA,SE,P),by="V4")

bfile <- "/camhpc/home/yhuang5/muscle_targets/regenie/susie_TTN/UKBB_maxHGS.TTN"
ld <- ld_matrix_local(variants = ttn_hg38$V4, bfile = bfile, plink_bin = "/camhpc/home/yhuang5/plink1.9/plink", with_alleles = T)
diag(ld) <- 1
ld[is.na(ld)] <- 0

snps <- data.frame("varid"=row.names(ld))
snps$varid <- as.character(snps$varid)
getrs <- function(x){return(strsplit(x,"_")[[1]][1])}
geta1 <- function(x){return(strsplit(x,"_")[[1]][2])}
geta2 <- function(x){return(strsplit(x,"_")[[1]][3])}
snps$rsid <- sapply(snps$varid,getrs)
snps$a1 <- sapply(snps$varid,geta1)
snps$a2 <- sapply(snps$varid,geta2)

ttn_hg38 <- merge(ttn_hg38,dplyr::select(snps,varid,V4=rsid,a1,a2),by="V4")
ttn_hg38 <- dplyr::filter(ttn_hg38,(ALLELE1==a1 & ALLELE0==a2) | (ALLELE1==a2 & ALLELE0==a1))

finalsnps <- ttn_hg38$varid
ld <- ld[finalsnps,finalsnps]
ttn_hg38$BETA.harm <- ifelse(ttn_hg38$ALLELE1==ttn_hg38$a2,ttn_hg38$BETA,-ttn_hg38$BETA)
ttn_hg38$Z.harm <- ttn_hg38$BETA.harm/ttn_hg38$SE
ttn_hg38$log10p <- -log10(ttn_hg38$P)
plot_ttn_hg38 <- dplyr::select(ttn_hg38,marker=V4,chr=V1,pos=V2,z=Z.harm)
plot_ttn_hg38$chr <- 2

chr <- 2
x.min <- 178500000
x.max <- 178900000
gene.region <- gassocplot2::genes_b38[(gassocplot2::genes_b38$chr==chr & !(gassocplot2::genes_b38$end<x.min) & !(gassocplot2::genes_b38$start>x.max)),]
gene.region$start[gene.region$start<x.min] <- x.min
gene.region$end[gene.region$end>x.max] <- x.max

plot_gene_two <- function (gene.region, chr, x.min, x.max, stack = FALSE) 
{
  small.gene <- (as.numeric(gene.region$end) - as.numeric(gene.region$start)) < 
    (x.max - x.min)/190
  gene.region$mid.point <- as.numeric(gene.region$start) + 
    (as.numeric(gene.region$end) - as.numeric(gene.region$start))/2
  gene.region$start[small.gene] <- gene.region$mid.point[small.gene] - 
    (x.max - x.min)/380
  gene.region$end[small.gene] <- gene.region$mid.point[small.gene] + 
    (x.max - x.min)/380
  genes.start.stop <- as.matrix(gene.region[, c("start", "end")])
  genes.df.pos <- data.frame(name = paste0("gene", rep(1:nrow(genes.start.stop), 
                                                       each = 2)), pos = as.vector(t(genes.start.stop)), y = (16 - 
                                                                                                                8 * rep(rep(1:2, each = 2), ceiling(nrow(genes.start.stop)/2))[1:(2 * 
                                                                                                                                                                                    nrow(genes.start.stop))]), stringsAsFactors = F)
  plot.pos <- ggplot(data = genes.df.pos, mapping = aes(x = pos, 
                                                        y = y)) + theme_bw() + xlab(paste0("Position on chromosome ", 
                                                                                           chr)) + ylab(" ") + scale_y_continuous(limits = c(-5, 
                                                                                                                                             17), breaks = c(8, 16), labels = c("      ", "      ")) + 
    theme(axis.title.y = element_text(vjust = 2), axis.title.x = element_text(vjust = -0.5), 
          axis.ticks.y = element_blank(), panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) + scale_x_continuous(limits = c(x.min, 
                                                                              x.max))
  if (stack == TRUE) {
    plot.pos <- plot.pos + theme(axis.title = element_text(size = 11), 
                                 axis.text = element_text(size = 9))
  }
  else {
    plot.pos <- plot.pos + theme(axis.title = element_text(size = 11), 
                                 axis.text = element_text(size = 9))
  }
  plot.genes <- plot.pos + geom_line(data = genes.df.pos, mapping = aes(x = pos, 
                                                                        y = y, group = name), colour = "blue4", size = 0.8)
  genes.df.mid.point <- data.frame(name = gene.region$gene, 
                                   x = as.numeric(gene.region$mid.point), y = (16 - 8 * 
                                                                                 rep(rep(1:2, each = 1), ceiling(nrow(gene.region)/2))[1:nrow(gene.region)] + 
                                                                                 3.6), stringsAsFactors = F)
  plot.genes <- plot.genes + geom_text(data = genes.df.mid.point, 
                                       mapping = aes(x = x, y = y, label = name), color = "black", 
                                       size = 9/.pt, fontface = 3)
  return(plot.genes)
}

p_gwas_genes <- plot_gene_two(gene.region,chr,x.min,x.max)

index_snp <- ttn_hg38$varid[ttn_hg38$P==min(ttn_hg38$P)]
ttn_hg38$ld <- ld[,index_snp]
ttn_hg38$r2 <- (ttn_hg38$ld)^2
ttn_hg38 <- mutate(ttn_hg38,R2 = case_when(r2>0.8 ~ "0.8-1.0",
                                           r2>0.6 & r2<=0.8 ~ "0.6-0.8",
                                           r2>0.4 & r2<=0.6 ~ "0.4-0.6",
                                           r2>0.2 & r2<=0.4 ~ "0.2-0.4",
                                           r2>0 ~ "0.0-0.2"))

res <- susieR::susie_rss(ttn_hg38$Z.harm,ld,L = 10,n = 366307)

p_gwas_snps <- ggplot(data = ttn_hg38,aes(x=V2,y=log10p,fill=R2)) + geom_point(pch=21, size=2) + scale_x_continuous(limits=c(178500000,178900000)) + scale_y_continuous(limits=c(0,18))
p_gwas_snps <- p_gwas_snps + scale_fill_manual(values=c("darkblue", "#66FF66", "#FFCC00", "#FF9933", "#CC3300"), drop=FALSE)
p_gwas_snps <- p_gwas_snps + geom_point(data=dplyr::filter(ttn_hg38,varid==index_snp), aes(x=V2,y=log10p), pch=23, colour="black", fill="purple", size=3)
p_gwas_snps <- p_gwas_snps + theme_pubr() + xlab("Chromosome 2 position") + ylab(expression("-log"["10"]*paste("(",italic("p"),")")))
p_gwas_snps <- p_gwas_snps + theme(axis.title.y=element_text(vjust=2.25, size=11), axis.text=element_text(size=9)) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
p_gwas_snps <- p_gwas_snps + theme(legend.text=element_text(size=7), legend.margin = margin(0.1,0.2,0.1,0.2, unit='cm'), legend.title=element_text(size=9), legend.background = element_rect(colour = "black")) + theme(panel.background=element_rect(fill=NA)) + theme(legend.position="bottom") + guides(fill = guide_legend(nrow = 1))
susie_sets <- rbind(data.frame("set"=rep("L1",length(res$sets$cs$L1)),"snp"=colnames(ld)[res$sets$cs$L1]),
                    data.frame("set"=rep("L2",length(res$sets$cs$L2)),"snp"=colnames(ld)[res$sets$cs$L2]),
                    data.frame("set"=rep("L3",length(res$sets$cs$L3)),"snp"=colnames(ld)[res$sets$cs$L3]))
label_points <- merge(ttn_hg38,susie_sets,by.x="varid",by.y="snp")
p_gwas_snps <- p_gwas_snps + geom_label_repel(data=label_points[label_points$varid==index_snp,],
                                              mapping=aes(x=V2,y=log10p,label=V4),fill="white",
                                              segment.color="black", size=7/.pt, point.padding=0.15, nudge_x=-0.5, nudge_y=5)
p_gwas_snps <- p_gwas_snps + geom_label_repel(data=label_points[label_points$set=="L1" & label_points$varid!=index_snp,],
                                              mapping=aes(x=V2,y=log10p,label=V4),fill="white",
                                              segment.color="black", size=6/.pt, point.padding=0.15, nudge_x=0.5, nudge_y=1,force=1.5,min.segment.length=0)
p_gwas_snps_lg <- g_legend(p_gwas_snps)
p_gwas_snps <- p_gwas_snps + theme(legend.position="none",axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
p_gwas <- grid.arrange(
  grobs = list(p_gwas_snps,p_gwas_genes,p_gwas_snps_lg),
  heights= c(5, 2.5, 1),
)

ggarrange(p_addlg,p_gwas,ncol = 1,nrow = 2,labels = c("a","b"),heights = c(1.5,1.1))
ggsave(filename = paste0(plot.gene.name,".variant_level.PTV.GWAS.v4.pdf"),
       device = "pdf",path = "/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/",width = 180,height = 210,units = "mm",dpi = 300)
