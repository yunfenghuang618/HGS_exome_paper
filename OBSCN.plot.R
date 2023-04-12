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

plot.gene.name <- "OBSCN"
chr <- unique(gene.anno$V1[gene.anno$gene.name==plot.gene.name])

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

chr_num <- strsplit(chr,"r")[[1]][2]
dis <- (range(sv_res$pos)[2]-range(sv_res$pos)[1])/50
dot.label <- data.frame("pos"=min(sv_res$pos)-25000,"beta"=ceiling(max(sv_res$beta))+10)

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

OBSCN_var_carrier <- dplyr::filter(sv_res,pos %in% c(228316782,228273387,228274842,228371292,228306921,228288084,228212169,228212170,228353063))

p_nolg <- ggplot(sv_res) + 
  geom_point(aes(x=pos,y=beta),shape=20,color="#619CFF") +
  geom_point(data = OBSCN_var_carrier,aes(x=pos,y=beta),shape=15,color="darkred") +
  geom_point(data = sv_res[sv_res$ClinVar=="ClinVar Pathogenic/Likely pathogenic",],aes(x=pos,y=beta),shape=18,color="#F8766D") +
  geom_point(data = dot.label,aes(x=pos,y=beta),shape=18,color="#F8766D") +
  geom_label(x = min(sv_res$pos)-25000+dis, y=ceiling(max(sv_res$beta))+10,label="ClinVar Pathogenic/Likely pathogenic",hjust="left",color="#F8766D",label.size = 0) +
  geom_line(aes(x=pos,y=beta), stat = 'smooth',method='loess', formula = 'y ~ x',alpha=1,color="midnightblue") +
  geom_ribbon(aes(x=pos, y=beta),stat='smooth', method = "loess", se=TRUE, alpha=0.1, formula = 'y ~ x',color=NA) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.6) +
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
p_addlg <- ggdraw(p_nolg) + draw_plot(lg2,0.1,0.1,0.2,0.2)

pdf(file=paste0("/mnt/depts/dept04/compbio/human_genetics/users/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/",plot.gene.name,".variant_level.PTV.v4.pdf"),bg = "white",width = 8.2,height = 11.6)
print(p_addlg)
dev.off()
