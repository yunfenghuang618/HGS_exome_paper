library(icd.data)
library(data.table)

get_icd10_chapter <- function(icd_code){
  my_letters <- toupper(letters[1:26])
  first_letter <- which(my_letters == substr(icd_code,1,1))
  first_number <- as.integer(substr(icd_code,2,3))
  if(first_letter == 21){
    ch_name <- "ICD10: Codes for special purposes"
    return(ch_name)
  }
  for(i in 1:length(icd10_chapters)){
    ch_name <- paste("ICD10",names(icd10_chapters[i]))
    ch <- data.frame(icd10_chapters[i])
    start <- ch[,1][1]
    stop <- ch[,1][2]
    start_letter <- which(my_letters == substr(start,1,1))
    start_number <- as.integer(substr(start,2,3))
    stop_letter <- which(my_letters == substr(stop,1,1))
    stop_number <- as.integer(substr(stop,2,3))
    if(start == "F01"){
      start_number = 0
    }
    if(stop == "T88"){
      stop_number <- 99
    }
    if(stop == "O9A"){
      stop_number <- 99
    }
    # letter is completly in between start and stop
    if(first_letter >= start_letter && first_letter <= stop_letter){
      return(ch_name)
    }
    # first letter completely matches start and stop letters, number is in between
    if(first_letter == start_letter && first_letter == stop_letter){
      if(first_number >= start_number && first_number <= stop_number){
        return(ch_name)
      }
    }
    # first letter matches start, number greater than start, first letter does not match stop
    if(first_letter == start_letter && first_letter != stop_letter){
      if(first_number >= start_number){
        return(ch_name)
      }
    }
    # first_letter does not match start, first letter matches stop, number less than stop
    if(first_letter != start_letter && first_letter == stop_letter){
      if(first_number <= stop_number){
        return(ch_name)
      }
    }
  }
  return("ICD10_not_found")
}

do_genes <- c("TTN")

for(gene in do_genes){
  print(gene)
  var_info <- read.csv("/home/cchen11/ukbb_cogexome/data/pheno_phewas_v1/Data_Dictionary_Showcase.csv")
  var_levels <- strsplit(as.character(var_info$Path)," > ")
  
  auto_cats <- read.csv("/home/cchen11/ukbb_cogexome/data/pheno_phewas_v1/phewas_pheno_categories.csv")
  
  phewas <- fread(paste("/home/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/phewas.maf0.001.",gene,".tsv",sep=""))
#  exc_i <- which(phewas$pheno %in% c("AGE_AT_DEATH","DATE_OF_DEATH","DATE_OF_DEATH_INT","CAUSE_DEATH_ICD10_CODES_INT"))
#  phewas <- phewas[-exc_i,]
  # get categories
  categories <- data.frame()
  unmatched <- c()
  for(i in 1:nrow(phewas)){
    pheno <- as.character(phewas$pheno[i])
    # f.n.0.0
    if(strsplit(pheno,"\\.")[[1]][1] == "f"){
      field_id <- strsplit(pheno,"\\.")[[1]][2]
      if(field_id %in% var_info$FieldID){
        all_levels <- var_levels[which(var_info$FieldID == field_id)][[1]]
        if(length(all_levels) == 2){
          categ <- as.character(var_levels[which(var_info$FieldID == field_id)][[1]][2])
        }
        if(length(all_levels) > 2){
          categ <- as.character(var_levels[which(var_info$FieldID == field_id)][[1]][2])
        }
        # categories <- append(categories,categ)
        categories <- rbind(categories,data.frame(pheno,categ))
        next()
      }
      else{
        #categ <- "Other"
        #categories <- append(categories,categ)
        #categories <- rbind(categories,data.frame(pheno,categ))
        next()
      }
    }
    if(strsplit(pheno,"_")[[1]][1] == "FH"){
      categ <- "Family history"
      #categories <- append(categories,categ)
      categories <- rbind(categories,data.frame(pheno,categ))
      next()
    }
    if(pheno %in% auto_cats$pheno){
      categ <- as.character(auto_cats$category[which(auto_cats$pheno == pheno)])
      if(categ == "ICD10"){
        # get chapter
        categ <- get_icd10_chapter(pheno)
      }
      # categories <- append(categories,categ)
      categories <- rbind(categories,data.frame(pheno,categ))
      next()
    }
    # categories <- append(categories, "Other")
    categ <- "Composite phenotypes"
    categories <- rbind(categories,data.frame(pheno,categ))
    unmatched <- append(unmatched,pheno)
  }
  
  
  phewas <- merge(categories,phewas,by = "pheno")
  phewas[,2] <- sapply(phewas[,2],as.character)
  phewas[,4] <- sapply(phewas[,4],as.character)
  
  # manual adjustments
  # combined cognitive function with cognitive function online
  phewas$categ[which(phewas$categ == "Cognitive function online")] <- "Cognitive function"
  # move categories with < 30 phenotypes to Other
  categs <- sort(as.character(unique(phewas$categ)))
  for(i in 1:length(categs)){
    categ <- categs[i]
    categ_i <- which(phewas$categ == categ)
    if(length(categ_i) < 30){
      phewas$categ[categ_i] <- "Other"
    }
  }
  
  phewas <- phewas[with(phewas,order(categ,long_pheno)),]
  # exclude physical activity measurements
  phewas <- phewas[-which(phewas$categ == "Physical activity measurement"),]
  categs <- sort(as.character(unique(phewas$categ)))
  # plot_cols <- rainbow(length(categs),s = 0.75)
  phewas$p <- as.numeric(phewas$p)
  pdf(paste("/home/yhuang5/muscle_targets/regenie/ukbb_exome_scripts/phewas_allpheno.",gene,".v2.pdf",sep=""),
      width = 18,
      height = 12)
  
  par(xpd = T)
  par(mar = c(20,4,4,4))
  if(gene == "LDLR"){
    phewas$p[which(phewas$pheno == "self_noncancer_1473.0")] <- 1.52774e-19
  }
  ymax <- max(-log10(phewas$p),na.rm=T)
  if(ymax == Inf){
    ymax <- 40
    phewas$p[which(phewas$p == 0)] <- 40
  }
  if(ymax < -log10(1e-5)){
    ymax <- 10
  }
  
  plot(c(1,nrow(phewas)),c(0,ymax),
       type = "n",
       axes = F,
       xlab = "",
       ylab = "-log10(p)",
       main = paste(gene,"PTV burden PheWAS"))
  y_pos <- -ymax/30
  xpoints <- c()
  for(i in 1:length(categs)){
    plot_col <- "#2573ba"
    if(i %% 2 == 0){
      plot_col <- "#6dad46"
    }
    categ <- categs[i]
    phewas_i <- phewas[which(phewas$categ == categ),]
    points(which(phewas$categ == categ),-log10(phewas_i$p),
           pch = 20,
           col = plot_col)
    xpoint <- median(which(phewas$categ == categ))
    text(xpoint,y_pos,
         categ,
         pos = 2,
         col = plot_col,
         srt = 45,
         cex = 0.75)
    #lines(c(xpoint,xpoint),c(0,y_pos),
    #      col = plot_col)
  }
  
  # label top 15 hits
  ygap <- ymax/10
  tophits <- phewas[with(phewas,order(p)),][1:15,]
  seen_phenos <- c()
  for(i in 1:nrow(tophits)){
    long_pheno <- tophits$long_pheno[i]
    pheno <- tophits$pheno[i]
    if(long_pheno %in% seen_phenos){
      next()
    }
    seen_phenos <- append(seen_phenos,long_pheno)
    xpos <- which(phewas$pheno == pheno)
    ypos <- -log10(tophits$p[i])
    if(i > 1){
      if(abs(-log10(tophits$p[i]) - -log10(tophits$p[i-1])) < ygap){
        if(abs(which(phewas$pheno == pheno) - which(phewas$pheno == tophits$pheno[i-1])) < 500){
          next()
        }
      }
    }
    if(tophits$p[i] > 1e-5){
      next()
    }
    text(xpos,ypos,long_pheno,
         pos = 4,
         cex = 1)
  }
  
  axis(2)
  par(xpd = F)
  abline(h = -log10(1e-5),
         col = "red",
         lty = "dashed")
  dev.off()
}
