library(dplyr)

Semi_impute_IBD <-function(haplotype,referencehaps,heatmapfile, return, usecore=T){
  
  if(is.null(usecore)|is.na(usecore)){usecore<-T}
  if(class(haplotype)!="data.frame"){
   haplotable<-read.table(haplotype, header=T)}else{haplotable <- haplotype}
  if(missing(heatmapfile)){heatmapfile <- "/tmp/heatmap_IBD.pdf"}
  if(missing(return)){return="unreliable"}
  
 hap <- haplotable[1,"hap"]
  
  haplotable <- haplotable[!duplicated(haplotable$pos),]
  
  
  referencetable <- subset(referencehaps, referencehaps$V2%in%haplotable$pos)
  haplotable <- subset(haplotable, haplotable$pos%in%referencetable$V2)
  
  
  haplotable <- haplotable[order(haplotable$pos),]
  haplotable$var <- as.character(haplotable$var)
  for(i in 1:nrow(haplotable)){
    haplotable[i,"var"] <- strsplit(as.character(haplotable[i,"ID"]), split=":")[[1]][2]
  }
  
  
  allelematching <- data.frame(happos=haplotable$pos, refpos=referencetable$V2, hapallele=haplotable$var, 
                               refref=as.character(referencetable$V3), refalt=as.character(referencetable$V4))
  
  allelematching$match <- ifelse(allelematching$hapallele==allelematching$refref, 0, 
                                 ifelse(allelematching$hapallele==allelematching$refalt, 1, 2))
  haplotable$binary <- allelematching$match
  
  if(usecore==T){
  corehap <- subset(haplotable, haplotable$class=="Double_link"&haplotable$round<21)
  if(nrow(corehap)<100){corehap <- subset(haplotable, haplotable$class=="Double_link")}
  unreliablehap <- subset(haplotable, haplotable$ID%in%corehap$ID==F)
  }
  
  if(usecore==F){
    corehap <- haplotable
    unreliablehap <- NULL
  }
  matcheable <- subset(referencetable, referencetable$V2%in%corehap$pos)
  matcheable[,2] <- as.numeric(as.character(matcheable[,2]))
  
  #reconstruct the IBD blocks from which the core haplotype is most likely constructed
  startingpoint <- 1
  
  binary <- corehap$binary
  IBD_build <- NULL
  unclear_SNPs <- NULL
  while(startingpoint < (length(binary)-4)){
    ## match the first 5 SNPs
    
    word <- binary[startingpoint:(startingpoint+4)]
    
    wordmatch <- NULL
    for(i in 1:ncol(matcheable)){
      column <- matcheable[,i]
      columnword <- column[startingpoint:(startingpoint+4)]
      if(length(which(word==columnword))==length(word)){
        wordmatch <- cbind(wordmatch, matcheable[,i])
        colnames(wordmatch)[ncol(wordmatch)] <- paste(i)
      }}
    
    if(is.null(wordmatch)){
      word <- binary[startingpoint:(startingpoint+3)]
      for(i in 1:ncol(matcheable)){
        column <- matcheable[,i]
        columnword <- column[startingpoint:(startingpoint+3)]
        if(length(which(word==columnword))==length(word)){
          wordmatch <- cbind(wordmatch, matcheable[,i])
          colnames(wordmatch)[ncol(wordmatch)] <- paste(i)
        }}
      
    }
    
    if(is.null(wordmatch)){
      startingpoint <- startingpoint+4
    }else{
      
      scoreframe <- data.frame(refnr=as.numeric(colnames(wordmatch)),firstmatch=startingpoint, matchlength=0, lastmatch=0)
      for(i in 1:nrow(scoreframe)){
        refnr <- scoreframe[i,1]
        column <- matcheable[,refnr]
        errors <- 0
        matchlength <- 0
        startcounting <- startingpoint+5
        lastmatch <- startingpoint+4
        while(errors<=0&startcounting<length(binary)){
          if(column[startcounting]==binary[startcounting]){
            matchlength <- matchlength+1
            lastmatch <- startcounting
            startcounting <- startcounting+1
          }else{errors <- errors+1
          startcounting <- startcounting+1}
        }
        
        scoreframe[i,3] <- matchlength
        scoreframe[i,4] <- lastmatch
        
      }
      
      scoreframe <- subset(scoreframe, scoreframe$matchlength==max(scoreframe$matchlength))
      
      if(nrow(scoreframe)>1){
        matched_rgn <- subset(referencehaps, referencehaps$V2>corehap[scoreframe[1,2],1]&referencehaps$V2<corehap[scoreframe[1,4],1])
        matched_rgn <- matched_rgn[,c(2,scoreframe$refnr)]
        matched_rgn$sum <- rowSums(matched_rgn[,2:ncol(matched_rgn)])
        
        unclear_SNPs <- c(unclear_SNPs, matched_rgn[matched_rgn$sum!=0&matched_rgn$sum!=ncol(matched_rgn)-2,1])
      }
      
      if(is.null(IBD_build)){
      IBD_build <- rbind(IBD_build, scoreframe[1,])
      startingpoint <- IBD_build[nrow(IBD_build),4]-4
      
      
      }else{
      if(scoreframe[1,1]!=IBD_build[nrow(IBD_build),"refnr"]){
        IBD_build <- rbind(IBD_build, scoreframe[1,])
        startingpoint <- scoreframe[1,4]-4}
        
      if(scoreframe[1,1]==IBD_build[nrow(IBD_build),"refnr"]){
         startingpoint <- startingpoint+10}}
      
    }}
  
  IBD_build$matchlength <- IBD_build$matchlength+5
 
  
unclear_SNPs <- unique(unclear_SNPs)
  refmatch <- NULL
  
  IBD_build<- subset(IBD_build, IBD_build$refnr!=0)
  for(i in 1:nrow(IBD_build)){
    block <- IBD_build[i,]
    if(block$refnr!=0&block$matchlength>6){
    start <- corehap[block$firstmatch,1]
    if(i<nrow(IBD_build)){
    end <- corehap[block$lastmatch-4,1]}
    if(i==nrow(IBD_build)){
      end <- corehap[block$lastmatch,1]}
     IBDblock <- referencehaps[referencehaps$V2>=start&referencehaps$V2<end,c(2,block$refnr)]
     colnames(IBDblock) <- c("pos", "cons")
    refmatch <- rbind(refmatch,IBDblock)}
  }
  
  refmatch <- subset(refmatch, refmatch$pos%in%unclear_SNPs==F)
  
  ### make a heatmap
    heatmapframe <- rbind(corehap[,c("pos", "binary")], unreliablehap[,c("pos", "binary")])
    heatmapframe <- heatmapframe[order(heatmapframe$pos),]
  heatmapframe$core <- ifelse(heatmapframe$pos%in%corehap$pos, heatmapframe$binary, 2)
  heatmapframe$unreliable <- ifelse(heatmapframe$pos%in%unreliablehap$pos, heatmapframe$binary, 2)
  
  heatmapframe$cons <- 2
  for(i in 1:nrow(heatmapframe)){
    row <- heatmapframe[i,]
    pos <- row$pos
    if(pos%in%refmatch$pos){heatmapframe[i,"cons"]<-refmatch[refmatch$pos==pos,"cons"]}
  }
  

 trefs <- referencehaps[referencehaps$V2%in%heatmapframe$pos,]
    for(i in 1:nrow(IBD_build)){
      block <- IBD_build[i,]
      if(block$refnr!=0&block$matchlength>6){
        start <- corehap[block$firstmatch,1]
        if(i<nrow(IBD_build)){
          end <- corehap[block$lastmatch-4,1]}
        if(i==nrow(IBD_build)){
          end <- corehap[block$lastmatch,1]}
          IBDblock <- trefs[trefs$V2>=start&trefs$V2<=end,c(2,block$refnr)]
        colnames(IBDblock) <- c("pos", block$refnr)
        IBDblock[,2] <- ifelse(IBDblock$pos%in%unclear_SNPs, 2, IBDblock[,2])
        missing <- data.frame(pos=heatmapframe[heatmapframe$pos%in%IBDblock$pos==F,1], V2=2)
        colnames(missing) <- c("pos", block$refnr)
        IBDblock <- rbind(IBDblock, missing)
        IBDblock <- IBDblock[order(IBDblock$pos),]
       
        
        heatmapframe <- cbind(heatmapframe, IBDblock[,2])
        colnames(heatmapframe)[ncol(heatmapframe)] <- block$refnr
  }}
 
 matr <- data.matrix(heatmapframe[,2:ncol(heatmapframe)]) 
 matr <- na.omit(matr)
 
 matching <- nrow(heatmapframe[(heatmapframe$unreliable==1&heatmapframe$cons==1)|(heatmapframe$unreliable==0&heatmapframe$cons==0),])
 outof <- nrow(heatmapframe[heatmapframe$unreliable!=2&heatmapframe$cons!=2,])
 
 row.names(matr) <- heatmapframe$pos
  colnames(matr)[1] <- "all"
  pdf(file=heatmapfile, useDingbats = F)
  heatmap.2(matr, breaks=c(0, 0.5,1.5, 2.5, 3.5, 4.5), col = c("blue", "green", "white", "orange", "red"),
            dendrogram = "none",
           Colv=F, Rowv=F,labRow=T, trace="none", 
           ylab=paste0("SNPs chr",trefs[1,1],":",min(rownames(matr)), "-", max(rownames(matr))), 
           cexCol=1, labCol=colnames(matr) ,            
           key=F, main=paste0("matching ", matching, "/",outof, "(", round((matching/outof)*100,2),"%)"))
  dev.off()
  
  heatmap.2(matr, breaks=c(0, 0.5,1.5, 2.5, 3.5, 4.5), col = c("blue", "green", "white", "orange", "red"),
            dendrogram = "none",
            Colv=F, Rowv=F,labRow=T, trace="none", 
            ylab=paste0("SNPs chr",trefs[1,1],":",min(rownames(matr)), "-", max(rownames(matr))), 
            cexCol=1, labCol=colnames(matr) ,            
            key=F, main=paste0("matching ", matching, "/",outof, "(", round((matching/outof)*100,2),"%)"))
  
  if(usecore==T){
  unreliablehap$IBD <- 2
  
  
  for(i in 1:nrow(unreliablehap)){
    pos <- unreliablehap[i,1]
    if(pos%in%refmatch$pos){
      unreliablehap[i,"IBD"] <- refmatch[refmatch$pos==pos,2]}
  }}
  
  
  if(return=="unreliable"){
  return(unreliablehap[unreliablehap$binary!=unreliablehap$IBD,])}
  
  if(return=="refhap"){
    refmatch <- subset(refmatch, refmatch$pos%in%unclear_SNPs==F)
    
    refmatch$ref <- "X"
    refmatch$alt <- "X"
    refmatch$core <- "X"
    for(i in 1:nrow(refmatch)){
      refmatch[i,3] <- as.character(referencehaps[referencehaps$V2==refmatch[i,1],3])
      refmatch[i,4] <- as.character(referencehaps[referencehaps$V2==refmatch[i,1],4])
      if(refmatch[i,1]%in%corehap$pos){
        refmatch[i,5] <- corehap[corehap$pos==refmatch[i,1], "var"]
      }}
      return(refmatch)
  }
  if(return=="matrix"){
    return(matr)}
  if(return=="matrix2"){
    return(matr2)}

}

# haplotype <- hap1
# 
# haplotypes <- read.table("~/resources/Small files/haps/HBB/final/family_HBB12_M_noimp.haps.txt", header=T)
#  hap1 <- haplotypes[haplotypes$hap==1,]
# 
# refmatch <- Semi_impute_IBD(hap1, HBB_referencehaps, return="refhap", usecore=T)
# plaatje2 <- Semi_impute_IBD(hap2, HBB_referencehaps, return="matrix", usecore=T)

# fams <- c(5,6,7,8,9,10,12)
# for(i in 1:length(fams)){
#   familynr <- fams[i]
# phaps <- read.table(paste0("~/resources/Small files/haps/HBB/final/family_HBB",familynr,"_D_noimp.haps.txt"), header=T)
# phap1 <- phaps[phaps$hap==1,]
# phap2 <- phaps[phaps$hap==2,]
# mhaps <- read.table(paste0("~/resources/Small files/haps/HBB/final/family_HBB",familynr,"_M_noimp.haps.txt"), header=T)
# mhap1 <- mhaps[mhaps$hap==1,]
# mhap2 <- mhaps[mhaps$hap==2,]
# 
# phap1_mat <- data.frame(Semi_impute_IBD(phap1, HBB_referencehaps,paste0("~/resources/pdfs/imp_heatmaps/family",familynr,"_D1_IBD.pdf"),return="matrix", usecore=T))
# phap2_mat <- data.frame(Semi_impute_IBD(phap2, HBB_referencehaps,paste0("~/resources/pdfs/imp_heatmaps/family",familynr,"_D2_IBD.pdf"),return="matrix", usecore=T))
# mhap1_mat <- data.frame(Semi_impute_IBD(mhap1, HBB_referencehaps,paste0("~/resources/pdfs/imp_heatmaps/family",familynr,"_M1_IBD.pdf"),return="matrix", usecore=T))
# mhap2_mat <- data.frame(Semi_impute_IBD(mhap2, HBB_referencehaps,paste0("~/resources/pdfs/imp_heatmaps/family",familynr,"_M2_IBD.pdf"),return="matrix", usecore=T))
# 
# matchingD1 <- nrow(phap1_mat[(phap1_mat$unreliable==1&phap1_mat$cons==1)|(phap1_mat$unreliable==0&phap1_mat$cons==0),])
# outofD1 <- nrow(phap1_mat[phap1_mat$unreliable!=2&phap1_mat$cons!=2,])
# percentageD1 <- round((matchingD1/outofD1)*100,1)
# 
# matchingD2 <- nrow(phap2_mat[(phap2_mat$unreliable==1&phap2_mat$cons==1)|(phap2_mat$unreliable==0&phap2_mat$cons==0),])
# outofD2 <- nrow(phap2_mat[phap2_mat$unreliable!=2&phap2_mat$cons!=2,])
# percentageD2 <- round((matchingD2/outofD2)*100,1)
# 
# matchingM1 <- nrow(mhap1_mat[(mhap1_mat$unreliable==1&mhap1_mat$cons==1)|(mhap1_mat$unreliable==0&mhap1_mat$cons==0),])
# outofM1 <- nrow(mhap1_mat[mhap1_mat$unreliable!=2&mhap1_mat$cons!=2,])
# percentageM1 <- round((matchingM1/outofM1)*100,1)
# 
# matchingM2 <- nrow(mhap2_mat[(mhap2_mat$unreliable==1&mhap2_mat$cons==1)|(mhap2_mat$unreliable==0&mhap2_mat$cons==0),])
# outofM2 <- nrow(mhap2_mat[mhap2_mat$unreliable!=2&mhap2_mat$cons!=2,])
# percentageM2 <- round((matchingM2/outofM2)*100,1)
# 
# data <- data.frame(names=c(paste0("D", familynr, ".1"),paste0("D", familynr, ".2"), paste0("M", familynr, ".1"),paste0("M", familynr, ".2")),
#                    matching=c(matchingD1,matchingD2,matchingM1,matchingM2),
#                    outof=c(outofD1,outofD2,outofM1,outofM2),
#                    percentages=c(percentageD1,percentageD2,percentageM1, percentageM2))
# if(exists("matchdata")==F){
#   matchdata <- data
# }else{matchdata <- rbind(matchdata, data)}
# 
# 
# }