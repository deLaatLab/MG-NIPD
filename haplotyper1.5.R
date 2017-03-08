library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)

####changelog
#1.1 added while-loop and iteration counter
#1.2 added seed-haplotype protocol
#1.3 -added auto-viewpoint selection
#     -added sif file creation from scratch
#.3: fixed bug, reduced first round stringency <<-- Improve!!!!
#- added weight/ambig recalculation
#clear previous data





## I assume there is a .hap file present in the directory with the same base name
#basename <- "HBB-D11_lenient_NS"
## the extension is used between the basename and the outputs
#extension <- "_manual"

## chromosome has to be manually specified
#chromosome <- "chr11"

haplotyper1.5 <- function(directory, basename, extension, chromosome){
  setwd(directory)
  badviewpoint <- 0
  viewpoint_override <- 0
  write <- TRUE 
  
  #define outputs
  haps <- paste(basename, ".hap", sep="")
  outputname <- paste(basename, extension, ".sif", sep="")
  hap1_name <- paste(basename, extension, ".hap1.txt", sep="" )
  hap2_name <- paste(basename, extension, ".hap2.txt", sep="" )
  
  haps <- read.table(haps)
  allhaps<- haps
  haps <- haps[haps$V2!=badviewpoint&haps$V4!=badviewpoint,]
  hapfile <- haps
  colnames(hapfile) <- c("var1", "pos1", "var2", "pos2", "weight")
  
  #########recalculate strenght and ambiguity with fixed haplotypes
  recalculate_weight <- function(hap){
    for (i in 1:nrow(hap)){
      ID <- as.character(hap[i,3])
      links <- subset(haps, haps$ID1==ID)
      haplinks <- subset(links, links$ID2%in%hap$ID)
      hap[i,4] <- sum(haplinks$V5) 
      
      return(hap)
    } 
    
  }
  
  
  recalculate_ambig <- function(hap, ambighap){
    for (i in 1:nrow(hap)){
      ID <- as.character(hap[i,3])
      links <- subset(haps, haps$ID1==ID)
      ambiglinks <- subset(links, links$ID2%in%ambighap$ID)
      hap[i,5] <- sum(ambiglinks$V5) 
      
      return(hap)
    } 
    
  }
  
  
  ### autodefine best viewpoint
  
  poslist <- hapfile[!duplicated(hapfile$pos1),]
  poslist <- poslist$pos1
  
  weight_table <- data.frame(pos=numeric(), weight=numeric())
  
  for(i in 1:length(poslist)){
    pos <- poslist[i]
    occurences <- subset(hapfile, hapfile$pos1==pos)
    weight <- sum(occurences$weight)
    line <- data.frame(pos=pos, weight=weight)
    weight_table[i,]<- line
  }
  weight_table <- weight_table[order(weight_table$weight, decreasing = T),]
  #weight_table <- weight_table[weight_table$pos<5200000&weight_table$pos<5300000,]
  
  viewpointpos <- weight_table[1,1]
  if(viewpoint_override>1){viewpointpos <- viewpoint_override}
  viewpointvars <- subset(hapfile, hapfile$pos1==viewpointpos)
  viewpointvars <- unique(viewpointvars$var1)
  hap1allele <- as.character(viewpointvars[1])
  hap2allele <- as.character(viewpointvars[2])
  
  
  iterations <- 25
  debug_iteration <- 4
  
  #############.hap to .sif
  allsif <- haps
  allsif$node1 <- paste(allsif$V2,allsif$V1, sep=":")
  allsif$link <- paste("pp")
  allsif$node2 <- paste(allsif$V4, allsif$V3, sep=":")
  allsif <- allsif[,6:8]
  allsif$used <- paste("X")
  allsif$linkID <- paste(allsif$node1, allsif$node2, sep=":")
  
  ##############
  
  hapfile$ID1 <- paste(hapfile[,2], hapfile[,1], sep=":")
  hapfile$ID2 <- paste(hapfile[,4], hapfile[,3], sep=":")
  
  IDlist <- unique(c(hapfile$pos1, hapfile$pos2))
  IDlist <- as.numeric(subset(IDlist, IDlist!=viewpointpos))
  
  hap1 <- data.frame(pos=viewpointpos, var=hap1allele, 
                     ID=paste(viewpointpos, hap1allele, sep=":"), 
                     weight=10000, ambig=0, class="VP", ratio=0)
  hap1$class <- as.character(hap1$class)
  
  hap2 <- data.frame(pos=viewpointpos, var=hap2allele, 
                     ID=paste(viewpointpos, hap2allele, sep=":"), 
                     weight=10000, ambig=0, class="VP", ratio=0)
  hap2$class <- as.character(hap2$class)
  
  loopnr <- 1
  
  loopstats <- data.frame(loop=numeric(), linked=numeric(), 
                          totalweight=numeric(), totalambig=numeric(),
                          bad=numeric(), isolated=numeric(),
                          selected=numeric())
  
  
  
  seedhap1 <- hap1
  seedhap2 <- hap2
  
  #hapfile<- hapfile[hapfile$pos1==5472343|hapfile$pos2==5472343,]
  
  ###########Iteration 1
  while(loopnr <= iterations){
    
    badsnps <- data.frame(pos=factor(), reason=character())
    
    for(i in 1:length(IDlist)){
      
      SNP <- IDlist[i]
      #SNP <- 117475328
      links <- subset(hapfile, hapfile$pos1==SNP)
      hap_links <- links[links$pos2%in%seedhap1$pos,]
      
      #check if both alleles represented
      
      #first if-else statements calls the ugly SNPs
      if(length(unique(links$ID1))==2){
        var1 <- unique(links$ID1)[1]
        var2 <- unique(links$ID1)[2]
        
        #score for possibility 1: allele/hap 1/1+2/2 and possibility 2: 1/2 2/1
        var1_links <- hap_links[hap_links$ID1%in%var1,]
        var2_links <- hap_links[hap_links$ID1%in%var2,]
        
        #Possibility1
        var1_hap1 <- var1_links[var1_links$ID2%in%seedhap1$ID,]
        var2_hap2 <- var2_links[var2_links$ID2%in%seedhap2$ID,]
        P1 <- sum(var1_hap1$weight, var2_hap2$weight)
        
        P1_links <- rbind(var1_hap1, var2_hap2)
        P1_links$linkID <- paste(P1_links$ID1, P1_links$ID2, sep=":")
        
        #possibility2
        var1_hap2 <- var1_links[var1_links$ID2%in%seedhap2$ID,]
        var2_hap1 <- var2_links[var2_links$ID2%in%seedhap1$ID,]
        
        P2 <- sum(var1_hap2$weight, var2_hap1$weight)
        
        P2_links <- rbind(var1_hap2, var2_hap1)
        P2_links$linkID <- paste(P2_links$ID1, P2_links$ID2, sep=":")
        
        var1_line <- data.frame(pos=SNP, var=var1_links[1,1], ID=var1)
        var2_line <- data.frame(pos=SNP, var=var2_links[1,1], ID=var2)
        
        #second if-else statement continues script if links with a haplotype have been found
        if(P1==0&P2==0){
          badline <- data.frame(pos=paste(SNP), reason="isolated")
          badsnps <- rbind(badsnps, badline)
        }else{
          #if var 1 is in hap 1
          if(P1>P2){
            Weight <- P1
            Ambig <- P2
            if(sum(var1_hap1$weight)>0&sum(var2_hap2$weight)>0) 
            {Class <- "Double_link"} else {Class <- "Single_link" } 
            var1_line$weight=Weight
            var1_line$ambig=Ambig
            var1_line$class=Class
            var1_line$ratio=Ambig/Weight
            
            var2_line$weight=Weight
            var2_line$ambig=Ambig
            var2_line$class=Class
            var2_line$ratio=Ambig/Weight
            
            hap1 <- rbind(hap1,var1_line)
            hap2 <- rbind(hap2, var2_line)
            
            allsif$used <- ifelse(allsif$linkID%in%P1_links$linkID==TRUE, loopnr, allsif$used)
            
          } 
          #if var 1 is in hap2
          else { 
            Weight <- P2
            Ambig <- P1
            if(sum(var1_hap2$weight)>0&sum(var2_hap1$weight)>0) 
            {Class <- "Double_link"} else {Class <- "Single_link" }    
            var1_line$weight=Weight
            var1_line$ambig=Ambig
            var1_line$class=Class
            var1_line$ratio=Ambig/Weight
            var2_line$weight=Weight
            var2_line$ambig=Ambig
            var2_line$class=Class
            var2_line$ratio=Ambig/Weight
            hap1<- rbind(hap1, var2_line)
            hap2<- rbind(hap2, var1_line)
            
            allsif$used <- ifelse(allsif$linkID%in%P2_links$linkID==TRUE, loopnr, allsif$used)}
        }
        
        #end of first if-else statement
      }
      else{ 
        badline <- data.frame(pos=paste(SNP), reason="single_allele")
        badsnps <- rbind(badsnps, badline)}
      #close the for-loop
    }
    
    loopline <- data.frame(loop=loopnr, linked=length(hap1[,1]), 
                           totalweight=sum(hap1$weight), totalambig=sum(hap1$ambig),
                           bad=length(badsnps[,1]), 
                           isolated=length(subset(badsnps, badsnps$reason=="isolated")[,1]),
                           selected=0)
    loopstats[loopnr,] <- loopline[1,]
    
    
    if(loopnr == debug_iteration){ 
      debughap1 <- hap1
      debughap2 <- hap2}
    
    ## add ratio criterion
    
    
    hap1 <- recalculate_weight(hap1)
    hap1 <- recalculate_ambig(hap1, hap2)
    hap2 <- recalculate_weight(hap2)
    hap2 <- recalculate_ambig(hap2, hap1)
    
    hap1[hap1$class=="VP",4]<-1000
    hap2[hap2$class=="VP",4]<-1000
    #first rounds: decreasing quantiles
    #last 3 rounds; all interactions with weight > 1 and no ambiguity
    #last 2 rounds: all interactions
    
    
    
    if(loopnr<5){
      #quantileweight <- quantile(hap1[2:length(hap1[,1]),4], (1-(loopnr*0.2)), names=F)}
      quantileweight <- 30-(5*loopnr)
      
      #hap1$ambig==0
      
      hap1 <- subset(hap1, hap1$weight>quantileweight&hap1$class!="Single_link")
      hap2 <- subset(hap2, hap2$weight>quantileweight&hap2$class!="Single_link")
      hap1 <- subset(hap1, hap1$ratio<0.1)
      hap2 <- subset(hap2, hap2$ratio<0.1)
      
    }
    
    if(loopnr<iterations-4){
      hap1 <- subset(hap1, hap1$class!="Single_link")
      hap2 <- subset(hap2, hap2$class!="Single_link")
      hap1 <- subset(hap1, hap1$ratio<0.3)
      hap2 <- subset(hap2, hap2$ratio<0.3)
    }
    
    
    
    loopstats[loopnr,7] <- length(hap1[,1])
    
    
    #update haps
    seedhap1 <- hap1
    seedhap2 <- hap2
    
    #update IDlist
    IDlist <- IDlist[!IDlist%in%seedhap1$pos]
    
    
    
    
    loopnr <- loopnr+1
    #close the while loop
    
  }
  
  
  #start working on bad snps
  #remove isolated SNPs
  badsnps <- subset(badsnps, badsnps$reason!="isolated")
  badsnps$reason <- as.character(badsnps$reason)
  singles <- length(badsnps$pos)
  
  for(i in 1:length(badsnps[,1])){
    SNP <- badsnps$pos[i]
    links <- subset(hapfile, hapfile$pos1==SNP)
    hap_links <- links[links$pos2%in%hap1$pos,]
    # 3 possibilities: 
    #the found variant links to 1 (P1) 
    #the found variant links to 2 (P2)
    #the found variant links to both
    #it doesn't link at all
    hap1_links <- hap_links[hap_links$ID2%in%hap1$ID,]
    hap2_links <- hap_links[hap_links$ID2%in%hap2$ID,]  
    
    P1 <- sum(hap1_links$weight)
    P2 <- sum(hap2_links$weight)
    
    if(P1>0&P2==0){
      hapline <- data.frame(pos=SNP, var=hap_links[1,1], 
                            ID=paste(SNP,hap_links[1,1], sep=":"),
                            weight=P1, ambig=P2, class="Single_linked_bad", ratio=P2/P1)
      hap1<- rbind(hap1, hapline)
      badsnps[badsnps$pos==SNP,2] <- "Single_linked"}
    else{
      if(P2>0&P1==0){
        hapline <- data.frame(pos=SNP, var=hap_links[1,1], 
                              ID=paste(SNP,hap_links[1,1], sep=":"),
                              weight=P2, ambig=P1, class="Single_linked_bad", ratio=P1/P2)
        hap2<- rbind(hap2, hapline)
        badsnps[badsnps$pos==SNP,2] <- "Single_linked"}
      else{ 
        if(P1>0&P2>0){badsnps[badsnps$pos==SNP,2] <- "Single_ambiguous"}
        else{badsnps[badsnps$pos==SNP,2] <- "Single_isolated"}
      }
    }
    
    #close the for loop
  }
  
  ### filter bad snps
  
  ambigs <- subset(hap1, hap1$ambig!=0)
  ambigs$ratio <- ambigs$ambig/ambigs$weight
  
  
  
  hap1 <- na.omit(subset(hap1, hap1$ratio<0.5))
  hap2 <- na.omit(subset(hap2, hap2$ratio<0.5))
  
  
  
  
  hap1_recalc <- recalculate_weight(hap1)
  hap1_recalc <- recalculate_ambig(hap1_recalc, hap2)
  hap2_recalc <- recalculate_weight(hap2)
  hap2_recalc <- recalculate_ambig(hap2_recalc, hap1)
  
  outputhap1 <- hap1[,1:2]
  outputhap2 <- hap2[,1:2]
  
  ranges.hap1 <- GRanges(seqnames=Rle(chromosome), ranges=as.numeric(outputhap1[,1]))
  bases.hap1 <- data.frame(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, ranges.hap1),stringsAscharacters=T)
  outputhap1$ref <- as.character(bases.hap1$value)
  outputhap1$hap1ref <- ifelse(outputhap1$var==outputhap1$ref, 1, 0)
  colnames(outputhap1) <- c("pos", "seq", "refallele", "hap1ref")
  
  ranges.hap2 <- GRanges(seqnames=Rle(chromosome), ranges=as.numeric(outputhap2[,1]))
  bases.hap2 <- data.frame(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, ranges.hap2),stringsAscharacters=T)
  outputhap2$ref <- as.character(bases.hap2$value)
  outputhap2$hap2ref <- ifelse(outputhap2$var==outputhap2$ref, 1, 0)
  colnames(outputhap2) <- c("pos", "seq", "refallele", "hap2ref")
  
  
  ######## write table output
  
  
  if(write==TRUE){
    write.table(hap1$ID, file=paste(hap1_name, sep="/"), row.names=F, col.names=F, quote=F)
    write.table(hap2$ID, file=paste(hap2_name, sep="/"), row.names=F, col.names=F, quote=F)
    saveRDS(outputhap1, file = paste(strsplit(hap1_name, split=".txt")[[1]][1], ".rds", sep=""))
    saveRDS(outputhap2, file = paste(strsplit(hap2_name, split=".txt")[[1]][1], ".rds", sep=""))
  }
  
  ####### make and write sif files
  
  haps$ID1 <- ifelse(haps$V2<haps$V4, paste(haps$V2, haps$V1, sep=":"), paste(haps$V4, haps$V3, sep=":"))
  haps$ID2 <- ifelse(haps$V2>haps$V4, paste(haps$V2, haps$V1, sep=":"), paste(haps$V4, haps$V3, sep=":"))
  haps$link <- paste(haps$ID1, haps$ID2)
  sif <- haps[!duplicated(haps$link),5:7]
  
  
  sif$link <- ifelse((sif$ID1%in%hap1$ID&sif$ID2%in%hap1$ID)|(sif$ID1%in%hap2$ID&sif$ID2%in%hap2$ID), "pp",
                     ifelse((sif$ID1%in%hap1$ID&sif$ID2%in%hap2$ID)|(sif$ID1%in%hap2$ID&sif$ID2%in%hap1$ID), "pd", "pr"))
  #Pd=ambiguous
  #Pr=unused
  
  
  sif <- sif[,c(1,2,4,3)]
  
  
  sif$ID1<-ifelse(sif$ID1%in%hap1$ID, paste(sif$ID1, "hap1", sep=":"),
                  ifelse(sif$ID1%in%hap2$ID, paste(sif$ID1, "hap2", sep=":"), 
                         paste(sif$ID1, "nohap", sep=":")))
  sif$ID2<-ifelse(sif$ID2%in%hap1$ID, paste(sif$ID2, "hap1", sep=":"),
                  ifelse(sif$ID2%in%hap2$ID, paste(sif$ID2, "hap2", sep=":"), 
                         paste(sif$ID2, "nohap", sep=":")))
  
  filtered_sif <- subset(sif, sif$link=="pp"|sif$link=="pd"&sif[,1]>5)
  
  filtered_sif <- filtered_sif[,2:4]
  
  write.table(filtered_sif, file=paste(outputname), col.names=F, row.names=F, quote=F)
  
  loopstats
  
}

#families <-c(6,7,8,9,10,11,12)
#for(i in 1:length(families)){
#  family <- families[i]
#  dir <- paste("~/resources/fastqs/HBB/fam", family, "/mapped_TLA/snp-output/", sep="")
#  Dbasename <- paste("HBB-D", family, "_lenient_NS", sep="")
#  Mbasename <- paste("HBB-M", family, "_lenient_NS", sep="")

#haplotyper1.5("~/resources/fastqs/DFSTLA/snp_output/", "DFS_NS_len", "ver.1.5", "chr7")
#  haplotyper1.5(dir, Mbasename, "ver.1.5", "chr11")
#}
directory <- "~/resources/fastqs/CYP300/combined/snp-output/"
basename <- "D6_CYP_all_q5_2515NS"
ext <- ""
chrom <- "chr6"
haplotyper1.5(directory, basename, ext, chrom)

