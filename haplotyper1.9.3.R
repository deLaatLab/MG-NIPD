

write <- T

##basename refers to the name of the output files, these will be written in the current working directory
##haplinks is the .hap file with links found in the input data
## chromosome refers to the chromsome to be analysed, format: "chr1".
haplotyper1.9 <- function(basename, haplinks, chromosome){
  
  badviewpoint <- 0
  viewpoint_override <- 0
if(is.null(write)){
  write <- T
}
  
  if(class(haplinks)!="data.frame"){
    haplinks <- read.table(haplinks)
  }
  
colnames(haplinks) <- c("var1", "pos1", "var2", "pos2", "weight")
#define outputs
outputname <- paste(basename, ".sif", sep=".")
hap1_name <- paste(basename,  ".hap1.txt", sep="" )
hap2_name <- paste(basename,  ".hap2.txt", sep="" )

haps<- haplinks
hapfile <- haplinks

haps <- haps[haps$pos1!=badviewpoint&haps$pos2!=badviewpoint,]
haps$used <- "X"
haps$ID1 <- paste(haps[,2], haps[,1], sep=":")
haps$ID2 <- paste(haps[,4], haps[,3], sep=":")
# haps$ID1 <- ifelse(haps[,2]<haps[,4], paste(haps[,2], haps[,1], sep=":"), paste(haps[,4], haps[,3], sep=":"))
# haps$ID2 <- ifelse(haps[,2]>haps[,4], paste(haps[,2], haps[,1], sep=":"), paste(haps[,4], haps[,3], sep=":"))

#########recalculate strenght and ambiguity with fixed haplotypes


recalculate_all <- function(hap1, hap2){
  
  
  allpositions <- rbind(hap1, hap2)
  allpositions <- allpositions[!duplicated(allpositions$pos),1]
  
  for (i in 1:length(allpositions)){ 
    position <- allpositions[i]
    searchID <- hap1[hap1$pos==position,3]
    
    if(nrow(hap1[hap1$pos==position,])==1){
      links <- subset(hapfile, hapfile$ID1==searchID)
      linkstohap1 <- nrow(links[links$ID2%in%hap1$ID,])
      weight_1 <- sum(links[links$ID2%in%hap1$ID,5])
      linkstohap2 <- nrow(links[links$ID2%in%hap2$ID,])
      ambig_1 <- sum(links[links$ID2%in%hap2$ID,5])
      
      
      plexity_1 <- linkstohap1
      amplexity_1 <- linkstohap2
    }else{
      weight_1 <- 0
      ambig_1 <- 0
      plexity_1 <- 0
      amplexity_1 <- 0}
    
    searchID2 <- as.character(hap2[hap2$pos==position,3])
    
    if(nrow(hap2[hap2$pos==position,])==1){
      links <- subset(hapfile, hapfile$ID1==searchID2)
      linkstohap2 <- nrow(links[links$ID2%in%hap2$ID,])
      weight_2 <- sum(links[links$ID2%in%hap2$ID,5])
      linkstohap1 <- nrow(links[links$ID2%in%hap1$ID,])
      ambig_2 <- sum(links[links$ID2%in%hap1$ID,5])
      plexity_2 <- linkstohap2
      amplexity_2 <- linkstohap1}else{
        weight_2 <- 0
        ambig_2 <- 0
        plexity_2 <- 0
        amplexity_2 <- 0}
    
    weight <- weight_1+weight_2
    ambig <- ambig_1+ambig_2
    plexity <- plexity_1+plexity_2
    amplexity <- amplexity_1+amplexity_2
    
    if(plexity_1!=0){
      hap1[hap1$pos==position,"plexity"] <- plexity
      hap1[hap1$pos==position,"amplexity"] <- amplexity
      hap1[hap1$pos==position,"weight"] <- weight
      hap1[hap1$pos==position,"ambig"] <- ambig}
    
    if(plexity_2!=0){
      hap2[hap2$pos==position,"plexity"] <- plexity
      hap2[hap2$pos==position,"amplexity"] <- amplexity
      hap2[hap2$pos==position,"weight"] <- weight
      hap2[hap2$pos==position,"ambig"] <- ambig}
    
    
  } 
  hap1$hap<-1
  hap2$hap<-2
  
  haplotypes <- rbind(hap1, hap2)
  haplotypes$ratio <- round(haplotypes$ambig/haplotypes$weight,3)
  return(haplotypes)
  
}



### autodefine best viewpoint

poslist <- hapfile[!duplicated(hapfile[,2]),2]

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

# weight_table <- weight_table[1:round(nrow(weight_table)/20,0),]
# middle <- min(weight_table$pos)+((max(weight_table$pos)-min(weight_table$pos))/2)
# 
# viewpointpos<- weight_table[which(abs(weight_table$pos-middle)==min(abs(weight_table$pos-middle))),"pos"]
viewpointpos <- weight_table[1,1]

if(viewpoint_override>1){viewpointpos <- viewpoint_override}
viewpointvars <- subset(haplinks, haplinks$pos1==viewpointpos)
viewpointvars <- unique(viewpointvars$var1)
hap1allele <- as.character(viewpointvars[1])
hap2allele <- as.character(viewpointvars[2])


iterations <- 25
debug_iteration <- 4

##############

hapfile$ID1 <- paste(hapfile[,2], hapfile[,1], sep=":")
hapfile$ID2 <- paste(hapfile[,4], hapfile[,3], sep=":")

IDlist <- unique(c(hapfile[,1], hapfile[,2]))
IDlist <- as.numeric(subset(IDlist, IDlist!=viewpointpos))

hap1 <- data.frame(pos=viewpointpos, var=hap1allele, 
                   ID=paste(viewpointpos, hap1allele, sep=":"), 
                   weight=10000,linkedweight=10000, ambig=0, class="VP", ratio=0, round=0, plexity=0, amplexity=0)
hap1$class <- as.character(hap1$class)

hap2 <- data.frame(pos=viewpointpos, var=hap2allele, 
                   ID=paste(viewpointpos, hap2allele, sep=":"), 
                   weight=10000,linkedweight=10000, ambig=0, class="VP", ratio=0, round=0,plexity=0,amplexity=0)
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
 # while(loopnr <= 21){

  
  badsnps <- data.frame(pos=factor(), reason=character())
  
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
      
      P1_plexity <- nrow(P1_links)
      
      #possibility2
      var1_hap2 <- var1_links[var1_links$ID2%in%seedhap2$ID,]
      var2_hap1 <- var2_links[var2_links$ID2%in%seedhap1$ID,]
      
      P2 <- sum(var1_hap2$weight, var2_hap1$weight)
      
      P2_links <- rbind(var1_hap2, var2_hap1)
      P2_links$linkID <- paste(P2_links$ID1, P2_links$ID2, sep=":")
      
      P2_plexity <- nrow(P2_links)
      var1_line <- data.frame(pos=SNP, var=var1_links[1,1], ID=var1)
      var2_line <- data.frame(pos=SNP, var=var2_links[1,1], ID=var2)
      
      #second if-else statement continues script if links with a haplotype have been found
      if(P1==0&P2==0){
        badline <- data.frame(pos=paste(SNP), reason="isolated")
        badsnps <- rbind(badsnps, badline)}else{
          #if var 1 is in hap 1
          if(P1>P2){
            Weight <- P1
            Ambig <- P2
            if(sum(var1_hap1$weight)>0&sum(var2_hap2$weight)>0) 
            {Class <- "Double_link"} else {Class <- "Single_link" } 
            var1_line$weight=Weight
            var1_line$linkedweight=Weight
            var1_line$ambig=Ambig
            var1_line$class=Class
            var1_line$ratio=Ambig/Weight
            var1_line$round=loopnr
            var1_line$plexity=P1_plexity
            var1_line$amplexity=P2_plexity
            
            var2_line$weight=Weight
            var2_line$linkedweight=Weight
            var2_line$ambig=Ambig
            var2_line$class=Class
            var2_line$ratio=Ambig/Weight
            var2_line$round=loopnr
            var2_line$plexity=P1_plexity
            var2_line$amplexity=P2_plexity
            
            hap1 <- rbind(hap1,var1_line)
            hap2 <- rbind(hap2, var2_line)}
          #if var 1 is in hap2
          if(P2>P1){
            #if var 1 is in hap2
            Weight <- P2
            Ambig <- P1
            if(sum(var1_hap2$weight)>0&sum(var2_hap1$weight)>0) 
            {Class <- "Double_link"} else {Class <- "Single_link" }    
            
            #allsif$used <- ifelse(allsif$linkID%in%P2_links$linkID==TRUE, loopnr, allsif$used)
            
            var1_line$weight=Weight
            var1_line$linkedweight=Weight
            var1_line$ambig=Ambig
            var1_line$class=Class
            var1_line$ratio=Ambig/Weight
            var1_line$round=loopnr
            var1_line$plexity=P2_plexity
            var1_line$amplexity=P1_plexity
            
            var2_line$weight=Weight
            var2_line$linkedweight=Weight
            var2_line$ambig=Ambig
            var2_line$class=Class
            var2_line$ratio=Ambig/Weight
            var2_line$round=loopnr
            var2_line$plexity=P2_plexity
            var2_line$amplexity=P1_plexity
            hap1<- rbind(hap1, var2_line)
            hap2<- rbind(hap2, var1_line)}}}else{         
              
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
  
  
  #first rounds: decreasing quantiles
  #last 3 rounds; all interactions with weight > 1 and no ambiguity
  #last 2 rounds: all interactions
  
  for(i in 1:nrow(hap1)){
  hap1[i,"var"] <- strsplit(as.character(hap1[i,"ID"]), split=":")[[1]][2]}
  for(i in 1:nrow(hap2)){
    hap2[i,"var"] <- strsplit(as.character(hap2[i,"ID"]), split=":")[[1]][2]}
  
  if(loopnr<5){
    #quantileweight <- quantile(hap1[2:length(hap1[,1]),4], (1-(loopnr*0.2)), names=F)}
    quantileweight <- 30-(5*loopnr)
    
    hap1 <- subset(hap1, hap1$ratio<0.1|hap1$class=="VP")
    hap2 <- subset(hap2, hap2$ratio<0.1|hap2$class=="VP")
    #hap1$ambig==0
    
    recalc_haps <- recalculate_all(hap1,hap2)
    hap1 <- recalc_haps[recalc_haps$hap==1,1:(ncol(recalc_haps)-1)]
    hap2 <- recalc_haps[recalc_haps$hap==2,1:(ncol(recalc_haps)-1)]
    
    
    hap1 <- subset(hap1, hap1$linkedweight>quantileweight&hap1$class!="Single_link"|hap1$class=="VP")
    hap2 <- subset(hap2, hap2$linkedweight>quantileweight&hap2$class!="Single_link"|hap2$class=="VP")
    
    
  }
  
  if(loopnr<iterations-4){
    hap1 <- subset(hap1, hap1$ratio<0.3)
    hap2 <- subset(hap2, hap2$ratio<0.3)
    
    recalc_haps <- recalculate_all(hap1,hap2)
    hap1 <- recalc_haps[recalc_haps$hap==1,1:(ncol(recalc_haps)-1)]
    hap2 <- recalc_haps[recalc_haps$hap==2,1:(ncol(recalc_haps)-1)]
    
    hap1 <- subset(hap1, hap1$class!="Single_link")
    hap2 <- subset(hap2, hap2$class!="Single_link")

  }else{
  
  recalc_haps <- recalculate_all(hap1,hap2)
  
  hap1 <- recalc_haps[recalc_haps$hap==1,1:(ncol(recalc_haps)-1)]
  hap2 <- recalc_haps[recalc_haps$hap==2,1:(ncol(recalc_haps)-1)]}
  
  loopstats[loopnr,7] <- length(hap1[,1])
  
  
  
  newhap1 <- hap1[hap1$pos%in%seedhap1$pos==F,]
  newhap2 <- hap1[hap2$pos%in%seedhap2$pos==F,]
  
 
  
  
  #update haps
  seedhap1 <- hap1
  seedhap2 <- hap2
  
  #update IDlist
  IDlist <- IDlist[!IDlist%in%seedhap1$pos]
  IDlist <- IDlist[!IDlist%in%seedhap2$pos]
  
  message(loopnr)
  
  loopnr <- loopnr+1
  #close the while loop
  
}





hap1$hap <- NULL
hap2$hap <- NULL

#start working on bad snps
#remove isolated SNPs
badsnps <- subset(badsnps, badsnps$reason!="isolated")
badsnps$reason <- as.character(badsnps$reason)
singles <- length(badsnps$pos)

hap1_b <- NULL
hap2_b <- NULL

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
  
  P1_plexity <- nrow(hap1_links)
  P2_plexity <- nrow(hap2_links)
  
  P1 <- sum(hap1_links$weight)
  P2 <- sum(hap2_links$weight)
  
  if(P1>0&P2==0){
    hapline <- data.frame(pos=SNP, var=hap_links[1,1], 
                          ID=paste(SNP,hap_links[1,1], sep=":"),
                          weight=P1, linkedweight=P1, ambig=P2, class="Single_linked_bad", 
                          ratio=P2/P1, round=26, plexity=P1_plexity, amplexity=P2_plexity)
    hap1_b<- rbind(hap1_b, hapline)
    badsnps[badsnps$pos==SNP,2] <- "Single_linked"
  }else{
    if(P2>0&P1==0){
      hapline <- data.frame(pos=SNP, var=hap_links[1,1], 
                            ID=paste(SNP,hap_links[1,1], sep=":"),
                            weight=P2,linkedweight=P2, ambig=P1, class="Single_linked_bad", 
                            ratio=P1/P2, round=26, plexity=P2_plexity, amplexity=P1_plexity)
      hap2_b<- rbind(hap2_b, hapline)
      badsnps[badsnps$pos==SNP,2] <- "Single_linked"}
    else{ 
      if(P1>0&P2>0){badsnps[badsnps$pos==SNP,2] <- "Single_ambiguous"}
      else{badsnps[badsnps$pos==SNP,2] <- "Single_isolated"}
    }
  }
  
  #close the for loop
}


### filter bad snps
  hap1_c <- rbind(hap1, hap1_b)
  hap2_c <- rbind(hap2, hap2_b)
  both_c <- recalculate_all(hap1 = hap1_c, hap2=hap2_c)
  both_c <- subset(both_c, both_c$ratio<0.5)
 
  hap1_b <- hap1_b[hap1_b$ID%in%both_c$ID,]
  hap2_b <- hap2_b[hap1_b$ID%in%both_c$ID,]
   

  ##mark the used links
  
  ##mark the used links
 
  hap1 <- rbind(hap1, hap1_b)
  hap2 <- rbind(hap2, hap2_b)
hap1$hap <- 1
hap2$hap <- 2

both <- recalculate_all(hap1 = hap1, hap2=hap2)
hap1 <- both[both$hap==1, 1:(ncol(both)-1)]
hap2 <- both[both$hap==2, 1:(ncol(both)-1)]

##### final fileter step
hap1 <- na.omit(subset(hap1, hap1$ratio<0.5))
hap2 <- na.omit(subset(hap2, hap2$ratio<0.5))

outputhap1 <- hap1[,1:2]
outputhap2 <- hap2[,1:2]



######## write table output


# if(write==TRUE){
#   write.table(hap1$ID, file=paste(hap1_name, sep="/"), row.names=F, col.names=F, quote=F)
#   write.table(hap2$ID, file=paste(hap2_name, sep="/"), row.names=F, col.names=F, quote=F)
#   saveRDS(outputhap1, file = paste(strsplit(hap1_name, split=".txt")[[1]][1], ".rds", sep=""))
#   saveRDS(outputhap2, file = paste(strsplit(hap2_name, split=".txt")[[1]][1], ".rds", sep=""))
# }

####### make and write network files


haps$node1added <- "X"
haps$node2added <- "X"


hap1_links <- NULL
for(i in 2:nrow(hap1)){
links <- NULL
  row <- hap1[i,]
connectables <- as.character(hap1[hap1$round<row$round,"ID"])
links <- haps[haps$ID2==row$ID&haps$ID1%in%connectables,]



if(is.na(links[1,1])==F){
links$used <- row$round
links$node2added <- row$round
for(j in 1:nrow(links)){
  links[j,"node1added"] <- hap1[hap1$ID==links[j,"ID1"],"round"]
  }
 if(exists("hap1_links")==F){
   hap1_links <- links
 }else{hap1_links <- rbind(hap1_links, links)}}else{
  hap1[i,"class"]<- "indirect"
  }
}

hap2_links <- NULL
for(i in 2:nrow(hap2)){
 links <- NULL
   row <- hap2[i,]
  connectables <- as.character(hap2[hap2$round<row$round,"ID"])
  links <- haps[haps$ID1==row$ID&haps$ID2%in%connectables,]
  
  if(is.na(links[1,1])==F){
    links$used <- row$round
    links$node2added <- row$round
    for(j in 1:nrow(links)){
      links[j,"node1added"] <- hap2[hap2$ID==links[j,"ID1"],"round"]
    }
    if(exists("hap2_links")==F){
      hap2_links <- links
    }else{hap2_links <- rbind(hap2_links, links)}}else{
      hap2[i,"class"]<- "indirect"
      }
}


hap1_links$hap <- 1
hap1_links$hap.1 <- 1
hap1_links$interaction <- "pp"
hap2_links$hap <- 2
hap2_links$hap.1 <- 2
hap2_links$interaction <- "pp"

all_links <- rbind(hap1_links, hap2_links)

both <- na.omit(both)

netwtable <- all_links[,c(7,8,13,6,9,10,11,12,2,4,5)]


netwtable$weight1 <- "X"
netwtable$weight2 <- "X"
netwtable$ambig1 <- "X"
netwtable$ambig2 <- "X"
netwtable$plexity1 <- "X"
netwtable$plexity2 <- "X"
netwtable$amplexity1 <- "X"
netwtable$amplexity2 <- "X"

both <- rbind(hap1, hap2)

for(i in 1:nrow(netwtable)){
  row <- netwtable[i,]
  hapdata1 <- both[both$ID==row$ID1,]
  if(nrow(hapdata1)==1){
    netwtable[i,"weight1"] <- hapdata1$weight
    netwtable[i,"ambig1"] <- hapdata1$ambig
    netwtable[i,"plexity1"] <- hapdata1$plexity
    netwtable[i,"amplexity1"] <- hapdata1$amplexity
    
  }
  
  hapdata2 <- both[both$ID==row$ID2,]
  if(nrow(hapdata2)==1){
    netwtable[i,"weight2"] <- hapdata2$weight
    netwtable[i,"ambig2"] <- hapdata2$ambig
    netwtable[i,"plexity2"] <- hapdata2$plexity
    netwtable[i,"amplexity2"] <- hapdata2$amplexity
  }
  
}

netwtable$edgeweight <- netwtable$weight
netwtable$weight <- NULL


colnames(netwtable) <- c("ID1", "ID2", "interaction","interactionused", "added", "added", "hap", "hap",
                         "pos", "pos","weight", "weight", "ambig", "ambig", "plexity", "plexity", "amplexity", "amplexity", 
                         "interactionweight")
hap1$hap <- 1
hap2$hap <- 2
haplotypes <- rbind(hap1, hap2)

list <- strsplit(x = as.character(haplotypes$ID), split=":")
haplotypes$var <- sapply(list,`[`,2)

if(write==T){
  write.table(netwtable, file=paste(basename,  "_network.txt", sep=""),quote=F, row.names = F, sep="\t")
  write.table(haplotypes, file=paste(basename,  ".haps.txt", sep=""),quote=F, row.names = F)}

return(haplotypes)
loopstats

}
# write=T
# haps <- haplotyper1.9(basename="M9_test", haplinks="~/resources/Small files/haps/HBB/HBB-M9_lenient_NS.hap", chromosome="chr11")

haplotyper1.9("family1_D_UMC_test", haplinks="~/resources/fastqs/family1/D1_nextseq/dbsnp_output/D1_dbsnp_Q20.hap", chromosome="chr7")

