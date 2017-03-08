#source("http://bioconductor.org/biocLite.R")
#biocLite("BSgenome.Hsapiens.UCSC.hg19")
#biocLite("GenomicRanges")
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicRanges)
library(beeswarm)
library(ggplot2)
library(gplots)
library(reshape2)



#########changelog#######
#1.4: added chromosome specification
#1.5: -correct use of processed dTLA data to using readcount-filtered pileup.
#     -added changelog
#     -made boxplot more general and more informative
#1.6: added beeswarm
#1.7: made type 3 SNP prediction correctly named 
# and more resistant against mis-calls by taking median instead of total
#1.7.1
#-added beeswarm export
#-added stats matrix
#1.7.2
#-added genotype prediction for type 1 and 3
#Start of the vcf processing towards a concise data table
#1.7.3: edited compensation to compensate with readcount, not allelecount
#1.7.4: added readcount recalculation to cope with Phred-filter
#1.7.5: fixed some bugs
#1.7.6: added sound
#1.9: auto file selection and controlscript added
#1.10: added result table and complete-dataset loop
#1.11: added CYP == T option
#1.12: fixed CYP == T option
#1.13: added dotplot for type 1's
#1.14: added type 3 barplot and fffilter
#1.15 added the successful picture of success
#1.16 added the the noninherited allele to success picture
#1.17 added a guessed inheritance figure

#1.19 Figure 4 table and SD/SE calculations
#1.20: added table outputs for GtM
#1.21: added complete manual input
#1.23: varaan compatibility and fix for the case of identical haps between parents.
#1.25: incorporated RHDO with class a and B snps
# -Added SNP list as output
# -start saving RData files
rm(list=ls())
h <- head

graphs<- T
write <- F
familynr <- 1


update_table <- F
#weirdshittable <- NULL

refbias <- 0.0
ffoverride <- 0.28
pHapoverride <- 0
fffilter <- 0.8
mhapoverride <- 0
#% deviating reads allowed in homozygous SNPs
homcall <- 10

#used only on TLA
minimumcoverage<- 15
# used in cfDNA
mincov_cfDNA <- 15
# used for CVS coverage
mincov_CVS <- 10

famID<- paste(sep="", "_CYP_", familynr)

wd <- paste("/data1/projects/nipd/pipeline/processed_data/family", famID, sep="")
setwd(wd)
root <- wd
fams <- famID
family <- fams
dad <- paste("family", famID, "_D_TLA_CYP", sep="")
mom <- paste("family", famID, "_M_TLA_CYP", sep="")
#cfDNA <- paste("/cfDNA/family", familynr, "_cfDNA_CFTR_Q0.vcf", sep="")
cfDNA <- paste("/cfDNA/family", famID, "_cfDNA_CYP.vcf", sep="")

CVSfile <- paste("~/RHDO/Standard/Family", familynr, "_CYP/CVS",familynr,"_CYP.vcf", sep="")


dadtla <- paste("/haplotypes/", dad, ".vcf", sep="")
momtla <- paste("/haplotypes/", mom, ".vcf", sep="")

mhap1 <- paste("/haplotypes/", mom, "_2515NS.hap1.txt", sep="")
mhap2 <-  paste("/haplotypes/", mom, "_2515NS.hap2.txt", sep="")
phap1 <- paste("/haplotypes/", dad, "_2515NS.hap1.txt", sep="")
phap2 <- paste("/haplotypes/", dad, "_2515NS.hap2.txt", sep="")


#################### start

hom_high <- 100-homcall
hom_low <- homcall

mhap1 <- paste(wd, mhap1, sep="/")
mhap2 <- paste(wd, mhap2, sep="/")
phap1 <- paste(wd, phap1, sep="/")
phap2 <- paste(wd, phap2, sep="/")

mhap1_file <- mhap1
mhap2_file <- mhap2
phap1_file <- phap1
phap2_file <- phap2



savebeeswarm <- FALSE


##################################Part 1######################

vcf.processing <- function(filename){
  input <- data.frame(read.table(paste(wd,filename,sep="/")))
  
  info <- input[,8]
  info <- gsub(pattern = "DP=", replacement = "0,", x=info) #get rid of the = signs
  info <- gsub(pattern = ";I16=", replacement = ",0,", x=info)
  infotest <- strsplit(info, split=",")
  
  parseInfo <- function(n,x)
  {return(as.numeric(strsplit(x[n],split=",")[[1]][c(2,4:7)]))}
  
  allInfo <- t(sapply(1:length(info),parseInfo,x=info))
  info <- as.data.frame(allInfo)
  
  chr <- input[,1]
  pos <- input[,2]
  refallele <- input[,4]
  altallele <- input[,5]
  readcount <- info[,1]
  reftot <- info[,2]+info[,3]
  alttot <- info[,4]+info[,5]
  
  Important.info <- data.frame(chr,pos,refallele,altallele,readcount=reftot+alttot,reftot,alttot)
  return(Important.info)
}

mtla <- vcf.processing(filename=momtla)
chromosome <- as.character(mtla[1,1])


dtla <- vcf.processing(filename=dadtla)

Important.info.cfdna <- vcf.processing(filename=cfDNA)
Important.info.cfdna$readcount <- Important.info.cfdna$reftot+Important.info.cfdna$alttot

#filter cfDNA for readcount
Important.info.cfdna <- subset(Important.info.cfdna, Important.info.cfdna[,5]>mincov_cfDNA)
#compensate cfDNA for reference bias
Important.info.cfdna[,6] <- ifelse(Important.info.cfdna$reftot>0, 
                                   Important.info.cfdna$reftot-(refbias*Important.info.cfdna$readcount), 
                                   Important.info.cfdna$reftot) 

Important.info.cfdna$alttot <- ifelse(Important.info.cfdna$alttot>0, 
                                      Important.info.cfdna$alttot+(refbias*Important.info.cfdna$readcount), 
                                      Important.info.cfdna$alttot)
Important.info.cfdna$readcount <- Important.info.cfdna$reftot+Important.info.cfdna$alttot
#Important.info.cfdna <- subset(Important.info.cfdna, Important.info.cfdna[,5]>20)

#all VCFs are now processed to contain useful information


################################################## Part 2  ################

#make SNP calls from pileup
fathercounts <- dtla
mothercounts <- mtla


#filter reads with low coverage
fathercounts <- subset(fathercounts, fathercounts[,5]>minimumcoverage)
mothercounts <- subset(mothercounts, mothercounts[,5]>minimumcoverage)
#filter for reads that are present in both datasets
common.coverage <- intersect(fathercounts[,2],mothercounts[,2])
fathercounts <- fathercounts[fathercounts[,2]%in%common.coverage,]
mothercounts <- mothercounts[mothercounts[,2]%in%common.coverage,]
##test <- data.frame(mothercounts[,2],fathercounts[,2])
#calculate the %ref statistic to call  homref,homalt or het
# %ref = (ref-alt)/readcount*100; round() function limits the number of decimals
fathercounts[,8]<- round(100*(fathercounts$reftot/fathercounts$readcount),1)
colnames(fathercounts)[8]<- "P%ref"
mothercounts[,8]<- round(100*(mothercounts$reftot/mothercounts$readcount),1)
colnames(mothercounts)[8]<- "M%ref"
#make a list with unified SNP positions and %ref of both parents
unified.list <- data.frame(mothercounts[,2],mothercounts[,8],fathercounts[,8])
colnames(unified.list) <- c("pos","M%ref","P%ref")
#now we can identify SNP types from this list; -100~-90 is homalt, 90-100 is homref
#-80 ~ 90 is het
#Type 2 SNPs are equal between parents, and thus, useless
#type 1 SNPs are SNPs where father and mother are homozygous for different alleles
#type 1.1: mother ref, father alt. 
#type 1.2: mother alt, father ref
type1.1 <- subset(unified.list, unified.list[,2]>hom_high & unified.list[,3]< hom_low)
type1.2 <- subset(unified.list, unified.list[,2]< hom_low & unified.list[,3]>hom_high)
#type 3 SNPs are SNPs where father is het and mother is hom.
#called later


############################## Part 3 ##############
#from all SNP counts in cfDNA, gather type 1.1 SNPs and 1.2 SNPs (no haplotypes needed)
covered1.1 <- intersect(type1.1[,1],Important.info.cfdna[,2])
covered1.2 <- intersect(type1.2[,1], Important.info.cfdna[,2])
#pull 1.1 positions from cfDNA and calculate ff
data1.1 <- match(covered1.1,Important.info.cfdna[,2])
data1.1 <- Important.info.cfdna[data1.1,]
data1.1[,8]<- round(data1.1[,7]/data1.1[,5],2)
colnames(data1.1)[8]<- "ff"
weirdsnps <- data1.1[data1.1$ff>0.5*fffilter,2]
data1.1 <- subset(data1.1, data1.1$ff<(0.5*fffilter))
#ff<-  round(2*sum(data1.1[,7])/sum(data1.1[,5]),2)

#do the same for 1.2 positions
data1.2 <- match(covered1.2,Important.info.cfdna[,2])
data1.2 <- Important.info.cfdna[data1.2,]
data1.2[,8]<- round(data1.2[,6]/data1.2[,5],2)
colnames(data1.2)[8]<- "ff"
weirdsnps <- c(weirdsnps, data1.2[data1.2$ff>0.5*fffilter,2])
data1.2 <- subset(data1.2, data1.2$ff<(0.5*fffilter))

fftype1 <- c(data1.1$ff, data1.2$ff)
fftype1reads <- sum(c(data1.1$readcount, data1.2$readcount))

#make fetal haplotype predictions
type_1 <- rbind(
  data.frame(chr=data1.1$chr, pos= data1.1$pos, refallele=data1.1$refallele, 
             Mref.inh= rep(1, each=length(data1.1$pos)), 
             Prefinh = rep(0,each=length(data1.1$pos)), reftot=data1.1$reftot, alttot=data1.1$alttot),
  data.frame(chr=data1.2$chr, pos= data1.2$pos, refallele=data1.2$refallele, 
             Mref.inh=rep(0, each=length(data1.2$pos)), 
             Prefinh=rep(1, each=length(data1.2$pos)),reftot=data1.2$reftot, alttot=data1.2$alttot)
)

write.table(type_1, file=paste(wd,"/family", famID,"_type1.txt", sep=""), quote=F, row.names=F)

######################################### Part 4#####################

fillin <- function(hap1,hap2){
  
  hap1.table <- read.table(hap1,sep=":")
  hap2.table <- read.table(hap2,sep=":")
  hap1.frame<- data.frame(hap1.table)
  hap1.frame[,1] <- as.numeric(hap1.frame[,1])
  hap1.frame <- hap1.frame[order(hap1.frame[,1]),]
  hap1.frame<- hap1.frame[!duplicated(hap1.frame[,1]),]
  hap2.frame<- data.frame(hap2.table)
  hap2.frame[,1] <- as.numeric(hap2.frame[,1])
  hap2.frame <- hap2.frame[order(hap2.frame[,1]),]
  hap2.frame <- hap2.frame[!duplicated(hap2.frame[,1]),]
  #find those SNPs that are found in both haplotypes and list with hap1 variant
  matched <- hap1.frame[hap1.frame[,1]%in%hap2.frame[,1],]
  ranges.matched <- GRanges(seqnames=Rle(chromosome), ranges=matched[,1])
  bases.matched <- data.frame(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, ranges.matched),stringsAscharacters=T)
  matched[,3]<- as.character(bases.matched[,3])
  matched[,2]<- as.character(matched[,2])
  matched[,4] <- ifelse(matched[,2] == matched[,3],1,0)
  matched<- data.frame(matched[,1],matched[,3:4])
  
  #Create a similar list for those SNPs that are found in only one of the two haplotypes
  #returns a list of positions that are haplotyped in hap1 and not in hap2
  unmatched12 <- setdiff(hap1.frame[,1],hap2.frame[,1]) 
  unmatched12 <- hap1.frame[hap1.frame[,1]%in%unmatched12,]
  if(length(unmatched12[,1])>0){
    ranges12 <- GRanges(seqnames=Rle(chromosome), ranges=unmatched12[,1])
    refbases12 <- data.frame(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, ranges12),stringsAscharacters=T)
    unmatched12[,3]<- as.character(refbases12[,3])
    unmatched12[,2]<- as.character(unmatched12[,2])
    unmatched12[,4] <- ifelse(unmatched12[,2] == unmatched12[,3],1,0)}
  #output an extra column where 1 means haplotype 1 is equal to the reference
  #do the same for haplotype 2, but output 0 if haplotype 2 is equal to the reference and 1 if it is not. 
  
  unmatched21 <- setdiff(hap2.frame[,1],hap1.frame[,1]) 
  unmatched21 <- hap2.frame[hap2.frame[,1]%in%unmatched21,]
  if(length(unmatched21[,1])>0){
    ranges21 <- GRanges(seqnames=Rle(chromosome), ranges=unmatched21[,1])
    refbases21 <- data.frame(getSeq(x=BSgenome.Hsapiens.UCSC.hg19, ranges21),stringsAscharacters=T)
    unmatched21[,3]<- as.character(refbases21[,3])
    unmatched21[,2]<- as.character(unmatched21[,2])
    unmatched21[,4] <- ifelse(unmatched21[,2]== unmatched21[,3],0,1)}
  
  if(length(unmatched21[,1])>0|length(unmatched12[,1])>0){
    allunmatched <- rbind(unmatched12,unmatched21)
    allunmatched <- data.frame(allunmatched[,1],allunmatched[,3:4])
    colnames(allunmatched)[1]<- "pos"
    colnames(matched)[1] <- "pos"
    unifiedhap1<- rbind(matched,allunmatched)}else{
      colnames(matched)[1] <- "pos"
      unifiedhap1 <- matched
    }
  
  unifiedhap1 <- unifiedhap1[order(unifiedhap1[,1]),]
  colnames(unifiedhap1)[2]<- "ref"
  colnames(unifiedhap1)[3]<- "hap1ref"
  return (unifiedhap1)
}

unifiedPhap1 <- fillin(hap1= phap1, hap2= phap2)
unifiedPhap2 <- fillin(hap1= phap2, hap2= phap1)
colnames(unifiedPhap2) <- c("pos", "ref", "hap2ref")

phap1_raw <- read.table(phap1_file, sep=":")
phap2_raw <- read.table(phap2_file, sep=":")
phap1_raw$V2 <- as.character(phap1_raw$V2)
phap2_raw$V2 <- as.character(phap2_raw$V2)

completePhap2 <- unifiedPhap2[,c(1,3)]
phap2_homs <- data.frame(pos=fathercounts$pos, Pref=fathercounts[,8])
## 1 is homref, 0 is homalt
phap2_homs$hap2ref <- ifelse(phap2_homs$Pref>hom_high, 1, 
                             ifelse(phap2_homs$Pref<hom_low, 0, NA))  
phap2_homs <- subset(phap2_homs, phap2_homs$hap2ref==1|phap2_homs$hap2ref==0)
phap2_homs <- phap2_homs[,c(1,3)]
completePhap2 <- rbind(completePhap2, phap2_homs)
completePhap2 <- completePhap2[order(completePhap2$pos),]

completePhap1 <- unifiedPhap1[,c(1,3)]
phap1_homs <- phap2_homs
colnames(phap1_homs) <- c("pos", "hap1ref")
completePhap1 <- rbind(completePhap1, phap1_homs)
completePhap1 <- completePhap1[order(completePhap1$pos),]

unifiedMhap1 <- fillin(hap1= mhap1, hap2=mhap2)
unifiedMhap2 <- fillin(hap1= mhap2, hap2= mhap1)
colnames(unifiedMhap2) <- c("pos", "ref", "hap2ref")
completeMhap2 <- unifiedMhap2[,c(1,3)]
mhap2_homs <- data.frame(pos=mothercounts$pos, Pref=mothercounts[,8])
mhap2_homs$hap2ref <- ifelse(mhap2_homs$Pref>hom_high, 1, 
                             ifelse(mhap2_homs$Pref<hom_low, 0, NA)) 
mhap2_homs <- subset(mhap2_homs, mhap2_homs$hap2ref==1|mhap2_homs$hap2ref==0)
mhap2_homs <- mhap2_homs[,c(1,3)]

completeMhap2 <- rbind(completeMhap2, mhap2_homs)
completeMhap2 <- completeMhap2[order(completeMhap2$pos),]

completeMhap1 <- unifiedMhap1[,c(1,3)]
mhap1_homs <- mhap2_homs
colnames(mhap1_homs) <- c("pos", "hap1ref")
completeMhap1 <- rbind(completeMhap1, mhap1_homs)
completeMhap1 <- completeMhap1[order(completeMhap1$pos),]

####################################Part5#######################

#extract data for type 3 SNPs

#paternal haplotype determination
#extract SNPs where hap 1 is reference
hap1ref <- subset(unifiedPhap1, unifiedPhap1$hap1ref==1)
hap1alt <- subset(unifiedPhap1, unifiedPhap1$hap1ref==0)
#for hap1ref, select SNPs where mother is hom alt
Momhomref <- subset(unified.list, unified.list[,2]>hom_high)
Momhomalt <- subset(unified.list, unified.list[,2]<hom_low)

hap1alt$momref <- Momhomref[match(hap1alt$pos,Momhomref$pos),2]
type3.1 <- na.omit(hap1alt)
hap1ref$momref <- Momhomalt[match(hap1ref$pos,Momhomalt$pos),2]
type3.2 <- na.omit(hap1ref)

type3.1[,5:7]<-Important.info.cfdna[match(type3.1$pos, Important.info.cfdna$pos),5:7]
type3.1 <- na.omit(type3.1)
type3.2[,5:7]<-Important.info.cfdna[match(type3.2$pos, Important.info.cfdna$pos),5:7]
type3.2 <- na.omit(type3.2)
type3.1$obs <- round(type3.1$alttot/type3.1$readcount, 3)
type3.2$obs <- round(type3.2$reftot/type3.2$readcount, 3)
weirdsnps <- c(weirdsnps, type3.1[type3.1$obs>(0.5*fffilter),1], type3.2[type3.2$obs<(0.5*fffilter),1])
type3.1 <- subset(type3.1, type3.1$obs<(0.5*fffilter))
type3.2 <- subset(type3.2, type3.2$obs<(0.5*fffilter))

Phap1reads <- sum(c(type3.1$readcount, type3.2$readcount))

#if hap1 is ref, hap 2 is alt
hap2ref <- subset(unifiedPhap1, unifiedPhap1$hap1ref==0)
hap2alt <- subset(unifiedPhap1, unifiedPhap1$hap1ref==1)

hap2ref$momref <- Momhomalt[match(hap2ref$pos,Momhomalt$pos),2]
type3.2b <- na.omit(hap2ref)
hap2alt$momref <- Momhomref[match(hap2alt$pos,Momhomref$pos),2]
type3.1b <- na.omit(hap2alt)

type3.2b[,5:7]<-Important.info.cfdna[match(type3.2b$pos, Important.info.cfdna$pos),5:7]
type3.2b <- na.omit(type3.2b)

type3.1b[,5:7]<-Important.info.cfdna[match(type3.1b$pos, Important.info.cfdna$pos),5:7]
type3.1b <- na.omit(type3.1b)

type3.1b$obs <- round(type3.1b$alttot/type3.1b$readcount, 3)
type3.2b$obs <- round(type3.2b$reftot/type3.2b$readcount, 3)
type3.1b <- subset(type3.1b, type3.1b$obs<(0.5*fffilter))
type3.2b <- subset(type3.2b, type3.2b$obs<(0.5*fffilter))


type3processing <- function(SNPgroup){
  empty <- nrow(SNPgroup)
  if(empty!=0){
    SNPgroup$hap2allele <- as.character(dtla[dtla$pos%in%SNPgroup$pos,4])
    SNPgroup$cfDNA_alt <- as.character(Important.info.cfdna[Important.info.cfdna$pos%in%SNPgroup$pos,4])
    strsplit(SNPgroup$cfDNA_alt, split=",")[[1]][1]
    for(i in 1:nrow(SNPgroup)){
      row <- SNPgroup[i,]
      hap2allele <- strsplit(row$hap2allele,split=",")[[1]][1]
      SNPgroup[i,9]<- hap2allele
      cfDNA_allele <- strsplit(row$cfDNA_alt, split=",")[[1]][1]
      SNPgroup[i,10]<- cfDNA_allele
      if(hap2allele!=cfDNA_allele){
        SNPgroup[i,8]<- 0
        SNPgroup[i,7]<- 0
      }
    }}
  return(SNPgroup)
}
type3.1 <- type3processing(type3.1)
type3.1b <- type3processing(type3.1b)

Phap2reads <- sum(c(type3.1b$readcount, type3.2b$readcount))

Phap2obs <- (sum(type3.1b[type3.1b$obs!=0,7])+sum(type3.2b[type3.2b$obs!=0,6]))/Phap2reads
Phap1obs <- (sum(type3.1[,7])+sum(type3.2[,6]))/Phap1reads
Phap2cnt <- (sum(type3.1b[type3.1b$obs!=0,7])+sum(type3.2b[type3.2b$obs!=0,6]))
Phap1cnt <- (sum(type3.1[,7])+sum(type3.2[,6]))


if(Phap1obs>Phap2obs){phap1inherited=TRUE} else {phap1inherited=FALSE}


if(is.na(Phap2obs)){
  if(Phap1obs>0.01){phap1inherited=TRUE}
  if(Phap1obs<0.01){phap1inherited=FALSE}
}
if(is.na(Phap1obs)){
  if(Phap2obs>0.01){phap1inherited=FALSE}
  if(Phap2obs<0.01){phap1inherited=TRUE}
}

if(pHapoverride==1)
{phap1inherited=TRUE}
if(pHapoverride==2)
{phap1inherited=FALSE}



#####predict fetal genotype for type 3 SNPs:
#if phap1 is inherited, type 3a SNPs are Het in child and type 3b are Hom
#if Phap1 is not inherited, type 3b snps are Het and type 3a are Hom
#create genotype predictions as if haplotypes are inherited

type_3a <- rbind(
  data.frame(chr=rep(chromosome, each=length(type3.1$pos)), pos= type3.1$pos, refallele=type3.1$ref, 
             Mref.inh=rep(1,each=length(type3.1$pos)), 
             Prefinh=rep(0, each=length(type3.1$pos))),
  data.frame(chr=rep(chromosome, each=length(type3.2$pos)), pos= type3.2$pos, refallele=type3.2$ref, 
             Mref.inh=rep(0,each=length(type3.2$pos)), 
             Prefinh=rep(1, each=length(type3.2$pos)))
)

type_3b <- rbind(
  data.frame(chr=rep(chromosome, each=length(type3.1b$pos)), pos= type3.1b$pos, refallele=type3.1b$ref, 
             Mref.inh=rep(1,each=length(type3.1b$pos)), 
             Prefinh=rep(0,each=length(type3.1b$pos))),
  data.frame(chr=rep(chromosome, each=length(type3.2b$pos)), pos= type3.2b$pos, refallele=type3.2b$ref, 
             Mref.inh=rep(0,each=length(type3.2b$pos)), 
             Prefinh=rep(1,each=length(type3.2b$pos)))
)

#if phap 1 is inherited, phap 2 is not-> invert Prefinh

if(phap1inherited==TRUE){ type_3b$Prefinh <- ifelse(type_3b$Prefinh==0, 1, 0)}else{type_3a$Prefinh <- ifelse(type_3a$Prefinh==0, 1, 0)}

type_3 <- rbind(type_3a, type_3b)

write.table(type_3, file=paste(wd,"/family", famID,"_type2.txt", sep=""), quote=F, row.names=F)

unifiedPhap1[,4]<- phap1inherited
unifiedPhap1[,5]<- ifelse(unifiedPhap1[,4]==TRUE & unifiedPhap1[,3]==0, 0,
                          ifelse(unifiedPhap1[,4]==TRUE & unifiedPhap1[,3]==1, 1,
                                 ifelse(unifiedPhap1[,4]==FALSE & unifiedPhap1[,3]==0, 1,
                                        ifelse(unifiedPhap1[,4]==FALSE & unifiedPhap1[,3]==1, 0,NA))))
inheritedPhap<- data.frame(unifiedPhap1[,1:2],unifiedPhap1[,5])
colnames(inheritedPhap)[3]<- "inh.ref"
#if phap1 is inherited; ref is inherited where phapref=1
#if phap1 is not inherited, ref is inherited where phapref=0

#additional ff determinations
type3a<- c(type3.1$obs, type3.2$obs)
type3b <- c(type3.1b$obs, type3.2b$obs)

if(phap1inherited==TRUE){fftype3<-type3a}else{fftype3<-type3b}

fftype1and3 <- c(fftype3,fftype1)

ff <- 2*(median(fftype1and3))
ffcount <- if(phap1inherited==TRUE){fftype1reads+Phap1reads}else{fftype1reads+Phap2reads}
if(ffoverride>0){
  ff <- ffoverride}


#############################Part6################
#### divide class 3 (maternal het) SNPs between class a and b
#unifiedMhap contains all heterozygous SNPs and only heterozygous snps
#class 3a: Inherited Phap is identical to mHap1
#class 3b: Inherited Phap is identical to mHap
#inherited Phap contains all het SNPs with known paternal inheritance, but not the paternal hom SNPs

unifiedMhap1$class <- 0
for( i in 1:nrow(unifiedMhap1)){
  pvar <- 2
  snp <- unifiedMhap1[i,1]
  paternalsnp <- inheritedPhap[inheritedPhap$pos==snp,]
  if(is.na(paternalsnp[1,1])){
    paternalsnp <- phap1_homs[phap1_homs$pos==snp,]
    if(is.na(paternalsnp[1,1])==FALSE){pvar <- paternalsnp[1,2]}
  }else{pvar <- paternalsnp[1,3]}
  
  if(pvar==0|pvar==1){
    if(unifiedMhap1[i,3]==pvar){unifiedMhap1[i,4]="A"}
    if(unifiedMhap1[i,3]!=pvar){unifiedMhap1[i,4]="B"}
  } 
}

unifiedMhap1 <- unifiedMhap1[unifiedMhap1$class!=0,]
RHDO_snps <- subset(Important.info.cfdna, Important.info.cfdna$pos%in%unifiedMhap1$pos)
coveredMhap1 <- subset(unifiedMhap1, unifiedMhap1$pos%in%RHDO_snps$pos)
RHDO_snps$class <- coveredMhap1$class
RHDO_snps$mhap1ref <- coveredMhap1$hap1ref
RHDO_snps$hap1_percentage <- 0
for(i in 1:nrow(RHDO_snps)){
  hap1ref <- RHDO_snps[i,9]
  hap1reads <- ifelse(hap1ref==1, RHDO_snps[i,6], RHDO_snps[i,7])
  percentage <- (hap1reads/RHDO_snps[i,5])*100
  RHDO_snps[i,10] <- percentage
}

weirdsnps <- c(weirdsnps, RHDO_snps[RHDO_snps$hap1_percentage>80|RHDO_snps$hap1_percentage<20,2])
## but which one is inherited?
## get netto counts per hap
hap1counts <-sum(RHDO_snps[RHDO_snps$mhap1ref==1,6])+sum(RHDO_snps[RHDO_snps$mhap1ref==0,7])
hap2counts <-sum(RHDO_snps[RHDO_snps$mhap1ref==1,7])+sum(RHDO_snps[RHDO_snps$mhap1ref==0,6])

RHDO_snps<- subset(RHDO_snps, RHDO_snps$hap1_percentage!=100&RHDO_snps$hap1_percentage!=0)
class_A <- RHDO_snps[RHDO_snps$class=="A",]
class_B <- RHDO_snps[RHDO_snps$class=="B",]
A_perc_med <- median(class_A$hap1_percentage)
B_perc_med <- median(class_B$hap1_percentage)

total_A <- (sum(class_A[class_A$mhap1ref==1,6])+sum(class_A[class_A$mhap1ref==0,7]))/sum(class_A$readcount)
total_B <- (sum(class_B[class_B$mhap1ref==1,6])+sum(class_B[class_B$mhap1ref==0,7]))/sum(class_B$readcount) 

if((nrow(class_A)<20)|is.na(class_A[1,1])==TRUE){noclassA<-TRUE}else{noclassA<-FALSE}
if((nrow(class_B)<20)|is.na(class_B[1,1])==TRUE){noclassB<-TRUE}else{noclassB<-FALSE}
#if both SNP classes are well represented
if(nrow(class_A)>19&nrow(class_B)>19){
  mhap_undeterminable <- FALSE
  
  ## if haplotype 1 is inherited, type A SNPs are overrepresented and type B snps are not
  ## overrepresentation should correspond (somewhat) to the ff.
  ## class A overrepresented and class B equal
  #closest to 50 is considered to be equal
  A_dev_50 <- A_perc_med-50
  B_dev_50 <- B_perc_med-50
  A_pos_dev_50 <- ifelse(A_dev_50>0, A_dev_50, A_dev_50*(-1))
  B_pos_dev_50 <- ifelse(B_dev_50>0, B_dev_50, B_dev_50*(-1))
  total_A <- (sum(class_A[class_A$mhap1ref==1,6])+sum(class_A[class_A$mhap1ref==0,7]))/sum(class_A$readcount)
  total_B <- (sum(class_B[class_B$mhap1ref==1,6])+sum(class_B[class_B$mhap1ref==0,7]))/sum(class_B$readcount)    
  ##class A equal, class B underrepresented
  if(A_pos_dev_50<B_pos_dev_50&B_dev_50<0){mhap1inherited <- FALSE}
  ##class B equal, class A underrepresented
  if(B_pos_dev_50<A_pos_dev_50&A_dev_50>0){mhap1inherited <- TRUE}
  
}


## if either snp class is underrepresented as with consanguineous parents, look only at the well represented ones. 
if(noclassA==TRUE&noclassB==FALSE){
  mhap_undeterminable <- TRUE
  ## no class A reads: many class B reads; the shared haplotype is mHap2
  ## We expect underrepresentation of hap 1 in class B SNPs if the shared haplotype is inherited
  ## Use knowledge of the fetal fraction to appreciate the meaning
  expect <- (100*ff)
  ## however we have to take variation into account
  expect <- -1*(expect/2)
  #calculate representation, positive means hap1 overrepresentation
  diff <- B_perc_med-50
  
  #if the difference if smaller than expected, hap 1 is probably not represented
  if(diff<expect|diff<(-3)){mhap1inherited<-FALSE}
  #if the difference is bigger than expected, hap1 is probably represented
  if(diff>expect){mhap1inherited<-TRUE}
}
if(noclassB==TRUE&noclassA==FALSE){
  mhap_undeterminable <- TRUE
  ## no class B reads: many class A reads; the shared haplotype is mHap1
  ## We expect overrepresentation of hap 1 in class B SNPs if the shared haplotype is inherited
  ## Use knowledge of the fetal fraction to appreciate the meaning
  ## A_perc_med is the median representation of hap 1 in class A snps, it should be 50+ff
  expect <- (100*ff)
  ## however we have to take variation into account
  expect <- expect/2
  #calculate representation, positive means hap1 overrepresentation
  diff<- A_perc_med-50
  #if the difference if smaller than expected, hap 1 is probably not represented
  if(diff<expect){mhap1inherited<-FALSE}
  #if the difference is bigger than expected, hap1 is probably represented
  if(diff>expect|diff>3){mhap1inherited<-TRUE}
}

if(mhapoverride==1){
  mhap1inherited <- TRUE} 
if(mhapoverride==2){
  mhap1inherited <- FALSE} 
plotdata <- data.frame(perc =RHDO_snps$hap1_percentage, class= RHDO_snps$class)
#beeswarm

if(savebeeswarm == TRUE){
  pdf(file=paste(wd, "beeswarm.pdf", sep="/"))
  
  beeswarm(perc~class,data=plotdata, 
           method ="center", 
           labels=c(paste("class A n=", nrow(plotdata[plotdata$class=="A",]), sep=""),
                    paste("class B n=", nrow(plotdata[plotdata$class=="B",]), sep="")),
           main ='Overrepresentation in cfDNA',
           ylab="RHDO (%)")
  #legend("top",legend=c("het","hom"), title="", pch=16, col=c(2,1))
  #legend("top",legend=c("",""), title="Ch. genotype", cex=0.7, bty='n')
  abline(h=50)
  dev.off()}

beeswarm(perc~class,data=plotdata, 
         method ="center", 
         labels=c(paste("class A n=", nrow(plotdata[plotdata$class=="A",]), sep=""),
                  paste("class B n=", nrow(plotdata[plotdata$class=="B",]), sep="")),
         main ='Overrepresentation in cfDNA',
         ylab="RHDO (%)")
#legend("top",legend=c("het","hom"), title="", pch=16, col=c(2,1))
#legend("top",legend=c("",""), title="Ch. genotype", cex=0.7, bty='n')
abline(h=50)


####### get total number of reads making up RHDO analysis


coveragecount <- median(RHDO_snps$readcount)

snptypedata <- c(
  sum(length(data1.1[,1]), length(data1.2[,1])),
  sum(length(type3.1[,1]), length(type3.2[,1]),length(type3.1b[,1]), length(type3.2b[,1])),
  sum(length(type3.1[,1]), length(type3.2[,1])),
  sum(length(type3.1b[,1]), length(type3.2b[,1])),
  nrow(RHDO_snps),
  coveragecount
)
stats <- data.frame(type= c("type 1", "type 2", "type 2a", "type 2b", "RHDO snps", "median coverage"), 
                    N=snptypedata)

data1.1$class <- "1a"
data1.2$class <- "1b"
SNPlist <- data.frame(chr=character(), pos=numeric(), class=character())
SNPlist <- rbind(SNPlist, rbind(data1.1[,c(1,2,9)], data1.2[,c(1,2,9)]))

type3.1$class <- "2a"
type3.1b$class <- "2b"
type3.2$class <- "2a"
type3.2b$class <- "2b"
chr <- as.character(mothercounts[1,1])
alltype2s <- rbind(rbind(type3.1[,c(1,2,3,4,5,6,7,8,11)], type3.1b[,c(1,2,3,4,5,6,7,8,11)]), rbind(type3.2, type3.2b))
alltype2s$chr <- chr

SNPlist <- rbind(SNPlist, alltype2s[,c(10,1,9)])
RHDO_A <- class_A
RHDO_B <- class_B
RHDO_A$class <- "3a"
RHDO_B$class <- "3b"

RHDOlist <- rbind(rbind(RHDO_A[,c(1,2,9,8,5,6,7,10)], RHDO_B[,c(1,2,9,8,5,6,7,10)]))
write.table(RHDOlist, file=paste(wd, "/family",famID,"_RHDO.txt", sep=""), quote=F, row.names=F)

### collect outlierSNPs
if(update_table==TRUE){
  if(exists("weirdshittable")==F){
    weirdshittable <- data.frame(pos=weirdsnps, fam=familynr, obs=1)}
  if(exists("weirdshittable")==T){
    weirdtable <- data.frame(pos=weirdsnps, fam=familynr, obs=1)
    weirdshittable <- rbind(weirdshittable, weirdtable)
    for(i in 1:nrow(weirdshittable)){
      pos <- weirdshittable[i,1]
      obs <- nrow(weirdshittable[weirdshittable$pos==pos,])
      weirdshittable[i,3] <- obs
    }
    
  }}

pdf(file=paste(sep="","family", fams, "_", "patinh.pdf"))
barplot(height = c((Phap1cnt/Phap1reads)*100, (Phap2cnt/Phap2reads)*100), col=c(rgb(0.6,0,0.6,0.5), rgb(0,0.8,0,0.5)),
        xlab=paste("phap1: ", Phap1cnt, "/", Phap1reads," n=", nrow(type3.1)+nrow(type3.2), " phap2: ",Phap2cnt, "/", Phap2reads,
                   " n=",nrow(type3.1b)+nrow(type3.2b),  sep=""),
        ylim=c(0,max(c((Phap1cnt/Phap1reads)*100, (Phap2cnt/Phap2reads)*100))+3),
        ylab= "paternal observations in cfDNA (%)")
dev.off()

##############################output##################
answer <- ifelse(phap1inherited==TRUE, 1, 2)
distance <- RHDO_snps[length(RHDO_snps[,2]),2]-RHDO_snps[1,2]



a <- paste("fetal fraction determination:")
b <-paste(ffcount, "reads used to establish a median ff of", ff)

c<- paste("paternal inheritance determination:")
d<- paste ("paternal haplotype 1 fraction:", round(Phap1obs, 3) , "determined from", Phap1reads, "reads")
e<-paste ("paternal haplotype 2 fraction:", round(Phap2obs,3) , "determined from", Phap2reads, "reads")
f<-paste("so we assumed the paternally inherited haplotype is number", answer)
g<-paste("we observed maternal haplotype 1", hap1counts,"times, and maternal haplotype 2", hap2counts, "times.")
H<-paste("over a distance of", round(distance/1000,1), "kb")

if(mhap_undeterminable==FALSE){
  k<-paste("class A snps represent haplotype 1 ", round(A_perc_med,2), "% and class B snps represent haplotype 1 ", round(B_perc_med,2), "%", sep="")
  L<- paste("so most likely inherited maternal haplotype is hap ", ifelse(mhap1inherited==TRUE, 1, 2), sep="")
}



if(mhap_undeterminable==TRUE&noclassA==TRUE){
  k<-paste("there are no class A SNPs and class B snps represent haplotype 1 ", round(B_perc_med,2), "%", sep="")
  L <- paste("the best guess so far would be haplotype ", if(mhap1inherited){1}else{2})
}

if(mhap_undeterminable==TRUE&noclassB==TRUE){
  k<-paste("there are no class B SNPs and class A snps represent haplotype 1 ", round(A_perc_med,2), "%", sep="")
  L <- paste("the best guess so far would be haplotype ", if(mhap1inherited){1}else{2})
}

if(mhap1inherited==TRUE){inheritedMhap <- completeMhap1}else{inheritedMhap <- completeMhap2}
if(phap1inherited==TRUE){inheritedPhap <- completePhap1}else{inheritedPhap <- completePhap2}
if(mhap1inherited==TRUE){noninheritedMhap <- completeMhap2}else{noninheritedMhap <- completeMhap1}
if(phap1inherited==TRUE){noninheritedPhap <- completePhap2}else{noninheritedPhap <- completePhap1}

inheritedPhap <- inheritedPhap[inheritedPhap$pos%in%inheritedMhap$pos,]
colnames(inheritedPhap) <- c("pos", "Pgeno")
inheritedPhap <- inheritedPhap[duplicated(inheritedPhap$pos)==FALSE,]
inheritedMhap <- inheritedMhap[inheritedMhap$pos%in%inheritedPhap$pos,]
inheritedMhap <- inheritedMhap[duplicated(inheritedMhap$pos)==FALSE,]
colnames(inheritedMhap) <- c("pos", "Mgeno")

inheritedMhap$Pgeno <- inheritedPhap$Pgeno
interestingSNPlist <- c(type_1$pos, type_3$pos, RHDO_snps$pos)

inheritedMhap <- inheritedMhap[inheritedMhap$pos%in%interestingSNPlist,]
inheritedMhap$var <- ifelse(inheritedMhap$Mgeno==0&inheritedMhap$Pgeno==0, "homalt",
                            ifelse(inheritedMhap$Mgeno==1&inheritedMhap$Pgeno==1, "homref", "het"))
childgenotype <- inheritedMhap

noninheritedPhap<- noninheritedPhap[noninheritedPhap$pos%in%childgenotype$pos,]
noninheritedPhap <- noninheritedPhap[duplicated(noninheritedPhap$pos)==FALSE,]
noninheritedMhap <- noninheritedMhap[noninheritedMhap$pos%in%childgenotype$pos,]

childgenotype$Pnon <- noninheritedPhap[,2]
childgenotype$Mnon <- noninheritedMhap[,2]

save.image(file = paste(wd, "/fam",familynr, "_", format(Sys.time(), "%Y-%m-%d_%H%M"),".RData", sep=""))

vcf.processing.v2 <- function(filename){
  file <- paste(filename, sep="")
  input <- data.frame(read.table(file))
  
  info <- input[,8]
  info <- gsub(pattern = "DP=", replacement = "0,", x=info) #get rid of the = signs
  info <- gsub(pattern = ";I16=", replacement = ",0,", x=info)
  infotest <- strsplit(info, split=",")
  
  parseInfo <- function(n,x)
  {return(as.numeric(strsplit(x[n],split=",")[[1]][c(2,4:7)]))}
  
  allInfo <- t(sapply(1:length(info),parseInfo,x=info))
  info <- as.data.frame(allInfo)
  
  chr <- input[,1]
  pos <- input[,2]
  refallele <- input[,4]
  altallele <- input[,5]
  readcount <- info[,1]
  reftot <- info[,2]+info[,3]
  alttot <- info[,4]+info[,5]
  
  
  Important.info <- data.frame(chr,pos,refallele,altallele,readcount,reftot,alttot)
  #Important.info$refp <- round(100*((Important.info$reftot-Important.info$alttot)/(Important.info$reftot+Important.info$alttot)),2)
  Important.info$refp <- round(100*(Important.info$reftot/Important.info$readcount),1)
  
  return(Important.info)
}

CVS <- vcf.processing.v2(CVSfile)
CVS <- subset(CVS, CVS$pos%in%childgenotype$pos)
CVS <- subset(CVS, CVS$readcount>mincov_CVS)


childgenotype$CVS <- "X"
childgenotype$CVSreads <- 0

for(i in 1:nrow(childgenotype)){
  pos <- childgenotype[i,1]
  child <- CVS[CVS$pos==pos,]
  if(nrow(child)==1){
    childvar <- ifelse(child$refp>hom_high, "homref",
                       ifelse(child$refp<hom_low, "homalt",
                              ifelse(child$refp>20&child$refp<80,"het", "unknown")))
    readcount <- child$readcount}
  
  if(nrow(child)==0){childvar<- "unknown"
  readcount <- 0}
  
  childgenotype[i,7] <- childvar
  childgenotype[i,8] <- readcount
  
}

figuredata <- childgenotype

figuredata$var <- ifelse(figuredata$var=="homref", 1,
                         ifelse(figuredata$var=="homalt", 0,
                                ifelse(figuredata$var=="het", 2, 3)))
figuredata$CVS <- ifelse(figuredata$CVS=="homref", 1,
                         ifelse(figuredata$CVS=="homalt", 0,
                                ifelse(figuredata$CVS=="het", 2, 3)))
figuredata$result <- ifelse(figuredata$CVS==3, 3,
                            ifelse(figuredata$CVS==figuredata$var, 3, 4))
figurematrix <- data.matrix(figuredata[,c(6,5,2,3,4,7,9)])

pdf(file=paste(wd, "/family", familynr,"_heatmap.pdf", sep=""))
heatmap.2(figurematrix, breaks=c(0, 0.5,1.5, 2.5, 3.5, 4.5), col = c("green", "blue", "orange", "white", "red"),
          Colv=F, Rowv=F,labRow=F, trace="none", xlab=paste(sep="",min(figuredata$pos), "-",max(figuredata$pos)), ylab=paste("SNPs n=", nrow(figurematrix), sep=""), 
          cexCol=1, labCol=c("Mnon", "Pnon", "Minh", "Pinh","Pred", "CVS", "res") ,            
          key=F, main=paste("family ",familynr, ifelse(phap1inherited==T, " phap1", " phap2"), 
                            ifelse(mhap1inherited==T, " mhap1", " mhap2"), sep=""))
dev.off()

heatmap.2(figurematrix, breaks=c(0, 0.5,1.5, 2.5, 3.5, 4.5), col = c("green", "blue", "orange", "white", "red"),
          Colv=F, Rowv=F,labRow=F, trace="none", xlab="", ylab=paste("SNPs n=", nrow(figurematrix), sep=""), 
          cexCol=1, labCol=c("Mnon", "Pnon", "Minh", "Pinh","Pred", "CVS", "res") ,            
          key=F, main=paste("family ",familynr, ifelse(phap1inherited==T, " phap1", " phap2"), 
                            ifelse(mhap1inherited==T, " mhap1", " mhap2"), sep=""))

save.image(file = paste(wd, "/fam",familynr, "_", format(Sys.time(), "%Y-%m-%d_%H%M"),".RData", sep=""))

stats

a
b
c
d
e
f
g
H
k
L

total_A
total_B
