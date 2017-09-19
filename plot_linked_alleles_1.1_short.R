#link the alleles to each other

###### standard_families
root <- "~/resources/fastqs/DFR/"

wd <- root
sifname <- paste("DFR_CFTR_dbsnp.sif", sep="")


plot.new()

sif <- read.table(paste(wd, sifname, sep="/"))
sif <- subset(sif, sif$V2=="pp")

links <- data.frame()
link1 <-strsplit(as.character(sif$V1), split=":")
link2 <-strsplit(as.character(sif$V3), split=":")
for(i in 1:length(link1)){
  links[i,1] <- as.numeric(link1[[i]][1])
  links[i,2] <- as.character(link1[[i]][2])
  links[i,3] <- as.numeric(link2[[i]][1])
  links[i,4] <- as.character(link2[[i]][2])
  links[i,5] <- as.character(link1[[i]][3])
}
haplotype1_links <- subset(links, links$V5=="hap1")
haplotype2_links <- subset(links, links$V5=="hap2")
#haplotype1_links[,1] <- as.numeric(gsub(haplotype1_links[,1], pattern = 508, replacement = 117199645))
#haplotype2_links[,1] <- as.numeric(gsub(haplotype2_links[,1], pattern = 508, replacement = 117199645))
l1 <- haplotype1_links
l2 <- haplotype2_links
#l1 <- read.delim(file1, h=F) ###haplotype 1
#l2 <- read.delim(file2, h=F) ###haplotype 2

#png(paste(wd, paste(parent, "_spiderplot.png", sep=""), sep="/"), wid=1200, hei=300)
pdf(paste(wd, paste(parent, "_spiderplot.pdf", sep=""), sep="/"))


nuc.pos <- unique(c(l1[,1], l1[,3], l2[,1], l2[,3]))
#####select rounded 5% and 95% quantile as size
size <- floor(diff(quantile(nuc.pos, c(0.05,0.95)))/1e3) 
#header <- paste(length(nuc.pos), " SNVs haplotyped over a distance of ", size, "kb (containing 90% of SNVs)", sep="")

#plot(range(c(l1[,c(1,3)],l2[,c(1,3)])), c(-3,3), type='n', axes=F, xlab="", ylab="", main=header)
#axis(1,seq(floor(min(l1[,c(1,3)])/20e3)*20e3, ceiling(max(l1[,c(1,3)])/20e3)*20e3, by=20e3))
plot(range(c(116790292,117510468)), c(-3,3), type='n', 
     axes=F, xlab="", ylab="")
axis(1, c(116790292, 117510468))
axis(3, c(117282526, 117282621,117199646, 117170953))

for(i in 1:nrow(l1)){
	xspline(c(l1[i,1], (l1[i,1]+l1[i,3])/2, l1[i,3]), c(1, 3,1), c(0,1,0), border=rgb(1,0,0,0.3), lwd=1)
}	

for(i in 1:nrow(l2)){
	xspline(c(l2[i,1], (l2[i,1]+l2[i,3])/2, l2[i,3]), c(-1, -3,-1), c(0,1,0), border=rgb(0,0,1,0.3), lwd=1)
}

col.list <- list(A = 'red', C = 'green', G = 'orange', T = 'blue', W= 'black', D='black')
#draw the top nucleotides
nuc1 <- l1[,1:2]
nuc2 <- l1[,3:4]
names(nuc2) <- names(nuc1)
nuc <- unique(rbind(nuc1, nuc2)) ####merge nucleotides of both links
nuc[,2] <- as.character(nuc[,2])
nuc <- nuc[order(nuc[,1]),] ######### order by coordinate
start.end <- range(nuc[,1])  ####### find range
i.pos <- seq(start.end[1], start.end[2], len=nrow(nuc))
text(i.pos, 0.3, nuc[,2], cex=0.4, col=unlist(col.list[nuc[,2]]))
segments(nuc[,1],0.9, i.pos, 0.5)

#draw the bottom nucleotides
nuc1 <- l2[,1:2]
nuc2 <- l2[,3:4]
names(nuc2) <- names(nuc1)
nuc <- unique(rbind(nuc1, nuc2))
nuc[,2] <- as.character(nuc[,2])
nuc <- nuc[order(nuc[,1]),]
start.end <- range(nuc[,1])
i.pos <- seq(start.end[1], start.end[2], len=nrow(nuc))
text(i.pos, -0.3, nuc[,2], cex=0.4, col=unlist(col.list[nuc[,2]]))
segments(nuc[,1],-0.9, i.pos, -0.5)

dev.off()

