#link the alleles to each other

###### standard_families


networkfile <- "~/resources/Small files/haps/HBB/M9_test_network.txt"


#plot.new()


plot_linked_alleles <- function(networkfile){
network <- read.table(networkfile, header=TRUE, stringsAsFactors = F)
#sif <- sif[,c(1,3,2,4:ncol(sif))]


cleaned <- subset(network, network$interaction=="pp")



haplotype1_links <- subset(network, network$hap==1)
haplotype2_links <- subset(network, network$hap==2)
rootsplit <- strsplit(networkfile, split="/")[[1]]
rootname <- strsplit(rootsplit[length(rootsplit)], split="_network.txt")[[1]][1]

pdf(paste0(dirname(networkfile), "/",rootname, "_spiderplot.pdf"))


nuc.pos <- unique(c(haplotype1_links$pos, haplotype1_links$pos.1, haplotype2_links$pos, haplotype2_links$pos.2))
#####select rounded 5% and 95% quantile as size
size <- floor(diff(quantile(nuc.pos, c(0.05,0.95)))/1e3) 
#header <- paste(length(nuc.pos), " SNVs haplotyped over a distance of ", size, "kb (containing 90% of SNVs)", sep="")

#plot(range(c(l1[,c(1,3)],l2[,c(1,3)])), c(-3,3), type='n', axes=F, xlab="", ylab="", main=header)
#axis(1,seq(floor(min(l1[,c(1,3)])/20e3)*20e3, ceiling(max(l1[,c(1,3)])/20e3)*20e3, by=20e3))
left_edge <- min(nuc.pos)
right_edge <- max(nuc.pos)

plot(range(c(left_edge,right_edge)), c(-3,3), type='n', 
     axes=F, xlab="", ylab="", main=rootname)
axis(1,c(left_edge, right_edge))
#axis(3, c(117282526, 117282621,117199646, 117170953))
abline(h=0)
for(i in 1:nrow(haplotype1_links)){
  dat <- haplotype1_links[i,]
	xspline(c(dat$pos, (dat$pos+dat$pos.1)/2, dat$pos.1), c(0, 3,0), c(0,1,0), border=rgb(1,0,0,0.3), lwd=1)
}	

for(i in 1:nrow(haplotype2_links)){
  dat <- haplotype2_links[i,]
  xspline(c(dat$pos, (dat$pos+dat$pos.1)/2, dat$pos.1),c(0, -3,0), c(0,1,0), border=rgb(0,0,1,0.3), lwd=1)
 
}



dev.off()
}

plot_linked_alleles("~/resources/Small files/haps/HBB/M9_test_network.txt")
