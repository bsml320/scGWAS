library("vegan")
library("corrplot")
library("fgsea")

setwd("D:/BIG/work/19-scGWAS")
source("source.trait.R")
local.path = "D:/BIG/work/19-scGWAS/"

info = c()
sig.genes = c()
sig.module = c()

input_header = ""
out_header = ""

ff = paste(input_header,"..shared_genes.scores.calibration.txt", sep="")
all.weight.data = read.table(ff, as.is=T, header=T)

ff = paste(input_header,".combined-z.txt", sep="")
dat = read.delim(ff, header=T, as.is=T)
dat = dat[!is.na(dat[,8]),]
dat$module_score_z = as.numeric(dat$module_score_z)
dat = dat[order(dat$cell_type, dat$module_score, decreasing=T),]
cells = sort(unique(dat$cell_type))

fr = paste(input_header, ".all.random_modules.txt", sep="")
datr = read.table(fr, header=F, as.is=T)

colfunc <- colorRampPalette(c("grey", "black"))
cc = colfunc( max( round(dat[,"module_score_z"]*100), na.rm=T ) )
	
### plot 1.
jpeg(paste(out_header,".jpg", sep=""), width=n*410, height=n*410, res=200)
n = ceiling(sqrt(length(cells)))
par(mfrow=c(n,n), mar=c(3,3,3,1))
if(n * (n-1) >= length(cells))par(mfrow=c(n-1,n), mar=c(3,3,3,1))

for(cell in cells){
	dat1 = dat[which(dat$cell_type==cell), ]
	dat1 = dat1[!is.na(dat1[,8]),]
	
	datrc = datr[which(datr$V8==cell), ]
	sapply(datrc[,1], function(u)length(strsplit(u, split=":")[[1]])) -> datr.size
	sapply(dat1[,1],  function(u)length(strsplit(u, split=":")[[1]])) -> dat1.size
	which(datr.size %in% dat1.size) -> ii
	datrc = datrc[ii, ]
	datr.size = datr.size[ii]
	
	zn = ze = zm = rep(0, nrow(datrc))
	for(ss in unique(datr.size)){
		which(datr.size == ss) -> ii
		zn[ii] = (datrc[ii, 5] - mean(datrc[ii, 5]))/sd(datrc[ii, 5])
		ze[ii] = (datrc[ii, 6] - mean(datrc[ii, 6]))/sd(datrc[ii, 6])
		zm[ii] = (datrc[ii, 4] - mean(datrc[ii, 4]))/sd(datrc[ii, 4])
	}
	which( (!is.na(zn)) & (!is.na(ze)) ) -> ii
	datrc = datrc[ii, ]
	zn = zn[ii]; ze = ze[ii]; zm = zm[ii]
	
	####################################################
	### estimate significance
	####################################################
	cn = quantile(zn[!is.na(zn)], probs=.95)
	ce = quantile(ze[!is.na(ze)], probs=.95)
	length( which(dat1[,8] > cn & dat1[,9] > ce) ) -> p1; 
	length( which(zn > cn & ze > ce) ) -> p2; 
	if(p1 > 3){
		p_hat = (p1 + p2)/(nrow(dat1) + nrow(datrc))
		p1 = p1/nrow(dat1)
		p2 = p2/nrow(datrc)
		z1 = (p1 - p2)/sqrt(p_hat * (1-p_hat) *(1/nrow(dat1) + 1/nrow(datrc)) )
	} else {
		z1 = 0
	}
	
	
	####################################################
	### define sig modules
	####################################################
	
	p = c()
	for(kk in 1:nrow(dat1))p = c(p, sum(zm > dat1[kk,"module_score_z"])/(length(zm) ))
	
	p_gwas = c()
	for(kk in 1:nrow(dat1))p_gwas = c(p_gwas, sum(zn > dat1[kk,"z_gwas"])/(length(zn) ))
	
	p_scrn = c()
	for(kk in 1:nrow(dat1))p_scrn = c(p_scrn, sum(ze > dat1[kk,"z_scrnaseq"])/(length(ze) ))

	sig.ii = which(p < 0.05 & p_gwas < 0.05 & p_scrn < 0.05)
	if(length(sig.ii) > 0){
		sig.module = rbind(sig.module, cbind(dat1[sig.ii, ], trait = trait, panel = panel, p=p[sig.ii], p_gwas = p_gwas[sig.ii], p_scrn = p_scrn[sig.ii]))
	}
	
	unique(unlist(sapply(dat1[sig.ii, 1], function(u)strsplit(u, split=":")[[1]]))) -> genes
	genes = as.vector(genes)
	genes = unique(genes)
	
	if(length(genes) > 5){
		sig.genes = rbind(sig.genes, c(  paste(sort(genes), collapse=":"), trait, cell) )
		weight.data = all.weight.data[all.weight.data[,10] == cell, ]
		
		gwas.data = weight.data[,8]; names(gwas.data) = weight.data[,1]
		res2 <- fgsea(pathways = list(A=genes), stats = gwas.data, minSize=5,maxSize=1000,nperm=10000, nproc=1)
		
		scrn.data = weight.data[,9]; names(scrn.data) = weight.data[,1]
		res4 <- fgsea(pathways = list(A=genes), stats = scrn.data, minSize=5,maxSize=1000,nperm=10000, nproc=1)
		
		info = rbind(info, c(trait, panel, cell, z1, length(genes), nrow(dat1), res2[1,c(2,4,5,7)], res4[1,c(2,4,5,7)] ))
	} else {
		info = rbind(info, c(trait, panel, cell, z1, length(genes), nrow(dat1), rep(0, 8)  ))
	}
	####################################################
	
	dot_col = rep("grey", nrow(dat1))
	which(dat1[,"module_score_z"] > 0) -> idx
	dot_col[idx] = cc[round(  1+dat1[idx,"module_score_z"]*100)]
	
	xmin = min(c(dat1[,"z_gwas"], zn), na.rm=T)
	xmax = max(c(dat1[,"z_gwas"], zn), na.rm=T)
	ymin = min(c(dat1[,"z_scrnaseq"], ze), na.rm=T)
	ymax = max(c(dat1[,"z_scrnaseq"], ze), na.rm=T)
	
	plot(0,0, col="white", xlab="", ylab="", 
	     main=paste(cell, "\nGWAS=", format( as.numeric(info[nrow(info), 9]), digits=3), "; scRNAseq=", format(  as.numeric(info[nrow(info), 13]), digits=3), "\nz1=", format(z1, digits=3), "", sep=""), 
		 pch=19, xlim=c(xmin, xmax), ylim=c(ymin, ymax) )
	points(zn, ze, col="lightblue", pch=19, cex=.7)
	points(dat1[, "z_gwas"], dat1[, "z_scrnaseq"], col=dot_col)
	points(dat1[sig.ii, "z_gwas"], dat1[sig.ii, "z_scrnaseq"], col="red", cex=1.5, pch=19)
	abline(v=quantile(ze[!is.na(ze)], probs=.95), lty=2, col="red")
	abline(h=quantile(zn[!is.na(zn)], probs=.95), lty=2, col="red")
	
	datc = rbind(cbind(zn, ze), cbind(dat1[, "z_gwas"], dat1[, "z_scrnaseq"])	 )
	ordiellipse(datc, groups=c(rep(1,nrow(datrc) ),rep(2, nrow(dat1) )),col=c(1:2), display = "sites", kind = "sd", label = F, conf = 0.95, lty=5,lwd=0.5)
	
	mtext("GWAS", 1, line=2)
	mtext("scRNAseq", 2, line=2)
}
dev.off()

write.table(info, file=paste(out_header, ".summary.info.txt", sep=""), row.names=F, quote=F, sep="\t")
write.table(sig.genes, file=paste(out_header, ".sig_genes.txt", sep=""), row.names=F, quote=F, sep="\t")
write.table(sig.module, file=paste(out_header, ".sig_module.txt", sep=""), row.names=F, quote=F, sep="\t")
