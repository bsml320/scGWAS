library("gplots")
library("pheatmap")

all_panels = read.table("Analysis/R1/0.summary/R1_scGWASv4.all_panels.txt", header=T, as.is=T)
traits.anno = read.table("traits.txt", as.is=T)

which(all_panels[,4] < 0) -> ii
all_panels = cbind(all_panels[ii, ], module_size=all_panels[-ii, "size"])

#all_panels = all_panels[which(all_panels[,"module_size"] > 10),]
print(dim(all_panels))

sig.res = c()
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "DER20CC"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "DER22"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Saunders"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Zeisel_Lvl5"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Lung10x"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Lung"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Madissoon_lung"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Heart"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Liver"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Madissoon_esophagus"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Pancreas_E-MTAB-5061"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Pancreas_GSE81547"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Pancreas_GSE81608"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Pancreas_GSE84133"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Pancreas_GSE85241"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "PBMC10k"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "Madissoon_spleen"),] )
sig.res = rbind(sig.res, all_panels[which(all_panels$panel == "VentoTormo"),] )

sig.res[which(sig.res$panel == "DER22"), "panel"] = "Brain_DER22"
sig.res[which(sig.res$panel == "DER20CC"), "panel"] = "Brain_DER20"
sig.res[which(sig.res$panel == "Saunders"), "panel"] = "Brain_Saunders"
sig.res[which(sig.res$panel == "Zeisel_Lvl5"), "panel"] = "Brain_Zeisel"
sig.res[which(sig.res$panel == "Lung"), "panel"] = "Lung_SS2"
sig.res[which(sig.res$panel == "Lung10x"), "panel"] = "Lung_10x"
sig.res[which(sig.res$panel == "Madissoon_lung"), "panel"] = "Madissoon_Lung"
sig.res[which(sig.res$panel == "Madissoon_esophagus"), "panel"] = "Esophagus"
sig.res[which(sig.res$panel == "Madissoon_spleen"), "panel"] = "Spleen"
sig.res[which(sig.res$panel == "VentoTormo"), "panel"] = "Decidua"


traits = sort(unique(sig.res$trait))
cells = (unique(  paste(sig.res$panel, sig.res$cell, sep=":")  ))
p.matrix = matrix(0, nrow=length(traits), ncol=length(cells))
rownames(p.matrix) = traits
colnames(p.matrix) = cells
for(cell in cells){
	for(trait in traits){
		which(  paste(sig.res$panel, sig.res$cell, sep=":") == cell & sig.res$trait == trait) -> ii
		if(length(ii) > 0){
			p.matrix[trait, cell] = sig.res[ii, "raw"]
		}
	}
}


sapply(cells, function(u){
	strsplit(u, split=":")[[1]][1]
}) -> panels

annotation_col = data.frame(
  panel = factor(unlist(panels))
)
rownames(annotation_col) = cells
head(annotation_col)

reds = colorRampPalette(c("red", "pink"))
greens = colorRampPalette(c("#66A61E", "lightgreen"))
rains = colorRampPalette(c("blue", "lightblue"))

ann_colors = list(
  panel = c(Brain_DER20 = reds(4)[1], Brain_DER22 = reds(4)[2], Brain_Saunders = reds(4)[3], Brain_Zeisel = reds(4)[4], 
            Lung_10x = greens(3)[1], Lung_SS2 = greens(3)[2], Madissoon_Lung = greens(3)[3], 
            Heart = "purple",Liver = "magenta", Esophagus = "cyan", 
			"Pancreas_E-MTAB-5061" = rains(5)[1], Pancreas_GSE81547 = rains(5)[2], Pancreas_GSE81608 = rains(5)[3], Pancreas_GSE84133 = rains(5)[4], Pancreas_GSE85241 = rains(5)[5],
			PBMC10k = "orange", Spleen = "mistyrose4", Decidua = "darkgoldenrod"
			)
)
head(ann_colors)

colfunc <- colorRampPalette(c("white", "black"))
cc = colfunc( 100*max( p.matrix ) )
cc[1:5] = "snow"

pheatmap(p.matrix, cluster_cols=F, cluster_rows=T, color=cc, annotation_col=annotation_col, labels_col=rep("", ncol(p.matrix) ), annotation_colors = ann_colors)

dev.off()

pdf("Analysis/R1/Figure3/Figure3A.scGWASv4_all.pdf", width=16, height=7)
pheatmap(p.matrix, cluster_cols=F, cluster_rows=T, color=cc, annotation_col=annotation_col, labels_col=rep("", ncol(p.matrix) ), annotation_colors = ann_colors)
dev.off()

