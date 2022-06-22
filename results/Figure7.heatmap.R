library("ggplot2")

col_1 = colorRampPalette(c("white", "red"))
col_red = c(rep("white", 16), col_1(23)) 

all_panels = read.table("results/R1_scGWASv4.all_panels.txt", header=T, as.is=T)

data2 = all_panels[all_panels[,2] == "Liver", ]

dat1 = data2[which(data2$adjust < 0),]
dat1[,4] = dat1[,5]

dat2 = data2[which(data2$adjust >= 0),]
dat2$NES = floor(dat2$adjust * 10 )/10
dat2_col = col_red[floor(dat2$NES*10)]

max_z = max(ceiling(data2$raw))
max_size = max(data2$size)

#######################################################################
  
g1 = ggplot(dat1,
       aes(x = cell, y = trait, fill = adjust, size=max_size ) ) +
       geom_point(shape = 21) + scale_x_discrete(drop=FALSE) + scale_y_discrete(drop=FALSE) + 
       scale_fill_gradient2(low="green4", high="lightgreen", limits=c(5, max_z) ) + 
       geom_point(data = dat2, mapping=aes(x=cell, y=trait, size=size * 0.8, fill= NES), colour=dat2_col ) + 
       theme(panel.grid.major = element_line(linetype = 2, color = "grey"), axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12), axis.text.y = element_text(size=12)  ) + 
       ylab("") + xlab("") + labs(fill='z (prop.\ntext)') 
       guides(size = guide_legend(title="# module\ngenes")   ) 

pdf("Figure7.Liver.pdf", width=5.9, height=6.9)
print(g1)
dev.off()

