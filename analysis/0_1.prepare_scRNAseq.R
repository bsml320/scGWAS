

dat = read.table("DER-22_Single_cell_expression_raw_UMI.tsv", as.is=T, header=T)

cell_type = colnames(dat)
sapply(cell_type, function(u)strsplit(u, split="\\.")[[1]][1]) -> ct

## calculate log2(cpm)
lib_size <- apply(dat[,-1], 2, sum)/1000000
log_cpm <- log2(t(apply(dat[,-1], 1, function(x){x/lib_size}))+1)
dim(log_cpm)
rownames(log_cpm) = dat[,1]
apply(log_cpm, 1, function(u)sum(u==0)/length(u)) -> rowCheck
log_cpm = log_cpm[which(rowCheck < 0.95),]

## calculate average expression per cell type
cell_type = sort(unique(ct[-1]))

avg_log_cpm = matrix(0, nrow=nrow(log_cpm), ncol=length(cell_type))
colnames(avg_log_cpm) = cell_type

for(k in 1:nrow(log_cpm)){
	tapply(log_cpm[k, ], ct[-1], mean) -> u
	avg_log_cpm[k, ] = u[cell_type]
}
rownames(avg_log_cpm) = rownames(log_cpm)

qua.cell_type = setdiff(colnames(avg_log_cpm), "NA" )
avg_log_cpm = avg_log_cpm[, qua.cell_type]

write.table(cbind(gene=rownames(avg_log_cpm), avg_log_cpm), file="DER-22.avg.tsv", row.names=F, quote=F, sep="\t")
