library("xlsx")

# create fc. with filtering >1.3 or <-1.3
#For each patient:
#1. Histology
#2. Mutations
#3. For each simplified pathway, the mean value of the filtered/unfiltered pathway
#4. Deciles-based score
#5. Unified score: any mutation = 10; No mutation = expression score

#	Take a look at the CNV file: if there is amplification, is it correlated with RNA?
#   microRNA: for each gene in the simplified pathways, there are the 5 miRNAs that target it. Perhaps correct 

#	Write an equation for a score that adjusts for miRNA

printf<-function(s,...) {
print(sprintf(s,...)); 
}

smart.fold_change<-function(x,cutoff) {
	threshold=log(cutoff,base=2);
	for (i in 1:ncol(x)) {
		x[abs(x[,i])>threshold,i]=c(NA);
	}
	return(x);
}
raw.mut.data= raw.mut.data=read.xlsx(file="mutations-31Jan2014.xls",as.is=T,head=T, 1 , stringsAsFactors=F);
rownames(raw.mut.data)=raw.mut.data[,"Patid"];

raw.mRNA.data=read.csv(file="mRNA_pathway-31Jan2014.csv",as.is=T,head=T);
row.names(raw.mRNA.data)=raw.mRNA.data[,1];
raw.mRNA.data=raw.mRNA.data[,-1];
raw.mRNA.data = raw.mRNA.data[-1,]
#raw.mRNA.data[1:5, 1:5]
#log(raw.mRNA.data[1:5, 1:5])
names1 = rownames(raw.mRNA.data)
names2 = colnames(raw.mRNA.data)
rownames(raw.mRNA.data) = 1:121
colnames(raw.mRNA.data) = 1:181
raw.mRNA.matrix = as.matrix(raw.mRNA.data[1:121, 1:181])
mRNA.matrix = matrix(as.numeric(raw.mRNA.matrix), nrow = 121, ncol = 181)
mRNA.data.log=log(mRNA.matrix,base=2);

rank.matrix=mRNA.data.log;
rownames(rank.matrix) = names1
colnames(rank.matrix) = names2

for (i in 1:ncol(rank.matrix)) {
	rank.matrix[,i]=rank(rank.matrix[,i]);
}	

groupings=sub(colnames(rank.matrix),pat="_.*",rep="",perl=T);
groupings.unique=unique(groupings);

mean.fc=NULL;
mean.fc.filtered=NULL;
median.fc=NULL;
mean.ranks=NULL;
median.ranks=NULL;
for (group in groupings.unique) {
	slice.ranks=data.frame(mRNA.data.log[,groupings==group]);
	slice.fc= data.frame(mRNA.data.log[,groupings==group]);
  
	mean.ranks=cbind(mean.ranks, apply(slice.ranks,1,mean));
	median.ranks=cbind(median.ranks,apply(slice.ranks,1,median));
  
	mean.fc=cbind(mean.fc,apply(slice.fc,1,mean));
	filtered.fc=smart.fold_change(slice.fc,1.3);
	mean.fc.filtered=cbind(mean.fc.filtered,apply(filtered.fc,1,function(x){mean(na.omit(x))}));
	median.fc=cbind(median.fc,apply(slice.fc,1,median));
}



colnames(mean.ranks)=groupings.unique;
colnames(median.ranks)=groupings.unique;
colnames(mean.fc)=groupings.unique;
colnames(mean.fc.filtered)=groupings.unique;
colnames(median.fc)=groupings.unique;

pdf(file="Combos_test_mean_rank.pdf");

for (col in groupings.unique) {
	order.for.plot=order(as.numeric(mean.fc[,col]));
#	plot(as.numeric(mean.ranks[order.for.plot,col]),main=sprintf("%s, mean-rank",col),xlab="Patient order",ylab="mean of ranks in pathway");
#	plot(as.numeric(median.ranks[order.for.plot,col]),main=sprintf("%s, median-rank",col),xlab="Patient order",ylab="median of ranks in pathway");
	plot(as.numeric(mean.fc[order.for.plot,col]),main=sprintf("%s, mean-log2 fold-change",col),xlab="Patient order",ylab="mean of ranks in pathway");
	plot(as.numeric(mean.fc.filtered[order(as.numeric(mean.fc.filtered[,col])),col]),main=sprintf("%s, mean-log2 fold-change with filter 1.3",col),xlab="Patient order (sorted)",ylab="mean of ranks in pathway");
#	plot(as.numeric(median.fc[order.for.plot,col]),main=sprintf("%s, median-log2 fold-change",col),xlab="Patient order",ylab="mean of ranks in pathway");
}

medians=apply(mRNA.data.log,2,median);
cvs=apply(mRNA.data.log,2,sd)/apply(mRNA.data.log,2,mean);
selected.columns=abs(medians)>=2;
boxplot(mRNA.data.log[,selected.columns]);
dev.off();

known.mutations=c("KRAS","EFGR","PIK3CA","BRAF","ERBB2","P53");
results.col.headings=c(
	"PID",
	"Histology",
	sprintf("W3.score.%s",groupings.unique),
	sprintf("mut.%s",known.mutations),
	sprintf("mRNA.%s",groupings.unique),
	sprintf("scoreR.%s",groupings.unique)
);

results.matrix=matrix(NA,ncol=length(results.col.headings),nrow=nrow(raw.mRNA.data)
                      ,dimnames=list(rownames(raw.mRNA.data),results.col.headings));
for (pid in rownames(raw.mRNA.data)) {
	printf("Starting to work on patient %s",pid);
	results.matrix[pid,"PID"]=pid;
	results.matrix[pid,"Histology"]=raw.mut.data[pid,"Histo"];
	for (mut in (known.mutations)) {
		cname=sprintf("mut.%s",mut);
		if (raw.mut.data[pid,mut]=="WT") {
			results.matrix[pid,cname]="";
		} else {
			results.matrix[pid,cname]=raw.mut.data[pid,mut];
		}
	}
	
	for (col in groupings.unique) {
		cname=sprintf("mRNA.%s",col);
		results.matrix[pid,cname]=mean.fc.filtered[pid,col];
	}
}

for (col in groupings.unique) {
	cname=sprintf("scoreR.%s",col);
	scores=cut(mean.fc.filtered[,col],10,labels=1:10);
	results.matrix[,cname]=scores;
}

for (col in groupings.unique) {
	score.cname=sprintf("W3.score.%s",col);
	mut.cname=sprintf("mut.%s",col);
	mRNA.cname=sprintf("scoreR.%s",col);
	if (sum(mut.cname %in% colnames(raw.mut.data))==0) {
		scores=-as.numeric(results.matrix[,mRNA.cname]);
	} else {
		scores=rep("10",nrow(results.matrix));
		scores[results.matrix[,mut.cname]==""]=results.matrix[results.matrix[,mut.cname]=="",mRNA.cname];
		
	}
	results.matrix[,score.cname]=scores;
}

write.csv(file="results.matrix.csv",results.matrix);






















