### calculate topology distances of the asex-tree and sex-tree to the unconstrained tree (this is done in R)###
setwd("/YOUR/DIRECTORY/WITH/TREES/")
library(ape)
PATH="/YOUR/DIRECTORY/WITH/TREES/"
files=list.files(PATH)
for (file in files){
	temp=read.tree(paste(PATH,file,sep="/"))
	out=dist.topo(temp,y=NULL,method="score")
	write.table(as.matrix(out),paste(PATH,file,"dist",sep=""))
}

### combine the results as one table (done in Unix) ###
sed -i -e "1d" *dist
awk 'NR == 1 || NR % 3 == 1' *dist > awked
cat awked | awk '{print$0,$3-$4}' > awksum

