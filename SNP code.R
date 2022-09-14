SNP=function(reg,filepath,outpath,leng){ #reg=region, filepath=input file path, outpath=output file path, leng=DNA seq length
library(plyr)
region=reg
for(k in 1:length(region))
{
  fn=paste(region[k],"-full-SNP.csv",sep="") 
  path=filepath[k]
  a=read.table(path,head=F,sep=",") 
  
  standard=a[2,1]
  len=leng
  L=nrow(a)
  p=seq(4,L,2)
  q=a[p,1] 
  judge=matrix(rep(0,len*2),ncol=2) 
  for (i in 1:len){
    count=0
    mutation=c()
    for(j in 1:length(q)){
      
      if(substring(q[j],i,i)!=substring(standard,i,i)){ 
        mutation=c(mutation,paste(substring(standard,i,i),"->",substring(q[j],i,i),sep="")) 
        count=count+1
      }
    }
    mutation1=count(mutation)
    muta=c()
    for (u in 1:nrow(mutation1)){muta=c(muta,paste(mutation1[u,1],mutation1[u,2]))}
    if(count/length(q)<0.01){
      judge[i,1]=0
    }
    else{judge[i,1]=count/length(q)}
    judge[i,2]=paste(muta,collapse=" ") 
  }
  
  
  result=data.frame(position=1:len,mutation_rate=judge[,1],mutation=judge[,2]) 
  output=result[which(result$mutation_rate>0.01),]
  write.table(output,paste(outpath,fn,sep=""),row.names=F,col.names=T,sep=",")
}
}