### Step1 SNP
step1=function(file_list){
  len=29409
  stand=read.table("WH.fas",head=F,sep=",")[2,1]
  for(i in 1:length(file_list)){
    a=read.table(paste0(file_list[i],".fas"),head=F,sep=",")
    l=nrow(a)
    p=seq(2,l,2)
    q=a[p,1]
    judge=rep(0,len)
    for(j in 1:length(q)){
      for(k in 1:len){
        judge[k]=judge[k]+(substring(q[j],k,k)!=substring(stand,k,k))
      }
    }
    result=data.frame(position=1:len,mutation=judge/length(q))
    mutate=result[which(result$mutation>=0.9),1]
    write.table(mutate,paste0(file_list[i],"-pos.txt"),row.names=F,col.names=F,sep=",")
    write.table(result$mutation,paste0(file_list[i],"-mut.csv"),row.names=F,col.names=T,sep=",")
  }
}


### Step2 
step2=function(file_list){
  nt=c("a","t","c","g","-")
  stand=read.table("WH.fas",head=F,sep=",")[2,1]
  fudong=0.01
  for(uu in 1:length(file_list)){
    pos1=read.table(paste0(file_list[uu],"-pos.txt"))[,1]
    pos2=read.table(paste0(file_list[uu],"-pos.txt"))[,1]
    a=read.table(paste0(file_list[uu],".fas"),head=F,sep=",")
    l=nrow(a)
    p=seq(2,l,2)
    q=a[p,1]
    name_max=paste0(file_list[uu],"-liansuo.csv")
    for(i in 1:(length(pos1)-1))
      for(j in (i+1):length(pos2)){
        mut_no=0
        mut_one=0
        mut_all=0
        p1=pos1[i]
        p2=pos2[j]
        m=matrix(0,nrow=5,ncol=5)
        for(k in 1:length(q)){
          pp1=(substring(q[k],p1,p1)!=substring(stand,p1,p1))
          pp2=(substring(q[k],p2,p2)!=substring(stand,p2,p2))
          if((pp1==0)&&(pp2==0)){
            mut_no=mut_no+1
          }
          else if((pp1==1)&&(pp2==1)){
            m1=substring(q[k],p1,p1)
            m2=substring(q[k],p2,p2)
            x=which(nt==m1)
            y=which(nt==m2)
            m[x,y]=m[x,y]+1
            mut_all=mut_all+1
          }
          else{
            mut_one=mut_one+1
          }
        }
        if((length(q)-mut_no)==0)next
        if(mut_one/(length(q)-mut_no)>fudong)next
        maxx=max(m)
        maxr=maxx/(length(q)-mut_no)
        s1=paste0(as.character(p1),"-",as.character(p2))
        for(ii in 1:4)
          for(jj in 1:4){
            if(m[ii,jj]>0){
              if((maxx==m[ii,jj])&&maxr>(1-fudong)){
                s=paste0(substring(stand,p1,p1),"->",nt[ii],"|",substring(stand,p2,p2),"->",nt[jj])
                write.table(s1,name_max,append=T,row.names=F,col.names=F,sep=",")
                write.table(s,name_max,append=T,row.names=F,col.names=F,sep=",")
                write.table(m[ii,jj],name_max,append=T,row.names=F,col.names=F,sep=",")
                write.table(maxr,name_max,append=T,row.names=F,col.names=F,sep=",")
              }
            }
          }
      }
  }
}


### Step3 
step3=function(file_list){
  for(uu in 1:length(file_list)){
    a=read.csv(paste0(file_list[uu],"-liansuo.csv"),header=F)[,1]
    mutation=read.csv(paste0(file_list[uu],"-mut.csv"))[,1]
    num=seq(3,length(a),4)
    num=a[num]
    uniq_num=unique(num)
    print(uniq_num)
    for(i in 1:length(uniq_num)){
      equal=which(uniq_num[i]==a)
      len=length(equal)*2
      site=rep("",len)
      mut=rep("",len)
      m=matrix(0,nrow=3,ncol=len)
      for(j in 1:(len/2)){
        ss=strsplit(a[equal[j]-2],"-",fixed = T)
        m[1,2*j-1]=ss[[1]][1]
        m[1,2*j]=ss[[1]][2]
        mm=strsplit(a[equal[j]-1],"|",fixed = T)
        m[2,2*j-1]=mm[[1]][1]
        m[2,2*j]=mm[[1]][2]
        
        m[3,2*j-1]=mutation[as.integer(ss[[1]][1])]
        m[3,2*j]=mutation[as.integer(ss[[1]][2])]
      }
      index=duplicated(m[1,])
      write.table(m[,!index],paste0(file_list[uu],"-unit.csv"),append=T,row.names=F,col.names=F,sep=",")
    }
  }
}




