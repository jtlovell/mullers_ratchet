parseCount_summary<-function(file, name){
  dat<-read.delim(file,sep=" ", head=F, stringsAsFactors=F)
  index<-grep("reference",dat$V1)
  all.out<-data.frame()
  for(i in 1:(length(index)-1)){
    i1<-index[i]
    i2<-index[i+1]
    counts<-dat[(i1+3):(i2-2),1:3]
    colnames(counts)<-c("branch","cat","count")
    morph<-as.character(dat[i1,1])
    totalcount<-as.numeric(as.character(dat[i1+1,4]))
    type<-dat[i1-1,]
    type<-as.character(apply(type[,1:5],1,function(x) paste(x, collapse="")))
    counts$type.name<-as.character(multigsub(as.numeric(as.character(key[,2])), as.character(key[,1]), counts$cat))
    out<-data.frame(name=name,morph,type,totalcount,counts, stringsAsFactors=F)
    all.out<-rbind(all.out,out)
  }
  return(all.out)
}