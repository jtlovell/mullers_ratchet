piesOnTrees<-function(dat=counts, pop=pop,apo,sex, neut="D4",cons="D0"){
  dat<-counts[counts$name==pop,]
  dat$newtot<-dat$totalcount
  dat$newtot[duplicated(dat$morph)]<-NA
  
  a1<-paste(apo,"hap1",sep="_")
  a2<-paste(apo,"hap2",sep="_")
  s1<-paste(sex,"hap1",sep="_")
  s2<-paste(sex,"hap2",sep="_")
  dat<-apply(dat,1,function(x){
    # part 1. reassign haplotype names, based on order in tree
    if(gregexpr(a1,x)[[2]][1]<gregexpr(a2,x)[[2]][1]){
      x<-gsub(a1,"apo1",x)
      x<-gsub(a2,"apo2",x)
    }else{
      x<-gsub(a2,"apo1",x)
      x<-gsub(a1,"apo2",x)
    }
    if(gregexpr(s1,x, fixed=T)[[2]][1]<gregexpr(s2,x)[[2]][1]){
      x<-gsub(s1,"sex1",x)
      x<-gsub(s2,"sex2",x)
    }else{
      x<-gsub(s2,"sex1",x)
      x<-gsub(s1,"sex2",x)
    }
    
    # part 2. replace brackets with parentheses so that morph can be read as a tree
    x<-gsub("[[]","z",x)
    x<-gsub("[]]","w0",x)
    x<-gsub("[[]","z",x)
    x<-gsub("[]]","w0",x)
    
    # part 3. reorder apo-sex within a group so that sex is first... makes it easier to parse trees
    x<-gsub("apo1,sex1","sex1,apo1", x)
    x<-gsub("apo2,sex2","sex2,apo2", x)
    x<-gsub("apo2,apo1","apo1,apo2", x)
    x<-gsub("sex2,sex1","sex1,sex2", x)
    x<-gsub("zapo1,apo2w0,zsex1,sex2w0","zsex1,sex2w0,zapo1,apo2w0", x)
    x<-gsub("zsex2,apo2w0,zsex1,apo1w0","zsex1,apo1w0,zsex2,apo2w0", x)
    x<-gsub("z","(",x)
    x<-gsub("w0",")",x)
    x
  })
  #regorganize dat
  dat<-data.frame(t(dat), stringsAsFactors=F)
  unique(dat$morph)
  unique(dat$branch)
  dat$totalcount<-as.numeric(dat$totalcount)
  dat$count<-as.numeric(dat$count)
  dat$newtot<-as.numeric(dat$newtot)
  
  #get totals for each tree type
  dat.tot<-ddply(dat,~morph, summarize, counts=sum(newtot, na.rm=T))
  
  #get totals/branch and snp type
  dat2<-dat[dat$type.name %in% c(cons,neut),]
  dat.sum<-ddply(dat2,~morph+type.name+branch, summarize, counts=sum(count))
  byfold<-cast(dat.sum, morph+branch~type.name, value="counts")
  
  # determine if trees are grouped or laddered
  byfold$treeshape<-ifelse(grepl("\\),\\(",byfold$morph), "grouped","ladder")
  
  # calculate the proportion conserved to non conserved
  byfold$prop.cons<-byfold[,cons]/(byfold[,cons]+byfold[,neut])
  byfold$prop.neut<-byfold[,neut]/(byfold[,cons]+byfold[,neut])
  
  par(mfrow=c(3,2), mar=c(2,1,2,1))
  
  for(i in unique(byfold$morph)){
    dat<-byfold[byfold$morph==i,]
    #make tree
    refRatio<-dat$D0[nchar(dat$branch)==25]/(dat$D0[nchar(dat$branch)==25]+dat$D4[nchar(dat$branch)==25])
    tre<-paste(pop,i,";",sep="")
    cat(tre, file = "ex.tre", sep = "\n")
    tre<-read.tree("ex.tre")
    cols<-tre$tip.label
    cols[grep("apo",cols)]<-"red"
    cols[grep("sex",cols)]<-"blue"
    cols[grep("ref",cols)]<-"black"
    tot<-as.numeric(apply(dat[,c(cons,neut)],1,sum))
    tot<-max(tot[-which(tot==max(tot))])
    plot(tre, tip.col=cols, cex=1.5)
    title(paste(pop,"  n.trees = ", dat.tot$counts[dat.tot$morph==i],sep=""), cex=4)
    legend("bottomleft",c("4D","0D"),col=c("cyan","red"), pch=19)
    if(tot==0) next
    #plot pies
    #first the terminal tips
    #internal branches
    if(dat$treeshape[1]=="ladder"){
      for(j in 1:4){
        lab<-tre$tip.label[j]
        props<-as.numeric(dat[dat$branch==lab,c("prop.cons","prop.neut")])
        edgeRatio<-props[1]/sum(props)
        sums<-sum(as.numeric(dat[dat$branch==lab,c(cons,neut)]))
        edgelabels(edge=(j+3), pie=props, cex=c(sums/tot+.5))
        if(sums>2){
          bt<-binom.test(x=as.numeric(dat[dat$branch==lab,c(cons,neut)]),n=sum,p=refRatio, alternative="greater")
          p<-ifelse(bt$p.value<0.05,"**","")
          edgelabels(edge=(j+3),text=p,bg=NULL,frame="none",adj=1.75, cex=1.5)
        }

      }
      ind1<-which(grepl(tre$tip.label[1],dat$branch) & grepl(tre$tip.label[4],dat$branch) & grepl(tre$tip.label[3],dat$branch) & grepl(tre$tip.label[2],dat$branch))
      ind2<-which(grepl(tre$tip.label[1],dat$branch) & grepl(tre$tip.label[3],dat$branch) & grepl(tre$tip.label[2],dat$branch))
      ind2<-ind2[ind2!=ind1]
      ind3<-which(grepl(tre$tip.label[1],dat$branch) & grepl(tre$tip.label[2],dat$branch))
      ind3<-ind3[!ind3 %in% c(ind1,ind2)]
      prop1<-as.numeric(dat[ind1,c("prop.cons","prop.neut")])
      prop2<-as.numeric(dat[ind2,c("prop.cons","prop.neut")])
      prop3<-as.numeric(dat[ind3,c("prop.cons","prop.neut")])
      sums1<-sum(as.numeric(dat[ind1,c(cons,neut)]))
      sums2<-sum(as.numeric(dat[ind2,c(cons,neut)])) 
      sums3<-sum(as.numeric(dat[ind3,c(cons,neut)]))
      edgelabels(edge=ind1, pie=prop1, cex=1)
      edgelabels(edge=ind2, pie=prop2, cex=c(sums2/tot+.5))
      if(sums2>2){
        bt2<-binom.test(x=as.numeric(dat[ind2,c(cons,neut)]),n=sums2,p=refRatio, alternative="greater")
        p2<-ifelse(bt2$p.value<0.05,"**","")
        edgelabels(edge=ind2,text=p2,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
      
      edgelabels(edge=ind3, pie=prop3, cex=c(sums3/tot+.5))
      if(sums3>2){
        bt3<-binom.test(x=as.numeric(dat[ind3,c(cons,neut)]),n=sums3,p=refRatio, alternative="greater")
        p3<-ifelse(bt3$p.value<0.05,"**","")
        edgelabels(edge=ind3,text=p3,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
     
    }else{
      lab1<-tre$tip.label[1]
      props1<-as.numeric(dat[dat$branch==lab1,c("prop.cons","prop.neut")])
      sum1<-sum(as.numeric(dat[dat$branch==lab1,c(cons,neut)]))
      edgelabels(edge=3, pie=props1, cex=(sum1/tot+.5))
      if(sum1>2){
        bt1<-binom.test(x=as.numeric(dat[dat$branch==lab1,c(cons,neut)]),n=sum1,p=refRatio, alternative="greater")
        p1<-ifelse(bt1$p.value<0.05,"**","")
        edgelabels(edge=3,text=p1,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
      lab2<-tre$tip.label[2]
      props2<-as.numeric(dat[dat$branch==lab2,c("prop.cons","prop.neut")])
      sum2<-sum(as.numeric(dat[dat$branch==lab2,c(cons,neut)]))
      edgelabels(edge=4, pie=props2, cex=c(sum2/tot+.5))
      if(sum2>2){
        bt2<-binom.test(x=as.numeric(dat[dat$branch==lab2,c(cons,neut)]),n=sum2,p=refRatio, alternative="greater")
        p2<-ifelse(bt2$p.value<0.05,"**","")
        edgelabels(edge=4,text=p2,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
      lab3<-tre$tip.label[3]
      props3<-as.numeric(dat[dat$branch==lab3,c("prop.cons","prop.neut")])
      sum3<-sum(as.numeric(dat[dat$branch==lab3,c(cons,neut)]))
      edgelabels(edge=6, pie=prop3, cex=(sum3/tot+.5))
      if(sum3>2){
        bt3<-binom.test(x=as.numeric(dat[dat$branch==lab3,c(cons,neut)]),n=sum3,p=refRatio, alternative="greater")
        p3<-ifelse(bt3$p.value<0.05,"**","")
        edgelabels(edge=6,text=p3,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
      lab4<-tre$tip.label[4]
      props4<-as.numeric(dat[dat$branch==lab4,c("prop.cons","prop.neut")])
      sum4<-sum(as.numeric(dat[dat$branch==lab4,c(cons,neut)]))
      edgelabels(edge=7, pie=props4, cex=sum4/tot+.5)
      if(sum4>2){
        bt4<-binom.test(x=as.numeric(dat[dat$branch==lab4,c(cons,neut)]),n=sum4,p=refRatio, alternative="greater")
        p4<-ifelse(bt4$p.value<0.05,"**","")
        edgelabels(edge=7,text=p4,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
      
      ind1<-which(grepl(tre$tip.label[1],dat$branch) & grepl(tre$tip.label[4],dat$branch) & grepl(tre$tip.label[3],dat$branch) & grepl(tre$tip.label[2],dat$branch))
      ind2<-which(grepl(tre$tip.label[1],dat$branch) & grepl(tre$tip.label[2],dat$branch))
      ind2<-ind2[ind2!=ind1]
      ind3<-which(grepl(tre$tip.label[3],dat$branch) & grepl(tre$tip.label[4],dat$branch))
      ind3<-ind3[!ind3 %in% c(ind1,ind2)]
      prop1<-as.numeric(dat[ind1,c("prop.cons","prop.neut")])
      prop2<-as.numeric(dat[ind2,c("prop.cons","prop.neut")])
      prop3<-as.numeric(dat[ind3,c("prop.cons","prop.neut")])
      sum1<-sum(as.numeric(dat[ind1,c(cons,neut)]))
      sum2<-sum(as.numeric(dat[ind2,c(cons,neut)]))
      sum3<-sum(as.numeric(dat[ind3,c(cons,neut)]))
      edgelabels(edge=1, pie=prop1,cex=1)
      
      edgelabels(edge=2, pie=prop2,cex=c(sum2/tot+.5))
      if(sum2>2){
        bt2<-binom.test(x=as.numeric(dat[ind2,c(cons,neut)]),n=sum2,p=refRatio, alternative="greater")
        p2<-ifelse(bt2$p.value<0.05,"**","")
        edgelabels(edge=2,text=p2,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
      edgelabels(edge=5, pie=prop3,cex=c(sum3/tot+.5))
      if(sum3>2){
        bt3<-binom.test(x=as.numeric(dat[ind3,c(cons,neut)]),n=sum3,p=refRatio, alternative="greater")
        p3<-ifelse(bt3$p.value<0.05,"**","")
        edgelabels(edge=5,text=p3,bg=NULL,frame="none",adj=1.75, cex=1.5)
      }
    }
  }
  return(list(statsByBranch=byfold,countsPerTree=dat.tot))
}


