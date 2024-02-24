library(data.table)
library(readxl)

metab<-data.table()

## collect the data from metaboliates of all charges
for(i in 1:4){
  metCols <- as.character(
    read_excel("Supplementary table 6.xlsx", 
               n_max = 1,sheet=i,skip=1,col_names = FALSE))
  temp.metab<-read_excel("Supplementary table 6.xlsx",
                 sheet=i,skip=2,col_names=metCols)
  temp.metab<-as.data.table(temp.metab)
  metab<-rbind(metab,temp.metab)
}

metab$pval_tau<-1
metab$FC_tau<-(-.Machine$double.xmax)
metab$pval_abeta<-1
metab$FC_abeta<-(-.Machine$double.xmax)

control<-c("Control_1","Control_2","Control_3")
case<-c("Tau_1","Tau_2","Tau_3"
)

for(i in 1:nrow(metab)){
  set(metab,
      i=i,
      j="FC_tau",
      value = tryCatch(log2(
        mean(as.numeric(
          metab[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              metab[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(metab,
      i=i,
      j="pval_tau",
      value = tryCatch(t.test(metab[i,control,with=FALSE],
                              metab[i,case,with=FALSE])$p.val,error=function(e){1}))
}

control<-c("Control_1","Control_2","Control_3")
case<-c("Abeta_1","Abeta_2","Abeta_3"
)

for(i in 1:nrow(metab)){
  set(metab,
      i=i,
      j="FC_abeta",
      value = tryCatch(log2(
        mean(as.numeric(
          metab[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              metab[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(metab,
      i=i,
      j="pval_abeta",
      value = tryCatch(t.test(metab[i,control,with=FALSE],
                              metab[i,case,with=FALSE])$p.val,error=function(e){1}))
}

metab$qval_tau<-p.adjust(metab$pval_tau,method="fdr")
metab$qval_abeta<-p.adjust(metab$pval_abeta,method="fdr")

metab[,"prize_qval":=min(qval_tau,qval_abeta),by=seq_len(nrow(metab))]


# piumet input format
piu.in<-data.table(id=metab$`m/z`[
  metab$prize_qval<0.1
],
charge=metab$Charge[
  metab$prize_qval<0.1
],
prize=-log10(metab$prize_qval[
  metab$prize_qval<0.1
]),stringsAsFactors=F)

#write file, then upload to PiuMet: https://fraenkel-nsf.csbi.mit.edu/piumet2/#
write.table(piu.in,"piumet_input_abeta_tau.txt",sep='\t',quote=F,
            row.names=F,col.names=F)
