#Human proteomics from Johnson et al. 2020
prot.dat<-fread("prot.dat_data/3.cleanDat.csv",stringsAsFactors=F)
colnames(prot.dat)[1]<-"id"
prot.dat$prot_id<-gsub("(.+)\\|.+","\\1",prot.dat$id)

#add metadata
meta<-fread(
  "Consensus.MASTER-Table 1.tsv",
  stringsAsFactors=F)

ad.id<-meta$`RAW File Name`[meta$Group%in%c("AD","AsymAD")]
control.id<-meta$`RAW File Name`[grep("Control",meta$Group)]

case<-ad.id[ad.id%in%colnames(prot.dat)]
control<-control.id[control.id%in%colnames(prot.dat)]

prot.dat$FC<-0
prot.dat$pval<-1

for(i in 1:nrow(prot.dat)){
  set(prot.dat,
      i=i,
      j="FC",
      value = tryCatch(log2(
        mean(2^(as.numeric(
          prot.dat[i,case,with=FALSE])),na.rm=TRUE))-log2(
            mean(2^(as.numeric(
              prot.dat[i,control,with=FALSE])),na.rm=TRUE)),error=function(e){1}))
  
  set(prot.dat,
      i=i,
      j="pval",
      value = tryCatch(t.test(prot.dat[i,control,with=FALSE],
                              prot.dat[i,case,with=FALSE])$p.val,error=function(e){1}))
}

prot.dat$qval<-p.adjust(prot.dat$pval,'fdr')

prot.comm.sig<-prot.dat[prot.dat$qval<0.0001,]

human.prot.commodity<-data.table(id=prot.comm.sig$prot_id,
                                 prize_val=-log10(prot.comm.sig$qval),
                                 source="human proteomics",
                                 magnitude=prot.comm.sig$FC,stringsAsFactors=F)

## drosophila proteomics and metabolomics
protCols <- as.character(
  read_excel("Supplementary table 4.xlsx", 
             n_max = 1,sheet=1,skip=6,col_names = FALSE))
protCols[12:20]<-c("Control 1","Control 2","Control 3",
                   "Abeta 1","Abeta 2","Abeta 3",
                   "tau 1","tau 2","tau 3")
dro.dat<-read_excel("Supplementary table 4.xlsx",
                    sheet=1,skip=7,col_names=protCols)
dro.dat<-as.data.table(dro.dat)

dro.dat$prot_name<-gsub(".+\\|(.+)\\|.+","\\1",dro.dat$`Protein ID (Uniprot)`)

# predicted orthologs from DIOPT
ortho.pred<-fread("auxiliary_files/ortholog_pred.xls",sep='\t',
                  stringsAsFactors=F)
ortho.pred<-ortho.pred[,"max_score":=max(`DIOPT Score`),by="Fly GeneID"]
ortho.pred.noDup<-ortho.pred[,
                             c("Fly GeneID","Fly Symbol","Human Symbol","Rank",
                               "DIOPT Score","max_score"),with=FALSE]
ortho.pred.noDup<-ortho.pred.noDup[-which(ortho.pred.noDup$Rank%in%c("low",
                                                                     "None"))]

human.dro<-merge(dro.dat,ortho.pred.noDup,by.x="prot_name",
                 by.y="Fly GeneID")
human.dro$pval_tau<-1
human.dro$FC_tau<-(-.Machine$double.xmax)
human.dro$pval_abeta<-1
human.dro$FC_abeta<-(-.Machine$double.xmax)

control<-c("Control 1","Control 2","Control 3")
case<-c("tau 1","tau 2","tau 3"
)

for(i in 1:nrow(human.dro)){
  set(human.dro,
      i=i,
      j="FC_tau",
      value = tryCatch(log2(
        mean(as.numeric(
          human.dro[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              human.dro[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(human.dro,
      i=i,
      j="pval_tau",
      value = tryCatch(t.test(human.dro[i,control,with=FALSE],
                              human.dro[i,case,with=FALSE])$p.val,error=function(e){1}))
}

control<-c("Control 1","Control 2","Control 3")
case<-c("Abeta 1","Abeta 2","Abeta 3"
)

for(i in 1:nrow(human.dro)){
  set(human.dro,
      i=i,
      j="FC_abeta",
      value = tryCatch(log2(
        mean(as.numeric(
          human.dro[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              human.dro[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(human.dro,
      i=i,
      j="pval_abeta",
      value = tryCatch(t.test(human.dro[i,control,with=FALSE],
                              human.dro[i,case,with=FALSE])$p.val,error=function(e){1}))
}

human.dro$qval_tau<-p.adjust(human.dro$pval_tau,method="fdr")
human.dro$qval_abeta<-p.adjust(human.dro$pval_abeta,method="fdr")

eqt<-as.data.table(read_excel("Supplementary_table3.xlsx"))
colnames(eqt)[which(colnames(eqt)=="Egene3")]<-"Gene"

#handle duplicates by choosing smallest qvalue
human.dro[,"prize_qval":=min(qval_tau,qval_abeta),by=seq_len(nrow(human.dro))]
human.dro.noDup<-human.dro[order(human.dro$`Human Symbol`,
                                 human.dro$prize_qval)]
human.dro<-human.dro.noDup[!duplicated(human.dro.noDup$`Human Symbol`),]

human.dro$orth_cap<-as.numeric(human.dro$`DIOPT Score`)/as.numeric(
  human.dro$max_score)

prot.abeta.commodity<-data.table(id=human.dro$`Human Symbol`[
  human.dro$qval_abeta<0.01
],
prize_val=-log10(human.dro$qval_abeta[
  human.dro$qval_abeta<0.01
]),
type="drosophila ABeta proteomics",
magnitude=human.dro$FC_abeta[
  human.dro$qval_abeta<0.01
],stringsAsFactors=F)

prot.commodity.tau<-data.table(id=human.dro$`Human Symbol`[
  human.dro$qval_tau<0.01
],
prize_val=-log10(human.dro$qval_tau[
  human.dro$qval_tau<0.01
]),
type="drosophila tau proteomics",
magnitude=human.dro$FC_tau[
  human.dro$qval_tau<0.01
],stringsAsFactors=F)

prot.commodity<-rbind(prot.abeta.commodity,
                      prot.commodity.tau)
prot.commodity<-prot.commodity[prot.commodity$id!=""]

prot.commodity<-prot.commodity[order(prot.commodity$prize_val),]
prot.commodity<-prot.commodity[!duplicated(prot.commodity$id),]

## output from piumet
metab<-fread("auxiliary_files/metab_qval_annotate.txt")

metab.sig<-metab[metab$qval_abeta<0.1|metab$qval_tau<0.1]
metab.sig<-metab.sig[grepl("^HMDB",metab.sig$`HMDB ID`),]
metab<-data.table(id=metab.sig$`HMDB ID`,
                  prize_val=(-log10(metab.sig$qval_tau)),
                  type="metabolomics",
                  magnitude=metab.sig$logFC_Tau,stringsAsFactors=F)

## eGene prizes
sink.commod<-data.table(id=eqt$Gene,
                        prize_val=-log10(eqt$FDR),
                        source="eGene",
                        magnitude=eqt$BETA_FE_discov,
                        stringsAsFactors=F)

#hits from Lohr et al. 2020
tau.prize<-fread("Table1.txt",
                 sep='\t',header=T,stringsAsFactors=F)

tau.prize<-tau.prize[tau.prize$`human gene`!=""]

screen.prize<-data.table(id=tau.prize$`human gene`,
                         prize_val=1,
                         source="neurodegeneration genes",
                         magnitude=0,
                         stringsAsFactors=F)

tau.screen<-fread("auxiliary_files/tau_modifiers.tsv",
                  stringsAsFactors=F)

tau.screen<-tau.screen[tau.screen$`human homolog`!=""]

tau.prize<-data.table(id=tau.screen$`human homolog`,
                      prize_val=1,
                      source="tau fly genetic screen",
                      magnitude=0,
                      stringsAsFactors=F)

## generate prizes for phosphoproteomics and upstream kinases
#characterize kinase activity on significantly altered
#phosphosites

## read in differential phosphorylation from 7_drosophila_proteomics_phosphoproteomics.R
phospho.gene<-fread("differential_drosophila_phosphoproteomics.txt")

upstream<-fread("auxiliary_files/upstream_kinases",
                stringsAsFactors=F)

phospho.gene.up<-merge(phospho.gene,upstream,by.x="Fly Symbol",
                       by.y="Fly_Symbol")

#test phosphosites
phospho.gene[phospho.gene$prize_qval<0.1&
               !duplicated(phospho.gene[,
                                        c("Fly Symbol","site_position_1")]),c(
                                          "Fly Symbol","Id","site_position_1"
                                        )][1:100,]

#limit to significant phosphosites
phospho.gene$FC<-0
phospho.gene$FC[phospho.gene$prize_qval[
  is.na(phospho.gene$prize_qval)==FALSE]==phospho.gene$qval_tau[
    is.na(phospho.gene$prize_qval)==FALSE]]<-
  phospho.gene$FC_tau[phospho.gene$prize_qval[
    is.na(phospho.gene$prize_qval)==FALSE]==phospho.gene$qval_tau[
      is.na(phospho.gene$prize_qval)==FALSE]]

phospho.gene$FC[phospho.gene$prize_qval[
  is.na(phospho.gene$prize_qval)==FALSE]==phospho.gene$qval_abeta[
    is.na(phospho.gene$prize_qval)==FALSE]]<-
  phospho.gene$FC_abeta[phospho.gene$prize_qval[
    is.na(phospho.gene$prize_qval)==FALSE]==phospho.gene$qval_abeta[
      is.na(phospho.gene$prize_qval)==FALSE]]

phospho.gene$source<-"tau drosophila phosphoproteomics"
phospho.gene$source[phospho.gene$qval_abeta==phospho.gene$prize_qval]<-"abeta drosophila phosphoproteomics"

phospho.gene.sig<-phospho.gene[phospho.gene$prize_qval<0.1&abs(phospho.gene$FC)>1]
phospho.gene.sig$`Human Symbol`<-gsub("MAPT","p-Tau",phospho.gene.sig$`Human Symbol`)

phospho.unique.fly<-phospho.gene.up[!duplicated(
  phospho.gene.up[,c("Fly Symbol","upstream_kinase")]),]

phospho.profile<-phospho.unique.fly[, lapply(.SD, median, na.rm=TRUE), 
                                    by=upstream_kinase, 
                                    .SDcols=c("control 1", "control 2", "control 3",
                                              "Abeta 1","Abeta 2","Abeta 3",
                                              "tau 1","tau 2","tau 3") ] 

# identify phosphosites correlated with the proteomic abundance of 
# upstream kinases
cor.df<-data.frame()
case<-c(
  "tau 1","tau 2","tau 3"
)

for(prot in unique(phospho.profile$upstream_kinase[
  phospho.profile$upstream_kinase%in%human.dro$`Human Symbol`
])){
  temp.df<-data.frame(kinase=prot,cor_val=cor(
    as.numeric(phospho.profile[phospho.profile$upstream_kinase%in%prot,
                               case,with=FALSE]),
    as.numeric(human.dro[human.dro$`Human Symbol`%in%prot,
                         case,with=FALSE]),method="spearman"
  ),source="tau upstream kinase")
  
  cor.df<-rbind(cor.df,temp.df)
}

case<-c("Abeta 1","Abeta 2","Abeta 3"
)

phospho.profile<-phospho.unique.fly[, lapply(.SD, median, na.rm=TRUE), 
                                    by=upstream_kinase, 
                                    .SDcols=c(
                                      "control 1", "control 2", "control 3",
                                              "Abeta 1","Abeta 2","Abeta 3",
                                              "tau 1","tau 2","tau 3"
                                    ) ] 

for(prot in unique(phospho.profile$upstream_kinase[
  phospho.profile$upstream_kinase%in%human.dro$`Human Symbol`
])){
  temp.df<-data.frame(kinase=prot,cor_val=cor(
    as.numeric(unlist(phospho.profile[phospho.profile$upstream_kinase%in%prot,
                                      case,with=FALSE])),
    as.numeric(unlist(human.dro[human.dro$`Human Symbol`%in%prot,
                                case,with=FALSE])),method="spearman"
  ),source="abeta upstream kinase")
  
  cor.df<-rbind(cor.df,temp.df)
}

phospho.profile$correlation_status<-"uncorrelated"

phospho.profile$correlation_status[
  phospho.profile$upstream_kinase%in%cor.df$kinase[
    cor.df$cor_val>0.4]]<-"correlated"

phospho.profile$correlation_status[
  phospho.profile$upstream_kinase%in%cor.df$kinase[
    cor.df$cor_val<(-0.4)]]<-"anticorrelated"

row.anno<-data.frame(
  correlation=phospho.profile$correlation_status)
rownames(row.anno)<-phospho.profile$upstream_kinase

cohort.anno<-gsub(" [0-9]","",c(control,case))
#cohort.anno<-gsub("Control","control",cohort.anno)

col.anno<-data.frame(disease_status=cohort.anno)
rownames(col.anno)<-gsub("Control","control",c(control,case))

phospho.profile<-phospho.profile[!duplicated(phospho.profile$upstream_kinase),
                                 cohort.anno,with=FALSE]

### generate prizes for omicsIntegrator ###
phos.prot<-phospho.gene.sig[!rev(duplicated(rev(phospho.gene.sig$`Human Symbol`))),]

phos.prize<-data.table(id=phos.prot$`Human Symbol`,
                       prize_val=(-log10(phos.prot$prize_qval)),
                       source=phos.prot$source,
                       magnitude=phos.prot$FC,
                       stringsAsFactors=F)

kin.prof<-phospho.unique.fly[, lapply(.SD, min, na.rm=TRUE), 
                             by=upstream_kinase, 
                             .SDcols=c("prize_qval","FC") ]

cor.kin<-cor.df$kinase[abs(cor.df$cor_val)>0.4]
kin.cor<-kin.prof[kin.prof$upstream_kinase%in%cor.kin,]

kin.prize<-data.table(id=kin.cor$upstream_kinase,
                      prize_val=(-log10(kin.cor$prize_qval)),
                      source="upstream kinase",
                      magnitude=kin.cor$FC,
                      stringsAsFactors=F)

phos.kin.prize<-rbind(phos.prize,kin.prize)

### lipidomics from Blanchard et al. ###
lipid<-fread("Lipidomics_analysis.tsv",sep='\t',
             stringsAsFactors=F)
case<-c("Sample 1","Sample 2","Sample 3","Sample 4",
        "Sample 5")
control<-c("Sample 7","Sample 8","Sample 9","Sample 10",
           "Sample 11","Sample 12")

lipid$female_43_FC<-1
lipid$female_43_pval<-1

lipid$male_43_FC<-1
lipid$female_43_pval<-1

lipid$combined_43_FC<-1
lipid$combined_43_pval<-1

for(i in 1:nrow(lipid)){
  set(lipid,
      i=i,
      j="combined_43_FC",
      value = tryCatch(log2(
        mean(as.numeric(
          lipid[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              lipid[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(lipid,
      i=i,
      j="combined_43_pval",
      value = tryCatch(t.test(lipid[i,control,with=FALSE],
                              lipid[i,case,with=FALSE])$p.val,error=function(e){1}))
}

case<-c("Sample 1","Sample 2","Sample 3")
control<-c("Sample 7","Sample 8","Sample 9")

for(i in 1:nrow(lipid)){
  set(lipid,
      i=i,
      j="female_43_FC",
      value = tryCatch(log2(
        mean(as.numeric(
          lipid[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              lipid[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(lipid,
      i=i,
      j="female_43_pval",
      value = tryCatch(t.test(lipid[i,control,with=FALSE],
                              lipid[i,case,with=FALSE])$p.val,error=function(e){1}))
}

case<-c("Sample 4",
        "Sample 5")
control<-c("Sample 10",
           "Sample 11","Sample 12")

for(i in 1:nrow(lipid)){
  set(lipid,
      i=i,
      j="male_43_FC",
      value = tryCatch(log2(
        mean(as.numeric(
          lipid[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              lipid[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(lipid,
      i=i,
      j="male_43_pval",
      value = tryCatch(t.test(lipid[i,control,with=FALSE],
                              lipid[i,case,with=FALSE])$p.val,error=function(e){1}))
}

lipid$combined_q<-p.adjust(lipid$combined_43_pval,method='fdr')
lipid$male_q<-p.adjust(lipid$male_43_pval,method='fdr')
lipid$female_q<-p.adjust(lipid$female_43_pval,method='fdr')

lipid[,'combined_pval':=min(`Comparison 3/4 F`,`Comparison 3/4M`,
                            `Comparison 3/4 ALL`),
      by=seq_len(nrow(lipid))]

lipid$source<-"combined APOE4 lipidomics"
lipid$source[
  lipid$combined_pval==lipid$`Comparison 3/4 F`]<-"female-specific APOE4 lipidomics"
lipid$source[
  lipid$combined_pval==lipid$`Comparison 3/4M`]<-"male-specific APOE4 lipidomics"

lipid$FC<-0

lipid$FC[lipid$combined_pval[
  is.na(lipid$combined_pval)==FALSE]==lipid$`Comparison 3/4 F`[
    is.na(lipid$combined_pval)==FALSE]]<-
  lipid$female_43_FC[lipid$combined_pval[
    is.na(lipid$combined_pval)==FALSE]==lipid$`Comparison 3/4 F`[
      is.na(lipid$combined_pval)==FALSE]]

lipid$FC[lipid$combined_pval[
  is.na(lipid$combined_pval)==FALSE]==lipid$`Comparison 3/4M`[
    is.na(lipid$combined_pval)==FALSE]]<-
  lipid$male_43_FC[lipid$combined_pval[
    is.na(lipid$combined_pval)==FALSE]==lipid$`Comparison 3/4M`[
      is.na(lipid$combined_pval)==FALSE]]

lipid$FC[lipid$combined_pval[
  is.na(lipid$combined_pval)==FALSE]==lipid$`Comparison 3/4 ALL`[
    is.na(lipid$combined_pval)==FALSE]]<-
  lipid$combined_43_FC[lipid$combined_pval[
    is.na(lipid$combined_pval)==FALSE]==lipid$`Comparison 3/4 ALL`[
      is.na(lipid$combined_pval)==FALSE]]

lipid.prize<-data.table(id=lipid$`Molecular Species`[
  lipid$combined_pval<0.05],
  prize=-log10(lipid$combined_pval[
    lipid$combined_pval<0.05
  ]),
  source=lipid$source[
    lipid$combined_pval<0.05
  ],
  magnitude=lipid$FC[
    lipid$combined_pval<0.05
  ])

lipid.prize<-lipid.prize[is.na(lipid.prize$id)==FALSE,]

#merge with HMDB IDs#
hmdb.id<-fread("auxiliary_files/metabolites-2020-07-20",sep=',',stringsAsFactors=F)
hmdb.id$NAME<-gsub("\\([0-9]+Z\\)","",hmdb.id$NAME)
hmdb.id$NAME<-gsub("\\([0-9]+Z,[0-9]+Z\\)","",hmdb.id$NAME)
hmdb.id$NAME<-gsub("\\([0-9]+Z,[0-9]+Z,[0-9]+Z\\)","",hmdb.id$NAME)
hmdb.id$NAME<-gsub("\\([0-9]+Z,[0-9]+Z,[0-9]+Z,[0-9]+Z\\)","",hmdb.id$NAME)



#format metabolites for merging
lipid.prize$id<-gsub("_","/",lipid.prize$id)
lipid.prize$id<-gsub("hexosylCer","GlcCer",
                     lipid.prize$id)

lipid.hmdb<-merge(lipid.prize,hmdb.id,by.x="id",
                  by.y="NAME")

lipid.prize.hmdb<-lipid.hmdb[,c("HMDB_ID","prize",
                                "source","magnitude")]

colnames(lipid.prize.hmdb)[1]<-"id"
colnames(lipid.prize.hmdb)[2]<-"prize_val"
lipid.prize.hmdb$id[grep("hexosyl",lipid.hmdb$id)]<-lipid.hmdb$id[
  grep("hexosyl",lipid.hmdb$id)
]

#combine all prize files
master.prize<-rbind(prot.commodity,
                    human.prot.commodity,
                    metab,
                    sink.commod,
                    screen.prize,
                    tau.prize,
                    phos.kin.prize,
                    lipid.prize.hmdb)

### add gwas hits ###

gwas<-fread("auxiliary_files/GWAS_hits.tsv",stringsAsFactors=F)
gwas.no.dup<-gwas[gwas$`Reported gene`%in%prizes$id==FALSE&
                    duplicated(gwas$`Reported gene`)==FALSE]
gwas.no.dup<-gwas.no.dup[-which(gwas.no.dup$`Reported gene`%in%c("","-"))]

gwas.prize<-data.frame(id=gwas.no.dup$`Reported gene`,
                       prize_val=1,source="GWAS locus",
                       magnitude=0)

gwas.prize$id<-gsub("\\s","",gwas.prize$id)
prize.add<-rbind(master.prize,gwas.prize)

#supplemental files from Corces et al. 2020
epi.loci<-fread("auxiliary_files/SupTable2-4 Locus Analysis-Table 1.tsv",
                stringsAsFactors=F)
ad.epi<-epi.loci[grep("AD",epi.loci$Disease),]
ad.epi$Affected_Gene<-gsub(" Locus","",ad.epi$Affected_Gene)

ad.epi.gene<-ad.epi[ad.epi$`Potential_Coding_Effect?`==TRUE&
                      ad.epi$Affected_Gene%in%prize.add$id==FALSE,]

gwas.loci<-fread("auxiliary_files/SupTable2-2 All GWAS SNPs-Table 1.tsv",
                 stringsAsFactors=F)

gwas.q.add<-merge(ad.epi,gwas.loci,by.x="Locus_Lead_SNP",
                  by.y="SNP_rsID")

gwas.q.add$GWAS_q<-p.adjust(as.numeric(gwas.q.add$GWAS_pvalue),method="bonferroni",
                            n=1e6)
gwas.q.add$GWAS_q[29]<-p.adjust(0.0000000398,method="bonferroni",
                                n=1e6)

coding.prize<-data.frame(id=gwas.q.add$Affected_Gene[
  gwas.q.add$`Potential_Coding_Effect?`==TRUE
],
prize_val=-log10(gwas.q.add$GWAS_q[
  gwas.q.add$`Potential_Coding_Effect?`==TRUE
]),source="accessible GWAS locus",
magnitude=0)

coding.prize$prize_val[coding.prize$prize_val==Inf]<-(-log10(.Machine$double.xmin))

prize.add.epi<-rbind(prize.add,coding.prize)

prize.add.epi$prize_val[is.na(prize.add.epi$prize_val)]<-0
for(sources in unique(prize.add.epi$source)){
  temp.df<-prize.add.epi[prize.add.epi$source==sources]
  if(all(prize.add.epi$prize_val[prize.add.epi$source==sources]!=1)){
    
    prize.add.epi$prize_val[prize.add.epi$source==sources]<-(
      temp.df$prize_val-min(temp.df$prize_val))/(max(temp.df$prize_val)-min(temp.df$prize_val))
  }
}

#weighting#

prize.add.epi$prize_val[
  prize.add.epi$source%in%c(
    "fly genetic screen","tau fly genetic screen")]<-prize.add.epi$prize_val[
      prize.add.epi$source%in%c(
        "fly genetic screen","tau fly genetic screen")]*0.95

prize.add.epi$prize_val[
  prize.add.epi$source%in%c(
    "accessible GWAS locus")]<-prize.add.epi$prize_val[
      prize.add.epi$source%in%c(
        "accessible GWAS locus")]*0.9

prize.add.epi$prize_val[
  prize.add.epi$source%in%c(
    "GWAS locus")]<-prize.add.epi$prize_val[
      prize.add.epi$source%in%c(
        "GWAS locus")]*0.85

prize.add.epi$prize_val[
  prize.add.epi$source%in%c(
    "human proteomics","abeta drosophila metabolomics",
    "tau drosophila metabolomics","male-specific APOE4 lipidomics",
    "female-specific APOE4 lipidomics")]<-prize.add.epi$prize_val[
      prize.add.epi$source%in%c(
        "human proteomics","abeta drosophila metabolomics",
        "tau drosophila metabolomics","male-specific APOE4 lipidomics",
        "female-specific APOE4 lipidomics")]*0.8

prize.add.epi$prize_val[
  prize.add.epi$source%in%c(
    "abeta drosophila phosphoproteomics",
    "tau drosophila phosphoproteomics","upstream kinase")]<-prize.add.epi$prize_val[
      prize.add.epi$source%in%c(
        "abeta drosophila phosphoproteomics",
        "tau drosophila phosphoproteomics","upstream kinase")]*0.75

prize.add.epi$prize_val[
  prize.add.epi$source%in%c(
    "drosophila tau proteomics","drosophila ABeta proteomics")]<-prize.add.epi$prize_val[
      prize.add.epi$source%in%c(
        "drosophila tau proteomics","drosophila ABeta proteomics")]*0.7

write.table(prize.add.epi,"Table_S7.txt",sep='\t',quote=F,row.names=F)


