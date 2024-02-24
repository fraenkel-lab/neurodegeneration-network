library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(data.table)
library(Hmisc)
library(tidytext)
library(readxl)

### phosphoproteomics
pyCols <- as.character(
  read_excel("Supplementary table 5.xlsx", 
             n_max = 1,sheet=3,skip=6,col_names = FALSE))
py<-read_excel("Supplementary table 5.xlsx",
               sheet=3,skip=7,col_names=pyCols)
py<-as.data.table(py)
pyCols[15:23]<-c("control 1","control 2","control 3","Abeta 1","Abeta 2",
                 "Abeta 3","tau 1","tau 2","tau 3")
py$group<-"phospho tyrosine"


stCols <- as.character(read_excel("Supplementary table 5.xlsx", 
  n_max = 1,sheet=1,skip=6,col_names = FALSE))
stCols[c(5,6,7)]<-c("motif1","motif2","motif3")
stCols[21:29]<-c("control 1","control 2","control 3","Abeta 1","Abeta 2",
                 "Abeta 3","tau 1","tau 2","tau 3")
st<-read_excel("Supplementary table 5.xlsx",
               sheet=1,skip=7,col_names=stCols)
st<-as.data.table(st)
st$group<-"phospho serine/threonine"

motif.list<-unique(c(py$motif1,st$motif1,st$motif2[
  is.na(st$motif2)==FALSE
],st$motif3[is.na(st$motif3)==FALSE],gg$motif1))

motif.df<-rbind(py,st,fill=TRUE)

motif.df$prot_name<-gsub(".+\\|(.+)\\|.+","\\1",motif.df$Id)

motif.df$pval_tau<-1
motif.df$FC_tau<-(-.Machine$double.xmax)
motif.df$pval_abeta<-1
motif.df$FC_abeta<-(-.Machine$double.xmax)

control<-c("control 1","control 2","control 3")
case<-c(#"Abeta 1","Abeta 2","Abeta 3",
  "tau 1","tau 2","tau 3"
)

#abeta comparisons
for(i in 1:nrow(motif.df)){
  set(motif.df,
      i=i,
      j="FC_tau",
      value = tryCatch(log2(
        mean(as.numeric(
          motif.df[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              motif.df[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(motif.df,
      i=i,
      j="pval_tau",
      value = tryCatch(t.test(motif.df[i,control,with=FALSE],
                              motif.df[i,case,with=FALSE])$p.val,error=function(e){1}))
}

control<-c("control 1","control 2","control 3")
case<-c("Abeta 1","Abeta 2","Abeta 3"
)

#tau comparisons
for(i in 1:nrow(motif.df)){
  set(motif.df,
      i=i,
      j="FC_abeta",
      value = tryCatch(log2(
        mean(as.numeric(
          motif.df[i,case,with=FALSE]),na.rm=TRUE))-log2(
            mean(as.numeric(
              motif.df[i,control,with=FALSE]),na.rm=TRUE)),error=function(e){1}))
  
  set(motif.df,
      i=i,
      j="pval_abeta",
      value = tryCatch(t.test(motif.df[i,control,with=FALSE],
                              motif.df[i,case,with=FALSE])$p.val,error=function(e){1}))
}

motif.df$qval_tau<-p.adjust(motif.df$pval_tau,method="fdr")
motif.df$qval_abeta<-p.adjust(motif.df$pval_abeta,method="fdr")

motif.df[,"prize_qval":=min(qval_tau,qval_abeta),by=seq_len(nrow(motif.df))]
motif.df<-motif.df[order(motif.df$prize_qval),]

#convert to human orthologs using DIOPT
ortho.pred<-fread("auxiliary_files/phosphoprot_ids.xls",
                  stringsAsFactors = F)
ortho.pred<-ortho.pred[,"max_score":=max(`DIOPT Score`),by="Fly GeneID"]
ortho.pred.noDup<-ortho.pred[,
                             c("Fly GeneID","Fly Symbol","Human Symbol","Rank",
                               "DIOPT Score","max_score"),with=FALSE]
ortho.pred.noDup<-ortho.pred.noDup[-which(ortho.pred.noDup$Rank%in%c("low",
                                                                     "None"))]

phospho.gene<-merge(motif.df,ortho.pred.noDup,by.x="prot_name",
                    by.y="Fly GeneID",allow.cartesian=TRUE)

## save output to file for 8_multi_omic_prize_generation.R
write.table(phospho.gene,"differential_drosophila_phosphoproteomics.txt",
            sep='\t',quote=F,row.names=F)

## all significant phosphoproteomics at FDR<0.1 heatmap
phospho.gene.heat.sig<-phospho.gene[!duplicated(phospho.gene[,c("Fly Symbol","group")])&
                                      (phospho.gene$qval_abeta<0.1|
                                         phospho.gene$qval_tau<0.1)]

## significant phosphorylation events in neurodegeneration modifiers
nd.mod<-fread("fly_genetic_screens/degeneration genes for MM.tsv",
              stringsAsFactors=F)
nd.mod<-nd.mod[!is.na(nd.mod$`human gene`)]
nd.mod<-nd.mod[-which(nd.mod$`human gene`%in%c(""," "))]

phospho.gene.nd<-phospho.gene[(phospho.gene$`Fly Symbol`%in%nd.mod$`Drosophila gene`|
                                 phospho.gene$`Human Symbol`%in%nd.mod$`human gene`)&
                                (phospho.gene$qval_abeta<0.1|
                                   phospho.gene$qval_tau<0.1)]

## filter ortholog mappings that are different in the phospho and screening data
orths.keep<-data.table(`Fly Symbol`=c("scrib","rdgB","Sod","RyR","Apc","santa-maria",
                                      "dop","Lmpt","NaCP60E"),
                       `Human Symbol`=c("SCRIB","PITPNM2","SOD1","RYR2","APC",
                                        "CD36","MAST1","FHL2","SCN8A"))

phospho.gene.nd.orth<-merge(phospho.gene.nd,orths.keep,by=c("Fly Symbol",
                                                            "Human Symbol"))

phospho.gene.nd.orth$phospho_label<-paste0(phospho.gene.nd.orth$`Fly Symbol`,
                                           "_ST",phospho.gene.nd.orth$site_position_1)
phospho.gene.nd.orth$phospho_label[
  !is.na(phospho.gene.nd.orth$site_position_2)]<-paste0(
    phospho.gene.nd.orth$`Fly Symbol`[!is.na(phospho.gene.nd.orth$site_position_2)],
    "_ST",phospho.gene.nd.orth$site_position_1[!is.na(phospho.gene.nd.orth$site_position_2)],
    "/ST",phospho.gene.nd.orth$site_position_2[!is.na(phospho.gene.nd.orth$site_position_2)])

phospho.gene.nd.orth$phospho_label[
  !is.na(phospho.gene.nd.orth$site_position3)]<-paste0(
    phospho.gene.nd.orth$`Fly Symbol`[!is.na(phospho.gene.nd.orth$site_position3)],
    "_ST",phospho.gene.nd.orth$site_position_1[!is.na(phospho.gene.nd.orth$site_position3)],
    "/ST",phospho.gene.nd.orth$site_position_2[!is.na(phospho.gene.nd.orth$site_position3)],
    "/ST",phospho.gene.nd.orth$site_position3[!is.na(phospho.gene.nd.orth$site_position3)])

abeta.phos.heat<-data.matrix(phospho.gene.nd.orth[,c("Abeta 1","Abeta 2","Abeta 3",
                                                     "control 1","control 2","control 3"),with=FALSE])
tau.phos.heat<-data.matrix(phospho.gene.nd.orth[,c("tau 1","tau 2","tau 3",
                                                   "control 1","control 2","control 3"),with=FALSE])

rownames(abeta.phos.heat)<-phospho.gene.nd.orth$phospho_label
rownames(tau.phos.heat)<-phospho.gene.nd.orth$phospho_label

row_ha=rowAnnotation(
  `Significant in Abeta flies`=ifelse(phospho.gene.nd.orth$qval_abeta<0.1,"Yes","No"),
  `Significant in Tau flies`=ifelse(phospho.gene.nd.orth$qval_tau<0.1,"Yes","No"),
  col=list(
    `Significant in Abeta flies`=c("Yes"="black",
                                   "No"="grey"),
    `Significant in Tau flies`=c("Yes"="black",
                                 "No"="grey"))
) 

col_ha_abeta=HeatmapAnnotation(`Fly model`=c(rep("Amyloid beta",3),
                                             rep("Control",3)),
                               col=list(`Fly model`=c("Amyloid beta"="#c47708",
                                                      "Control"="#4d4d4d")),
                               show_annotation_name = FALSE)

col_ha_tau=HeatmapAnnotation(`Fly model`=c(rep("Tau",3),
                                           rep("Control",3)),
                             col=list(`Fly model`=c("Tau"="#4a4270",
                                                    "Control"="#4d4d4d")))

#scale data relative to control
scaled.abeta<-t(scale(t(abeta.phos.heat)))
scaled.tau<-t(scale(t(tau.phos.heat)))

## get hierarchical clustering of the whole data and apply it to the heatmaps
total.dat<-data.matrix(phospho.gene.nd.orth[,c("Abeta 1","Abeta 2","Abeta 3",
                                               "tau 1","tau 2","tau 3",
                                               "control 1","control 2","control 3"),with=FALSE])

#dist mat of scaled data
euclidean<-dist(t(scale(t(total.dat))))
hcluster<-hclust(euclidean)

pdf("Figure_3E.pdf",
    width=9,height=5.56)
p1<-Heatmap(scaled.abeta[hcluster$order,],
            name="Z-scored Protein abundance",
            #col=f2,
            top_annotation=col_ha_abeta,
            left_annotation=row_ha,
            show_row_names = FALSE,
            show_column_names=FALSE,
            column_names_rot = 60,
            cluster_rows = FALSE,
            cluster_columns = FALSE)


p2<-Heatmap(scaled.tau[hcluster$order,],
            #col=f2,
            top_annotation=col_ha_tau,
            show_row_names = TRUE,
            show_column_names=FALSE,
            column_names_rot = 60,
            #cluster_rows = TRUE,
            cluster_columns = FALSE,
            #row_names_gp = gpar(fontsize = 20),
            #column_names_gp = gpar(fontsize = 16))
            show_heatmap_legend=FALSE)

heat.list<-p1+p2
heat.list
dev.off()

### proteomics
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

#predicted orthologs from DIOPT
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
case<-c(#"tau 1","tau 2","tau 3"
  "Abeta 1","Abeta 2","Abeta 3"
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

prot.sig<-human.dro[!duplicated(human.dro[,c("Fly Symbol")])&
                      (human.dro$qval_abeta<0.1|
                         human.dro$qval_tau<0.1)]

## overlap of differential proteomics with ND modifiers
prot.nd<-human.dro[(human.dro$`Fly Symbol`%in%nd.mod$`Drosophila gene`|
                      human.dro$`Human Symbol`%in%nd.mod$`human gene`)&
                     (human.dro$qval_abeta<0.1|
                        human.dro$qval_tau<0.1)]

## filter ortholog mappings that are different in the phospho and screening data
orths.keep<-data.table(`Fly Symbol`=c("scrib","rdgB","Sod","RyR","Apc","santa-maria",
                                      "dop","Lmpt","NaCP60E","lam","gish","CG8245",
                                      "Egfr","Appl","Rala","msps","bel",
                                      "CG9934","Cont","epsilonCOP"),
                       `Human Symbol`=c("SCRIB","PITPNM2","SOD1","RYR2","APC",
                                        "CD36","MAST1","FHL2","SCN8A","LMNB1",
                                        "CSNK1G3","TMEM53","ERBB3","APP","RALA",
                                        "CKAP5","DDX3X","UBE4B","CNTN6","COPE"))

prot.nd.orth<-merge(prot.nd,orths.keep,by=c("Fly Symbol","Human Symbol"))
abeta.prot.heat<-data.matrix(prot.nd.orth[,c("Abeta 1","Abeta 2","Abeta 3",
                                             "Control 1","Control 2","Control 3"),with=FALSE])
tau.prot.heat<-data.matrix(prot.nd.orth[,c("tau 1","tau 2","tau 3",
                                           "Control 1","Control 2","Control 3"),with=FALSE])

rownames(abeta.prot.heat)<-prot.nd.orth$`Fly Symbol`
rownames(tau.prot.heat)<-prot.nd.orth$`Fly Symbol`

row_ha=rowAnnotation(
  `Significant in Abeta flies`=ifelse(prot.nd.orth$qval_abeta<0.1,"Yes","No"),
  `Significant in Tau flies`=ifelse(prot.nd.orth$qval_tau<0.1,"Yes","No"),
  col=list(
    `Significant in Abeta flies`=c("Yes"="black",
                                   "No"="grey"),
    `Significant in Tau flies`=c("Yes"="black",
                                 "No"="grey"))
) 

col_ha_abeta=HeatmapAnnotation(`Fly model`=c(rep("Amyloid beta",3),
                                             rep("Control",3)),
                               col=list(`Fly model`=c("Amyloid beta"="#c47708",
                                                      "Control"="#4d4d4d")),
                               show_annotation_name = FALSE)

col_ha_tau=HeatmapAnnotation(`Fly model`=c(rep("Tau",3),
                                           rep("Control",3)),
                             col=list(`Fly model`=c("Tau"="#4a4270",
                                                    "Control"="#4d4d4d")))

#scale data relative to control
scaled.abeta<-t(scale(t(abeta.prot.heat)))
scaled.tau<-t(scale(t(tau.prot.heat)))

## get hierarchical clustering of the whole data and apply it to the heatmaps
total.dat<-data.matrix(prot.nd.orth[,c("Abeta 1","Abeta 2","Abeta 3",
                                       "tau 1","tau 2","tau 3",
                                       "Control 1","Control 2","Control 3"),with=FALSE])

#dist mat of scaled data
euclidean<-dist(t(scale(t(total.dat))))
hcluster<-hclust(euclidean)

pdf("Figure_3D.pdf",
    width=9,height=5.56)
p1<-Heatmap(scaled.abeta[hcluster$order,],
            name="Z-scored Protein abundance",
            #col=f2,
            top_annotation=col_ha_abeta,
            left_annotation=row_ha,
            show_row_names = FALSE,
            show_column_names=FALSE,
            column_names_rot = 60,
            cluster_rows = FALSE,
            cluster_columns = FALSE)


p2<-Heatmap(scaled.tau[hcluster$order,],
            #col=f2,
            top_annotation=col_ha_tau,
            show_row_names = TRUE,
            show_column_names=FALSE,
            column_names_rot = 60,
            #cluster_rows = TRUE,
            cluster_columns = FALSE,
            #row_names_gp = gpar(fontsize = 20),
            #column_names_gp = gpar(fontsize = 16))
            show_heatmap_legend=FALSE)

heat.list<-p1+p2
heat.list
dev.off()

### metabolomics, after running PiuMet
metab<-fread("auxiliary_files/metab_qval_annotate.txt")

metab.sig<-metab[metab$qval_abeta<0.1|metab$qval_tau<0.1]


## make combined upset plots depicting proteomic overlap, metabolomic overlap, phosphoproteomic overlap
prot.up.list<-list(`Up in amyloid beta flies`=prot.sig$`Fly Symbol`[
  prot.sig$qval_abeta<0.1&prot.sig$FC_abeta>1],
  `Down in amyloid beta flies`=prot.sig$`Fly Symbol`[
    prot.sig$qval_abeta<0.1&prot.sig$FC_abeta<1],
  `Up in tau flies`=prot.sig$`Fly Symbol`[
    prot.sig$qval_tau<0.1&prot.sig$FC_tau>1],
  `Down in tau flies`=prot.sig$`Fly Symbol`[
    prot.sig$qval_tau<0.1&prot.sig$FC_tau<1])

phos.up.list<-list(`Up in amyloid beta flies`=phospho.gene.heat.sig$`Fly Symbol`[
  phospho.gene.heat.sig$qval_abeta<0.1&phospho.gene.heat.sig$FC_abeta>1],
  `Down in amyloid beta flies`=phospho.gene.heat.sig$`Fly Symbol`[
    phospho.gene.heat.sig$qval_abeta<0.1&phospho.gene.heat.sig$FC_abeta<1],
  `Up in tau flies`=phospho.gene.heat.sig$`Fly Symbol`[
    phospho.gene.heat.sig$qval_tau<0.1&phospho.gene.heat.sig$FC_tau>1],
  `Down in tau flies`=phospho.gene.heat.sig$`Fly Symbol`[
    phospho.gene.heat.sig$qval_tau<0.1&phospho.gene.heat.sig$FC_tau<1])

metab.up.list<-list(`Up in amyloid beta flies`=metab.sig$Compound[
  metab.sig$qval_abeta<0.1&metab.sig$FC_abeta>1],
  `Down in amyloid beta flies`=metab.sig$Compound[
    metab.sig$qval_abeta<0.1&metab.sig$FC_abeta<1],
  `Up in tau flies`=metab.sig$Compound[
    metab.sig$qval_tau<0.1&metab.sig$FC_tau>1],
  `Down in tau flies`=metab.sig$Compound[
    metab.sig$qval_tau<0.1&metab.sig$FC_tau<1])

prot.comb.size<-comb_size(make_comb_mat(prot.up.list))
phos.comb.size<-comb_size(make_comb_mat(phos.up.list))
metab.comb.size<-comb_size(make_comb_mat(metab.up.list))

pdf("drosophila_multiomics/plots/Figure_3C.pdf",
    width=10,height=7.5)
p1<-UpSet(make_comb_mat(prot.up.list), column_title = "Drosophila proteomics",
          comb_order = order(-prot.comb.size),
          row_names_gp = gpar(fontsize = 12, col = "black"),
          column_names_gp = gpar(fontsize = 12, col="black"),
          top_annotation = upset_top_annotation(make_comb_mat(prot.up.list),
                                                axis_param = list(gp = gpar(fontsize = 12))),
          right_annotation = upset_right_annotation(make_comb_mat(prot.up.list),
                                                    axis_param = list(gp = gpar(fontsize = 12)))
) +
  UpSet(make_comb_mat(phos.up.list), column_title = "Drosophila phosphoproteomics",
        comb_order = order(-phos.comb.size),
        row_names_gp = gpar(fontsize = 12, col = "black"),
        top_annotation = upset_top_annotation(make_comb_mat(phos.up.list),
                                              axis_param = list(gp = gpar(fontsize = 12)),
                                              show_annotation_name = FALSE),
        right_annotation = upset_right_annotation(make_comb_mat(phos.up.list),
                                                  axis_param = list(gp = gpar(fontsize = 12))),
        column_names_gp = gpar(fontsize = 12, col="black")) +
  UpSet(make_comb_mat(metab.up.list), column_title = "Drosophila metabolomics",
        comb_order = order(-metab.comb.size),
        top_annotation = upset_top_annotation(make_comb_mat(metab.up.list),
                                              axis_param = list(gp = gpar(fontsize = 12)),
                                              show_annotation_name = FALSE),
        right_annotation = upset_right_annotation(make_comb_mat(metab.up.list),
                                                  axis_param = list(gp = gpar(fontsize = 12))),
        row_names_gp = gpar(fontsize = 12, col = "black"),
        column_names_gp = gpar(fontsize = 12, col="black"))
p1
dev.off()

## heatmap of proteomics and phosphoproteomics again, but by fold change and 
## stars for the significance
prot.fc.nd.mat<-data.matrix(prot.nd.orth[!duplicated(prot.nd.orth$`Fly Symbol`),
                                         c("FC_abeta","FC_tau")])
rownames(prot.fc.nd.mat)<-prot.nd.orth$`Fly Symbol`[!duplicated(prot.nd.orth$`Fly Symbol`)]
colnames(prot.fc.nd.mat)<-c("Amyloid beta","Tau")

prot.qv.nd.mat<-data.matrix(prot.nd.orth[!duplicated(prot.nd.orth$`Fly Symbol`),
                                         c("qval_abeta","qval_tau")])
rownames(prot.qv.nd.mat)<-prot.nd.orth$`Fly Symbol`[!duplicated(prot.nd.orth$`Fly Symbol`)]
colnames(prot.qv.nd.mat)<-c("Amyloid beta","Tau")

pdf("Figure_3D.pdf",
    width=9,height=5.56)
p1<-Heatmap(t(prot.fc.nd.mat),
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(t(prot.qv.nd.mat)[i, j] < 0.1) {
                grid.text(paste0(sprintf("%0.2f",t(prot.fc.nd.mat)[i,j]),"*"),x, y,
                          rot=0)
              } else if(t(prot.qv.nd.mat)[i, j] < 0.01) {
                grid.text(paste0(sprintf("%0.2f",t(prot.fc.nd.mat)[i,j]),"**"),x, y,
                          rot=0)
              }else{
                grid.text(sprintf("%0.2f",t(prot.fc.nd.mat)[i,j]),x,y,
                          rot=0)
              }},
            name="Log2 fold change",
            #col=f2,
            show_row_names = TRUE,
            show_column_names=TRUE,
            row_names_gp=grid::gpar(fontsize=20),
            column_names_gp=grid::gpar(fontsize=20),
            column_names_rot = 60,
            cluster_rows = FALSE,
            cluster_columns = TRUE)
p1

dev.off()

## and now phos
phos.fc.nd.mat<-data.matrix(phospho.gene.nd.orth[,
                                                 c("FC_abeta","FC_tau")])
rownames(phos.fc.nd.mat)<-phospho.gene.nd.orth$phospho_label
colnames(phos.fc.nd.mat)<-c("Amyloid beta","Tau")

phos.qv.nd.mat<-data.matrix(phospho.gene.nd.orth[,
                                                 c("qval_abeta","qval_tau")])
rownames(phos.qv.nd.mat)<-phospho.gene.nd.orth$phospho_label
colnames(phos.qv.nd.mat)<-c("Amyloid beta","Tau")

pdf("Figure_3E.pdf",
    width=9,height=5.56)
p1<-Heatmap(t(phos.fc.nd.mat),
            cell_fun = function(j, i, x, y, w, h, fill) {
              if(t(phos.qv.nd.mat)[i, j] < 0.1) {
                grid.text(paste0(sprintf("%0.2f",t(phos.fc.nd.mat)[i,j]),"*"),x, y,
                          rot=0)
              } else if(t(phos.qv.nd.mat)[i, j] < 0.01) {
                grid.text(paste0(sprintf("%0.2f",t(phos.fc.nd.mat)[i,j]),"**"),x, y,
                          rot=0)
              }else{
                grid.text(sprintf("%0.2f",t(phos.fc.nd.mat)[i,j]),x,y,
                          rot=0)
              }},
            name="Log2 fold change",
            #col=f2,
            show_row_names = TRUE,
            show_column_names=TRUE,
            row_names_gp=grid::gpar(fontsize=20),
            column_names_gp=grid::gpar(fontsize=20),
            column_names_rot = 60,
            cluster_rows = FALSE,
            cluster_columns = TRUE)
p1

dev.off()

## for proteomics perform pathway enrichments for the + tau -abeta, +tau +abeta,
## -tau -abeta and -tau + abeta groups
sig.list.prot<-list(prot.sig[prot.sig$FC_tau>0&prot.sig$FC_abeta>0&
                               prot.sig$qval_abeta<0.1&
                               prot.sig$qval_tau<0.1],
                    prot.sig[prot.sig$FC_tau>0&prot.sig$FC_abeta<0&
                               prot.sig$qval_abeta<0.1&
                               prot.sig$qval_tau<0.1],
                    prot.sig[prot.sig$FC_tau<0&prot.sig$FC_abeta>0&
                               prot.sig$qval_abeta<0.1&
                               prot.sig$qval_tau<0.1],
                    prot.sig[prot.sig$FC_tau<0&prot.sig$FC_abeta<0&
                               prot.sig$qval_abeta<0.1&
                               prot.sig$qval_tau<0.1])
groups<-c("Up tau, up abeta","Up tau, down abeta",
          "Down tau, up abeta","Down tau, down abeta")
pathway.res<-data.table()

for(i in 1:length(sig.list.prot)){
  gostres <- gost(query = sig.list.prot[[i]]$`Fly Symbol`,
                  organism = "dmelanogaster", ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = TRUE,
                  user_threshold = 0.05,correction="g_SCS",
                  domain_scope = "annotated", 
                  numeric_ns = "", as_short_link = FALSE)
  
  results<-gostres$result
  enrich.k.r<-results
  
  enrich.k.r<-enrich.k.r[order(enrich.k.r$p_value),]
  enrich.k.r<-enrich.k.r[enrich.k.r$term_size<500&
                           enrich.k.r$intersection_size>3,]
  print(enrich.k.r[enrich.k.r$source=="GO:BP",c(
    "term_name","term_size","intersection")][1:10,])
  print(enrich.k.r[enrich.k.r$source=="REAC",c(
    "term_name","term_size","intersection")][1:10,])
  
  supptable.enrich<-enrich.k.r[grep("GO:BP",enrich.k.r$source),]
  if(nrow(supptable.enrich)>0){
    supptable.enrich$overlapping_genes<-unlist(supptable.enrich$intersection)
    supptable.enrich$group<-groups[i]
    pathway.res<-rbind(pathway.res,supptable.enrich)
  }
}

pathway.res$logp<-(-log10(pathway.res$p_value))
top.10.by.group<-pathway.res[order(pathway.res),head(.SD,10),by="group"]
top.10.by.group$logp<-(-log10(top.10.by.group$p_value))

pathway.res %>%
  group_by(group) %>%
  top_n(10) %>%
  ungroup %>%
  mutate(group = as.factor(group),
         term_name = reorder_within(term_name, logp, group)) %>%
  ggplot(aes(y=term_name, x=logp
  )) + 
  geom_col()+
  facet_wrap(~group,scales="free")+
  scale_y_reordered() +
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text=element_text(size=16))+
  xlab('Negative log10 FDR-adjusted p-value')+
  ylab('enriched pathways')
ggsave("Supplementary_figure_4G.pdf",
       width=20,height=12.5)

## check enrichment for those nd modifiers diff abundant in proteomics
gostres <- gost(query = prot.nd.orth$`Fly Symbol`[prot.nd.orth$qval_abeta<0.1],
                organism = "dmelanogaster", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = TRUE,
                user_threshold = 0.05,correction="g_SCS",
                domain_scope = "annotated", 
                numeric_ns = "", as_short_link = FALSE)

results<-gostres$result

enrich.k.r<-results[grep("REAC|GO",results$source)
                    ,]

enrich.k.r<-enrich.k.r[order(enrich.k.r$p_value),]
enrich.k.r<-enrich.k.r[enrich.k.r$term_size<500&
                         enrich.k.r$intersection_size>3,]
print(enrich.k.r[enrich.k.r$source=="GO:BP",c(
  "term_name","term_size","intersection")][1:10,])
enrich.k.r$logp<-(-log10(enrich.k.r$p_value))

ggplot(enrich.k.r,aes(y=reorder(term_name,
                                logp), x=logp
)) + 
  geom_col()+
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text=element_text(size=16))+
  xlab('Negative log10 FDR-adjusted p-value')+
  ylab('enriched pathways')
ggsave("Supplementary_figure_4H.pdf",
       width=20,height=12.5)

## make volcano plots
phospho.gene$phospho_label<-paste0(phospho.gene$`Fly Symbol`,
                                   "_ST",phospho.gene$site_position_1)
phospho.gene$phospho_label[
  !is.na(phospho.gene$site_position_2)]<-paste0(
    phospho.gene$`Fly Symbol`[!is.na(phospho.gene$site_position_2)],
    "_ST",phospho.gene$site_position_1[!is.na(phospho.gene$site_position_2)],
    "/ST",phospho.gene$site_position_2[!is.na(phospho.gene$site_position_2)])

phospho.gene$phospho_label[
  !is.na(phospho.gene$site_position3)]<-paste0(
    phospho.gene$`Fly Symbol`[!is.na(phospho.gene$site_position3)],
    "_ST",phospho.gene$site_position_1[!is.na(phospho.gene$site_position3)],
    "/ST",phospho.gene$site_position_2[!is.na(phospho.gene$site_position3)],
    "/ST",phospho.gene$site_position3[!is.na(phospho.gene$site_position3)])

df.list<-list(human.dro,phospho.gene,metab)

for(i in 1:length(df.list)){
  test.df<-df.list[[i]]
  
  if (i==1){
    ## tau plot
    test.df$label<-"no significant change"
    test.df$label[test.df$qval_tau<0.1&
                    test.df$FC_tau<0]<-"significantly downregulated"
    test.df$label[test.df$qval_tau<0.1&
                    test.df$FC_tau>0]<-"significantly upregulated"
    
    ggplot(test.df, aes(y=-log10(qval_tau), x=FC_tau
    )) + 
      geom_point(aes(color=label))+
      theme_classic()+
      scale_color_manual(values=c("no significant change"="black",
                                  "significantly downregulated"="blue",
                                  "significantly upregulated"="red"))+
      theme(axis.text=element_text(size=16,color="black"),
            axis.title=element_text(size=16),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16))+
      # geom_text_repel(data=test.df[test.df$qval_tau<0.1],aes(label=`Fly Symbol`),
      #                 size=5,max.overlaps = 10)+
      #guides(fill = guide_legend(override.aes = aes(text = "")))+
      geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
      ylab('Negative log10 FDR-adjusted p-value')+
      xlab('log2 fold change')
    ggsave("Supplementary_figure_4B.pdf",
           width=7,height=5)
    
    ## abeta plot
    test.df$label<-"no significant change"
    test.df$label[test.df$qval_abeta<0.1&
                    test.df$FC_abeta<0]<-"significantly downregulated"
    test.df$label[test.df$qval_abeta<0.1&
                    test.df$FC_abeta>0]<-"significantly upregulated"
    
    ggplot(test.df, aes(y=-log10(qval_abeta), x=FC_abeta
    )) + 
      geom_point(aes(color=label))+
      theme_classic()+
      scale_color_manual(values=c("no significant change"="black",
                                  "significantly downregulated"="blue",
                                  "significantly upregulated"="red"))+
      theme(axis.text=element_text(size=16,color="black"),
            axis.title=element_text(size=16),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16))+
      # geom_text_repel(data=test.df[test.df$qval_abeta<0.1],aes(label=`Fly Symbol`),
      #                 size=5,max.overlaps = 10)+
      #guides(fill = guide_legend(override.aes = aes(text = "")))+
      geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
      ylab('Negative log10 FDR-adjusted p-value')+
      xlab('log2 fold change')
    ggsave("Supplementary_figure_4A.pdf",
           width=7,height=5)
    
  }else if(i==2){
    test.df$label<-"no significant change"
    test.df$label[test.df$qval_tau<0.1&
                    test.df$FC_tau<0]<-"significantly downregulated"
    test.df$label[test.df$qval_tau<0.1&
                    test.df$FC_tau>0]<-"significantly upregulated"
    
    ggplot(test.df, aes(y=-log10(qval_tau), x=FC_tau
    )) + 
      geom_point(aes(color=label))+
      theme_classic()+
      scale_color_manual(values=c("no significant change"="black",
                                  "significantly downregulated"="blue",
                                  "significantly upregulated"="red"))+
      theme(axis.text=element_text(size=16,color="black"),
            axis.title=element_text(size=16),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16))+
      # geom_text_repel(data=test.df[test.df$qval_tau<0.1],aes(label=phospho_label),
      #                 size=5,max.overlaps = 10)+
      #guides(fill = guide_legend(override.aes = aes(text = "")))+
      geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
      ylab('Negative log10 FDR-adjusted p-value')+
      xlab('log2 fold change')
    ggsave("Supplementary_Figure_4D.pdf",
           width=7,height=5)
    
    ## abeta plot
    test.df$label<-"no significant change"
    test.df$label[test.df$qval_abeta<0.1&
                    test.df$FC_abeta<0]<-"significantly downregulated"
    test.df$label[test.df$qval_abeta<0.1&
                    test.df$FC_abeta>0]<-"significantly upregulated"
    
    ggplot(test.df, aes(y=-log10(qval_abeta), x=FC_abeta
    )) + 
      geom_point(aes(color=label))+
      theme_classic()+
      scale_color_manual(values=c("no significant change"="black",
                                  "significantly downregulated"="blue",
                                  "significantly upregulated"="red"))+
      theme(axis.text=element_text(size=16,color="black"),
            axis.title=element_text(size=16),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16))+
      # geom_text_repel(data=test.df[qval_abeta<0.1],aes(label=phospho_label),
      #                 size=5,max.overlaps = 10)+
      #guides(fill = guide_legend(override.aes = aes(text = "")))+
      geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
      ylab('Negative log10 FDR-adjusted p-value')+
      xlab('log2 fold change')
    ggsave("Supplementary_Figure_4C.pdf",
           width=7,height=5)
    
  }else if(i==3){
    test.df$label<-"no significant change"
    test.df$label[test.df$qval_tau<0.1&
                    test.df$FC_tau<0]<-"significantly downregulated"
    test.df$label[test.df$qval_tau<0.1&
                    test.df$FC_tau>0]<-"significantly upregulated"
    
    ggplot(test.df, aes(y=-log10(qval_tau), x=FC_tau
    )) + 
      geom_point(aes(color=label))+
      theme_classic()+
      scale_color_manual(values=c("no significant change"="black",
                                  "significantly downregulated"="blue",
                                  "significantly upregulated"="red"))+
      theme(axis.text=element_text(size=16,color="black"),
            axis.title=element_text(size=16),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16))+
      # geom_text_repel(data=test.df[test.df$qval_tau],aes(label=Metabolite),
      #                 size=5,max.overlaps = 10)+
      #guides(fill = guide_legend(override.aes = aes(text = "")))+
      geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
      ylab('Negative log10 FDR-adjusted p-value')+
      xlab('log2 fold change')
    ggsave("Supplementary_Figure_4F.pdf",
           width=7,height=5)
    
    ## abeta plot
    test.df$label<-"no significant change"
    test.df$label[test.df$qval_abeta<0.1&
                    test.df$FC_abeta<0]<-"significantly downregulated"
    test.df$label[test.df$qval_abeta<0.1&
                    test.df$FC_abeta>0]<-"significantly upregulated"
    
    ggplot(test.df, aes(y=-log10(qval_abeta), x=FC_abeta
    )) + 
      geom_point(aes(color=label))+
      theme_classic()+
      scale_color_manual(values=c("no significant change"="black",
                                  "significantly downregulated"="blue",
                                  "significantly upregulated"="red"))+
      theme(axis.text=element_text(size=16,color="black"),
            axis.title=element_text(size=16),
            legend.text=element_text(size=16),
            legend.title=element_text(size=16))+
      #geom_text_repel(data=test.df[test.df$qval_abeta<0.1],aes(label=Metabolite),
      #                size=5,max.overlaps = 10)+
      #guides(fill = guide_legend(override.aes = aes(text = "")))+
      geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
      ylab('Negative log10 FDR-adjusted p-value')+
      xlab('log2 fold change')
    ggsave("Supplementary_Figure_4E.pdf",
           width=7,height=5)
  }
}