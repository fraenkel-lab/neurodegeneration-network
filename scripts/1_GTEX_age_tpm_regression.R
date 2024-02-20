library(data.table)
library(mgcv)
library(ggplot2)
library(lme4)
library(lmerTest)
library(enrichplot)

#run hypergeometric test
hypergeom_test<-function(x,m,n,k){
  # x: number of successes
  # m: number of successes in population
  # n: size of population
  # k: number of draws
  return(1-phyper(x-1,m,n-m,k))
}

mixed_model_func<-function(d){
  m <- lmer(z ~age_death+Sex+PMI+tissue_origin+Latino+(1|ID),
            data=d,REML=FALSE)
  
  p<-c(summary(m)[10]$coefficients["age_death","Pr(>|t|)"],
       summary(m)[10]$coefficients["age_death","Estimate"],
       summary(m)[10]$coefficients["age_death","t value"])
  return(p)
}

gam_func<-function(d){
  m <- gam(z ~ age_death+Sex+PMI+Latino,data=d)
  p <- c(summary(m)[[4]][[2]],summary(m)[[1]][[2]],summary(m)[[3]][[2]])
  return(p)
}

# download the tpm data for all brain tissues in the GTEx project, found at
# https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
# for GTEx V8
direc<-"GTEX_brain_tx/tpm/"
files<-dir(direc,recursive=TRUE,pattern="tpm")

gtex.expression<-fread(paste0(direc,files[1]))

#uniprot key for limiting to just protein-coding genes
key<-fread("uniprot_key.txt")
key<-key[!duplicated(key$`UniProtKB Gene Name symbol`)&
           !key$`UniProtKB Gene Name symbol`%in%c("")]

for(i in 1:length(files)){
  tmp<-fread(paste0(direc,files[i]))
  
  tmp.lim<-tmp[,-which(colnames(tmp)%in%colnames(gtex.expression)),
               with=FALSE]
  
  gtex.expression<-cbind(gtex.expression,tmp.lim)
}

# write concatenated file for all brain tissues
fwrite(gtex.expression,"GTEX_brain_tx/tpm/all_brain_cpm.txt",sep='\t',
       quote=F,row.names=F)

gtex.expression<-fread("GTEX_brain_tx/tpm/all_brain_cpm.txt")

#closed-access phenotyping data from the GTEx project
gtex.metadat<-fread("GTEX_closed_meta/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")
gtex.sample.meta<-fread("GTEX_closed_meta/SampleAttributesDS.txt")

tissue.ids<-unique(gtex.sample.meta$SMTSD[gtex.sample.meta$SMTS=="Brain"])
tissue.ids<-c(tissue.ids,"all brain")

# perform analysis for each tissue, supplementary figure
for(tissue in tissue.ids){
  gtex.sample.meta<-fread("GTEX_closed_meta/SampleAttributesDS.txt")
  
  if(tissue%in%c("all brain")){
    gtex.sample.meta<-gtex.sample.meta[gtex.sample.meta$SMTS=="Brain"
    ]
    
    subj.dict<-fread("GTEX_closed_meta/phs000424.v8.pht002741.v8.p2.GTEx_Sample.MULTI.txt")
    
    gtex.metadat.merge<-merge(gtex.metadat,subj.dict[
      subj.dict$SAMPID%in%colnames(gtex.expression),c("SUBJID","SAMPID")],
      by="SUBJID")
    
    gtex.metadat.sample<-merge(gtex.metadat.merge,gtex.sample.meta,
                               by="SAMPID")
    
    gtex.expression.unique<-gtex.expression[,!duplicated(
      colnames(gtex.expression))&
        colnames(gtex.expression)%in%
        gtex.metadat.sample$SAMPID,
      with=FALSE]
    
    gtex.metadat.match<-gtex.metadat.sample[match(colnames(gtex.expression.unique),
                                                  gtex.metadat.sample$SAMPID),]
    
    gtex.mat<-data.table(t(gtex.expression.unique))
    colnames(gtex.mat)<-gtex.expression$Description
    
    gtex.mat[is.na(gtex.mat)]<-1e-6
    
    #for all tissues, use mixed effects model
    gam.p <- apply(gtex.mat, 2, function(z){
      #input dataframe
      d <- data.frame(z=z, age_death=gtex.metadat.match$AGE,
                      Sex=gtex.metadat.match$SEX,
                      PMI=gtex.metadat.match$SMTSISCH/60,
                      tissue_origin=factor(gtex.metadat.match$SMTSD),
                      Latino=factor(gtex.metadat.match$ETHNCTY),
                      ID=factor(gtex.metadat.match$SUBJID),
                      stringsAsFactors=FALSE)
      
      p <- tryCatch(
        mixed_model_func(d),error=function(e){
          data.frame(pval=1,Estimate=0,t_value=0)
        })
      p
      
    })
    
    de.table<-data.table()
    cor.table<-data.table()
    for(i in 1:(ncol(gam.p))){
      
      #mixed-effects model outputs
      temp.table<-data.table(Gene=colnames(gam.p)[i],
                             p_val=unlist(gam.p[1,i]),
                             coefficient=unlist(gam.p[2,i]),
                             t_statistic=unlist(gam.p[3,i]),
                             stringsAsFactors=F)
      
      de.table<-rbind(temp.table,de.table)
    }
  }
  gtex.sample.meta<-gtex.sample.meta[gtex.sample.meta$SMTS=="Brain"&
                                     gtex.sample.meta$SMTSD==tissue
  ]
  
  subj.dict<-fread("GTEX_closed_meta/phs000424.v8.pht002741.v8.p2.GTEx_Sample.MULTI.txt")
  
  gtex.metadat.merge<-merge(gtex.metadat,subj.dict[
    subj.dict$SAMPID%in%colnames(gtex.expression),c("SUBJID","SAMPID")],
    by="SUBJID")
  
  gtex.metadat.sample<-merge(gtex.metadat.merge,gtex.sample.meta,
                             by="SAMPID")
  
  gtex.expression.unique<-gtex.expression[,!duplicated(
    colnames(gtex.expression))&
      colnames(gtex.expression)%in%
      gtex.metadat.sample$SAMPID,
    with=FALSE]
  
  gtex.metadat.match<-gtex.metadat.sample[match(colnames(gtex.expression.unique),
                                                gtex.metadat.sample$SAMPID),]
  
  gtex.mat<-data.table(t(gtex.expression.unique))
  colnames(gtex.mat)<-gtex.expression$Description
  
  gtex.mat[is.na(gtex.mat)]<-1e-6
  
  #generalized additive model by tissue
  gam.p <- apply(gtex.mat, 2, function(z){
    d <- data.frame(z=z, age_death=gtex.metadat.match$AGE,
                    Sex=gtex.metadat.match$SEX,
                    PMI=gtex.metadat.match$SMTSISCH/60,
                    tissue_origin=factor(gtex.metadat.match$SMTSD),
                    Latino=factor(gtex.metadat.match$ETHNCTY),
                    ID=factor(gtex.metadat.match$SUBJID),
                    stringsAsFactors=FALSE)
    
    p <- gam_func(d)
    p
                            
                          })
  
  de.table<-data.table()
  cor.table<-data.table()
  for(i in 1:(ncol(gam.p))){
    
    
    #GAM outputs
    temp.table<-data.table(Gene=colnames(gam.p)[i],
                           p_val=gam.p[1,i],
                           coefficient=gam.p[2,i],
                           t_statistic=gam.p[3,i],
                           stringsAsFactors=F)
    
    de.table<-rbind(temp.table,de.table)
  }
  
  de.table.sub<-de.table[order(de.table$p_val),]
  
  de.table.sub<-de.table.sub[!is.na(de.table.sub$p_val)&
                               !duplicated(de.table.sub$Gene)&
                               de.table.sub$Gene%in%key$`Gene name`
                             ,]
  
  de.table.sub$padj<-p.adjust(de.table.sub$p_val,method="fdr")
  de.table.sub$eff_size<-de.table.sub$t_stat*sqrt(var(gtex.metadat.match$AGE))
  
  write.table(de.table.sub,
              paste0("GTEX_brain_tx/outs/",
                     tissue,"_age_genes_tpm_prot_coding_PMI_mixed.txt"),
              sep='\t',quote=F,row.names=F)
}

## heatmap of which nd modifiers are associated with aging by tissue type
gene.table<-data.table(Gene=gtex.expression$Description)

direc<-"GTEX_brain_tx_outs/"
files<-dir(direc,recursive=FALSE,pattern="_mixed.txt")
tissue.ids<-gsub("Brain - ([A-Za-z ()]+)_.*","\\1",files)
tissue.ids[c(2,7,10,12)]<-c("Combined Brain Tissues","Anterior Cingulate Cortex",
                              "Frontal Cortex","Nucleus Accumbens","Spinal Cord")

for(i in 2:length(files)){
  temp.tissue<-fread(paste0(direc,files[i]))
  
  if(is.null(temp.tissue$coefficient)){
    temp.tissue$coefficient<-temp.tissue$t_stat
  }
  
  # define age associated effect size as a coefficient greater than 0.1,
  # indicating a change of 1TPM/decade
  temp.nd<-temp.tissue[temp.tissue$padj<1e-1&abs(temp.tissue$coefficient)>0.1
                       ,
                       c("Gene","coefficient")]
  
  colnames(temp.nd)[2]<-tissue.ids[i]
  gene.table<-merge(gene.table,temp.nd,by="Gene",all.x=TRUE)
}

nd.mat<-as.matrix(gene.table[,2:ncol(gene.table),with=FALSE])
nd.mat[is.na(nd.mat)]<-0

nd.mat<-nd.mat[,!colSums(nd.mat)==0]

rownames(nd.mat)<-gene.table$Gene

#nd.mat<-nd.mat[rownames(nd.mat)%in%nd.mod$`human gene`,]

breaks<-c(seq(min(nd.mat,na.rm=TRUE),0,length.out=ceiling(5001/2)+1),
          seq(max(nd.mat,na.rm=TRUE)/5001,max(nd.mat,na.rm=TRUE),
              length.out=floor(5001/2)))

# breaks<-c(seq(min(nd.mat,na.rm=TRUE),0,length.out=500))
cols<-colorRampPalette(c("blue","white","red"))(5001)

n <- 26
qual_col_pals<-brewer.pal.info[brewer.pal.info$category == 'qual'&
                                 brewer.pal.info$colorblind==TRUE,]
col_vector<-unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

nd.mat.plot<-nd.mat[which(rowSums(nd.mat)!=0)#&
                    #row.anno=="age-associated splicing"
                    , ]
heatmap_fig<-pheatmap(nd.mat.plot,
                      cluster_rows = TRUE,
                      show_rownames=FALSE,
                      breaks=breaks,
                      #annotation_col=col.annotation,
                      #annotation_colors = anno.cols,
                      #annotation_row=row.anno,
                      #annotation_col=col.anno,
                      scale="none",
                      color=cols,
                      fontsize_row = 10, fontsize_col = 16,
                      show_colnames=TRUE,cluster_cols=TRUE)

save_pheatmap_png <- function(x, filename, width=2000, height=1750, res = 150) {
  #png(filename, width = width, height = height, res = res)
  pdf(filename)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(heatmap_fig,
                  "Supplementary_Figure_1.pdf")

## after files are written re-read regression results
de.table.sub<-fread(
  "GTEX_brain_tx/outs/all_brain_age_genes_tpm_prot_coding_PMI_mixed.txt")

# for supplementary figure 5 plot regression of HNRNPA2B1 expression with age
# in frontal cortex
gene.set<-c("HNRNPA2B1")
plot.df<-data.table()
gam.list<-data.table()

for(gene in gene.set){
  test.plot.df<-data.frame(Gene_expression=unlist(gtex.mat[,gene,
                                                           with=FALSE]), 
                           Age=gtex.metadat.match$AGE,
                           Sex=gtex.metadat.match$SEX,
                           PMI=gtex.metadat.match$SMTSISCH/60,
                           tissue_origin=factor(gtex.metadat.match$SMTSD),
                           ethnicity=factor(gtex.metadat.match$ETHNCTY),
                           sample=factor(gtex.metadat.match$SUBJID),
                           stringsAsFactors=FALSE)
  
  gam.out<-lmer(Gene_expression~Age+Sex+PMI+tissue_origin+
                  ethnicity+(1|sample),data=test.plot.df,REML=TRUE)
  
  gam.summary<-summary(gam.out)
  
  temp.dt<-data.table(Gene=gene,
                      b=gam.summary$coefficients["(Intercept)","Estimate"],
                      coef=gam.summary$coefficients["Age","Estimate"])
  
  test.plot.df$Gene<-gene
  
  test.plot.df<-merge(temp.dt,test.plot.df,by="Gene")
  plot.df<-rbind(test.plot.df,plot.df)
  
}


ggplot(plot.df,aes(x=Age,y=Gene_expression))+
  geom_point(aes(x=Age,y=Gene_expression),alpha=0.2)+
  geom_abline(aes(intercept=b,slope=coef),
              size=1,color="red")+
  ylab('Gene Expression (TPM)')+
  theme_classic()+
  theme(axis.text=element_text(size=24),
        axis.text.x=element_text(size=24),
        axis.text.y=element_text(size=24),
        axis.title=element_text(size=24),
        strip.text = element_text(size=24))+
  facet_wrap(~Gene,scales="free_y")
ggsave("Supplementary_Figure_5.pdf")

## plot mean expression of neurodegeneration modifiers for figure 2B

# first re-read the expression for all brain tissue
gtex.sample.meta<-fread("GTEX_closed_meta/SampleAttributesDS.txt")

gtex.sample.meta<-gtex.sample.meta[gtex.sample.meta$SMTS=="Brain"
]

subj.dict<-fread("GTEX_closed_meta/phs000424.v8.pht002741.v8.p2.GTEx_Sample.MULTI.txt")

gtex.metadat.merge<-merge(gtex.metadat,subj.dict[
  subj.dict$SAMPID%in%colnames(gtex.expression),c("SUBJID","SAMPID")],
  by="SUBJID")

gtex.metadat.sample<-merge(gtex.metadat.merge,gtex.sample.meta,
                           by="SAMPID")

gtex.expression.unique<-gtex.expression[,!duplicated(
  colnames(gtex.expression))&
    colnames(gtex.expression)%in%
    gtex.metadat.sample$SAMPID,
  with=FALSE]

gtex.metadat.match<-gtex.metadat.sample[match(colnames(gtex.expression.unique),
                                              gtex.metadat.sample$SAMPID),]

gtex.mat<-data.table(t(gtex.expression.unique))
colnames(gtex.mat)<-gtex.expression$Description

gtex.mat[is.na(gtex.mat)]<-1e-6
#calculate geometric means
log.exp<-data.matrix(log(gtex.mat[,
                                  colnames(gtex.mat)%in%nd.mod$`human gene`,
                                  with=FALSE]))

log.exp[!is.finite(log.exp)] <- 0
log.exp.t<-data.table(log.exp)
log.exp.t$patid<-unlist(gtex.metadat.match$SUBJID)
log.exp.sampmean<-log.exp.t[,lapply(.SD,mean),by="patid",
                            .SDcols=colnames(log.exp.t)[1:(NCOL(log.exp.t)-1)]]
log.exp.all<-log(gtex.mat)

log.exp.all<-data.matrix(log.exp.all)
log.exp.all[!is.finite(log.exp.all)] <- 0
log.exp.all.t<-data.table(log.exp.all)
log.exp.all.t$patid<-unlist(gtex.metadat.match$SUBJID)
log.exp.all.sampmean<-log.exp.t[,lapply(.SD,mean),by="patid",
                                .SDcols=colnames(log.exp.all.t)[
                                  1:(NCOL(log.exp.all.t)-1)]]

plot.df<-data.table(Gene_expression=exp(rowMeans(log.exp.samp[,2:ncol(log.exp)])), 
                    Age=gtex.metadat.match$AGE,
                    Sex=gtex.metadat.match$SEX,
                    PMI=gtex.metadat.match$SMTSISCH/60,
                    tissue_origin=factor(gtex.metadat.match$SMTSD),
                    ethnicity=factor(gtex.metadat.match$ETHNCTY),
                    Sample=factor(gtex.metadat.match$SUBJID),
                    stringsAsFactors=FALSE)

gam.out<-lmer(Gene_expression~Age+Sex+PMI+tissue_origin+
                ethnicity+(1|Sample),data=plot.df,REML=FALSE)

### for GAM
gam.out<-gam(Gene_expression~Age+Sex+PMI,data=plot.df,REML=FALSE)

gam.summary<-summary(gam.out)

plot.df.all<-data.table(Gene_expression=exp(rowMeans(log.exp.all)), 
                        Age=gtex.metadat.match$AGE,
                        Sex=gtex.metadat.match$SEX,
                        PMI=gtex.metadat.match$SMTSISCH/60,
                        tissue_origin=factor(gtex.metadat.match$SMTSD),
                        ethnicity=factor(gtex.metadat.match$ETHNCTY),
                        Sample=factor(gtex.metadat.match$SUBJID),
                        stringsAsFactors=FALSE)

#for mixed effect model
gam.out.all<-lmer(Gene_expression~Age+Sex+PMI+tissue_origin+
                    ethnicity+(1|Sample),data=plot.df.all,REML=FALSE)

gam.summary.all<-summary(gam.out.all)

# get confidence interval (seemingly standard error of mean)
gam.summary.all.ci<-merTools::predictInterval(gam.out.all, as.data.frame(plot.df))

gam.summary.ci<-merTools::predictInterval(gam.out, as.data.frame(plot.df))

age.gam.summary<-as.data.frame(effects::effect("Age",gam.out))
age.gam.summary$group<-"Neurodegeneration Modifiers"
age.gam.all.summary<-as.data.frame(effects::effect("Age",gam.out.all))
age.gam.all.summary$group<-"All Genes"

full.plot.df<-rbind(plot.df,plot.df.all)

ggplot(full.plot.df,aes(x=Age,y=Gene_expression))+
  geom_point(alpha=0.2)+
  geom_ribbon(data=age.gam.summary,
              aes(ymin = lower, ymax = upper, y=fit), alpha = .15) +
  geom_ribbon(data=age.gam.all.summary,
              aes(ymin = lower, ymax = upper,y=fit), alpha = .15) +
  geom_line(data=age.gam.summary,aes(Age,fit,fill=group,
                                     group=group),
            size=1,color="#FF7F0E",show.legend=TRUE)+
  geom_line(data=age.gam.all.summary,aes(Age,fit,fill=group,
                                         group=group),
            size=1,color="#1F77B4",show.legend=TRUE)+
  scale_colour_manual("Gene Set",
                      breaks = c("Nerodegeneration Genes","All Genes"),
                      values = c("#1F77B4", "#FF7F0E"))+
  facet_wrap(~group,scales="free_y")+
  # scale_fill_manual("Gene Set",
  #                     breaks = c("Nerodegeneration Modifiers","All Genes"),
  #                     values = c("#1F77B4", "#FF7F0E"))+
  ylab('Gene Expression (TPM)')+
  xlab('Age (Years)')+
  theme_classic()+
  theme(axis.text=element_text(size=24,color="black"),
        axis.text.x=element_text(size=24,color="black"),
        axis.text.y=element_text(size=24,color="black"),
        axis.title=element_text(size=24,color="black"),
        strip.text = element_text(size=24,color="black"),
        legend.text = element_text(size=20,color="black"),
        legend.title = element_text(size=22,color="black"))
ggsave("Figure_2B.pdf",
       height=12)

#permutation test of signficance of ND mod gene set
n_perm<-10000

test.stat.p<-gam.summary$coefficients["Age","Pr(>|t|)"]
test.stat.coef<-gam.summary$coefficients["Age","Estimate"]
stat.df<-data.table()
for(i in 1:n_perm){
  #gene.set<-sample(key$`Gene name`,length(unique(nd.mod$`human gene`)))
  
  random.idx<-sample(seq(1,nrow(gtex.metadat.match)),nrow(gtex.metadat.match))
  plot.df<-data.table(Gene_expression=unlist(rowMeans(gtex.mat[,
                                                               colnames(gtex.mat)%in%nd.mod$`human gene`,
                                                               with=FALSE])), 
                      Age=gtex.metadat.match$AGE[random.idx],
                      Sex=gtex.metadat.match$SEX,
                      PMI=gtex.metadat.match$SMTSISCH/60,
                      tissue_origin=factor(gtex.metadat.match$SMTSD),
                      ethnicity=factor(gtex.metadat.match$ETHNCTY),
                      Sample=factor(gtex.metadat.match$SUBJID),
                      stringsAsFactors=FALSE)
  
  tmp.gam<-lmer(Gene_expression~Age+Sex+PMI+tissue_origin+
                  ethnicity+(1|Sample),data=plot.df,REML=FALSE)
  
  tmp.gam.summary<-summary(tmp.gam)
  temp.stat<-data.table(t_value=tmp.gam.summary$coefficients["Age","Estimate"],
                        p_value=tmp.gam.summary$coefficients["Age","Pr(>|t|)"],stringsAsFactors=F)
  
  stat.df<-rbind(temp.stat,stat.df)
}

t.empirical<-length(which(abs(stat.df$t_value)>=abs(test.stat.coef)))
p.empirical<-length(which(stat.df$p_value<=test.stat.p))

## generate GSEA results for all brain expression (table 1), Figure 2C
nd.mod<-fread("Table1.txt")
nd.mod<-nd.mod[!is.na(nd.mod$`human gene`)]
nd.mod<-nd.mod[-which(nd.mod$`human gene`%in%c(""," "))]

gene.set<-nd.mod$`human gene`
nd.df<-data.table(ID="Neurodegeneration Screen Hits",Gene=gene.set)
em2 <- GSEA(ranks, TERM2GENE = nd.df)


gseaplot2(em2, geneSetID = 1, title = "Neurodegeneration Modifiers",
          color="black",base_size=16,pvalue_table=FALSE)
ggsave("Figure_2C")

#plot enrichment of neurodegeneration modifiers in each tissue
direc<-"GTEX_brain_tx/outs/tissue_age_assoc/"
files<-dir(direc,recursive=FALSE,pattern="_mixed.txt")

#elminate putamen since there are no age-associated 
tissue.files<-files[-grep("Putamen",files)]

deg.dt<-data.table()
for(i in 1:ncol(nd.mat.plot)){
  temp.de<-fread(paste0(direc,files[i]))
  if(is.null(temp.de$coefficient)){
    temp.de$coefficient<-temp.de$t_stat
  }
  out.p<-hypergeom_test(length(which(nd.mod$`human gene`%in%
                                       temp.de$Gene[
                                         temp.de$padj<0.1&
                                           abs(temp.de$coefficient)>0.1
                                       ])),
                        nrow(temp.de[!duplicated(temp.de$Gene)&temp.de$padj<1e-1&
                                       abs(temp.de$coefficient)>0.1]),
                        nrow(key),nrow(nd.mod))
  
  temp.pop<-data.table(tissue=colnames(nd.mat)[i],
                       group="All Protein-Coding Genes",
                       DEG_proportion=nrow(temp.de[
                         !duplicated(temp.de$Gene)&temp.de$padj<1e-1&
                           abs(temp.de$coefficient)>0.1])/nrow(key),
                       p_val=out.p,
                       stringsAsFactors=F)
  pop.ci<-binom.test(nrow(temp.de[
    !duplicated(temp.de$Gene)&temp.de$padj<1e-1&abs(temp.de$coefficient)>0.1]),
    nrow(key),temp.pop$DEG_proportion)$conf.int
  
  temp.pop$upper_CI<-pop.ci[2]
  temp.pop$lower_CI<-pop.ci[1]
  
  temp.nd<-data.table(tissue=colnames(nd.mat)[i],
                      group="Neurodegeneration Modifiers",
                      DEG_proportion=length(which(nd.mod$`human gene`%in%
                                                    temp.de$Gene[
                                                      temp.de$padj<0.1&
                                                        abs(temp.de$coefficient)>0.1
                                                    ]))/nrow(nd.mod),
                      p_val=out.p,
                      stringsAsFactors=F)
  
  pop.ci<-binom.test(length(which(nd.mod$`human gene`%in%
                                    temp.de$Gene[
                                      temp.de$padj<0.1&
                                        abs(temp.de$coefficient)>0.1
                                    ])),
                     nrow(nd.mod),temp.nd$DEG_proportion)$conf.int
  
  temp.nd$upper_CI<-pop.ci[2]
  temp.nd$lower_CI<-pop.ci[1]
  
  deg.dt<-rbind(temp.pop,temp.nd,deg.dt)
  print(c(colnames(nd.mat)[i],out.p))
  
}

deg.dt$tissue<-gsub(" \\(basal ganglia\\)","",deg.dt$tissue)
deg.dt<-deg.dt[!deg.dt$tissue%in%c("Spinal Cord","Putamen")]
deg.dt$tissue_factor<-factor(deg.dt$tissue,levels=c("Amygdala",
                                                    "Caudate",
                                                    "Anterior Cingulate Cortex",
                                                    "Nucleus Accumbens",
                                                    "Cerebellar Hemisphere","Cerebellum",
                                                    "Frontal Cortex","Hypothalamus",
                                                    "Combined Brain Tissues","Hippocampus",
                                                    "Cortex","Substantia nigra"))
ggplot(deg.dt,
       aes(y=DEG_proportion, x=tissue_factor,
           fill=group)) + 
  geom_bar(stat="identity",position="dodge")+
  #geom_boxplot(outlier.size=0,coef=0) +
  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #              geom = "point", width = 0.5)+
  scale_colour_manual(values = c("#ff7f0e", "#1f77b4"))+
  #stat_boxplot(geom="errorbar") + 
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),
                stat="identity",position="dodge")+
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.text.x=element_text(size=16,color="black",
                                 angle=45,hjust=1),
        axis.title=element_text(size=16),
        legend.text = element_text(size=24), 
        legend.title = element_text(size=16))+
  xlab("Tissue")+
  ylab("Proportion of Age-Associated Genes")
ggsave("Figure 2D.pdf",width=15)