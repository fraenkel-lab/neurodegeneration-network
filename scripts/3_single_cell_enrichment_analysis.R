library("data.table")
library("Seurat")
library("GenomeInfoDb")
library("GenomeInfoDbData")
library("monocle3")
library("harmony")
library("cowplot")
library("SeuratWrappers")
library("ggplot2")
library("GSVA")
library("gprofiler2")
library("reactome.db")
library("fgsea")
library("ordinal")
library("lsa")
library("SeuratDisk")
library("Matrix")

hypergeom_test<-function(x,m,n,k){
  # x: number of successes
  # m: number of successes in population
  # n: size of population
  # k: number of draws
  return(1-phyper(x-1,m,n-m,k))
}

### read in single-nucleus RNA-seq data from synapse (syn2580853)

sc.dat<-Read10X(data.dir="ROSMAP/",gene.column=1)

metadat<-fread("ROSMAP/meta_clin_tech.txt",stringsAsFactors=F)
print(unique(metadat$broad.cell.type))
filt.gene<-fread("ROSMAP/filtered_gene_row_names.txt",stringsAsFactors=F,
                 header=F)

#annotate
rownames(metadat)<-gsub("\\.","-",metadat$TAG)
alz <- CreateSeuratObject(counts = sc.dat, project = "alz",
                          meta.data=metadat)

## batch-correct excitatory neurons
Idents(alz)<-alz$broad.cell.type
alz.ex<-subset(alz,idents=c("Ex"))

## remove doublets inferred from scrublet in 2_single_cell_preprocess_doublet_filter.py
doublets<-fread("doublet_cells.txt",sep='\t',stringsAsFactors=F)

Idents(alz.ex)<-alz.ex$TAG
alz.ex<-alz.ex[,-which(alz.ex$TAG%in%doublets$TAG)]

## run preprocessing functions
alz.ex@assays$Raw_RNA<-alz.ex@assays$RNA
alz.ex<-NormalizeData(alz.ex,verbose = FALSE)
alz.ex<-FindVariableFeatures(alz.ex, selection.method = "mean.var.plot")
alz.ex<-ScaleData(alz.ex, verbose = FALSE)

alz.ex<-RunPCA(alz.ex,pc.genes = alz.ex@var.genes, npcs = 20, verbose = FALSE)

## batch correct with Harmony
alz.ex <- RunHarmony(alz.ex,group.by.vars=c("msex","pmi","Subject",
                                  "sequencingBatch"), plot_convergence = TRUE)

## identify clusters overrepresented by AD patients or healthy controls

## cluster with leiden clustering, from parameter sweep choose resolution=0.1
alz.ex<-FindClusters(alz.ex,algorithm=4,res=0.3)

## determine if clusters are overrepresented by AD patients
alz.ex@meta.data$cell_type<-paste0(alz.ex@meta.data$broad.cell.type,
                                   alz.ex@meta.data$Cluster)

alz.ex@meta.data$AD_status<-"AD"
alz.ex@meta.data$AD_status[
  alz.ex@meta.data$pathology.group=="no-pathology"]<-"Control"

alz.fetch<-FetchData(alz.ex,c("AD_status","cell_type"))
overrep.df<-data.table()
for(clust in unique(alz.ex$cell_type)){
  for(status in as.character(c("AD","Control"))){
    print(c(length(which(alz.fetch$AD_status%in%c(status)&
                           alz.fetch$cell_type%in%c(as.character(clust)))),
            length(which(alz.fetch$AD_status%in%c(status))),
            nrow(alz.fetch),
            length(which(alz.fetch$cell_type%in%c(as.character(clust))))))
    hyper.p<-hypergeom_test(length(which(alz.fetch$AD_status%in%c(status)&
                                alz.fetch$cell_type%in%c(as.character(clust)))),
                            length(which(alz.fetch$AD_status%in%c(status))),
                            nrow(alz.fetch),
                  length(which(alz.fetch$cell_type%in%c(as.character(clust)))))
    
    temp.dt<-data.table(cluster=clust,AD_status=status,pval=hyper.p)
    print(temp.dt)
    
    overrep.df<-rbind(temp.dt,overrep.df)
  }
}

overrep.df$padj<-p.adjust(overrep.df$pval,method="fdr")

## plot enriched clusters and average expression of screen hits
options(repr.plot.width = 15, repr.plot.height = 12)

##calculate average expression
nd.mod<-fread("Table1.txt")
alz.ex@meta.data$`Neurodegeneration Modifiers`<-0

exp.dat<-data.table(FetchData(alz,c(rownames(alz),"cell_type")))
ex.avg<-rowMeans(exp.dat[grep("Ex",exp.dat$cell_type),
                           colnames(exp.dat)%in%nd.mod$`human gene`,
                           with=FALSE])

alz.ex@meta.data$`Neurodegeneration Modifiers`[
  grep("Ex",alz.ex@meta.data$cell_type)]<-ex.avg

alz.ex$AD_enrich<-"Not Enriched"

alz.ex$AD_enrich[alz.ex$cell_type%in%
                   overrep.df$cluster[overrep.df$AD_status=="AD"&
                                        overrep.df$padj<0.1&
                                        grepl("^Ex",overrep.df$cluster)]]<-"Enriched for AD"
alz.ex$AD_enrich[alz.ex$cell_type%in%
                   overrep.df$cluster[overrep.df$AD_status=="Control"&
                                        overrep.df$padj<0.1&
                                        grepl("^Ex",overrep.df$cluster)]]<-"Enriched for Control"

alz.ex$AD_enrich<-as.factor(alz.ex$AD_enrich)

plot1<-DimPlot(alz.ex,group.by=c("AD_enrich"),raster=FALSE,cols=c("#FF7F0E","#1F77B4",
                                                                  "Gray"))& 
  scale_colour_manual(values= c("#FF7F0E","#1F77B4","Gray"),
                      labels=c("Enriched for AD","Enriched for Control","Not Enriched")) & 
  labs(color="Disease Enrichment")+
  theme(axis.text=element_text(size=16,color="black"))+
  theme(axis.title=element_text(size=16)) + 
  theme(legend.text = element_text(size=20)) +
  
  theme(legend.title = element_text(size=20)) &
  xlab("UMAP 1") &
  ylab("UMAP 2") 

plot2<-FeaturePlot(alz.ex,c("Neurodegeneration screen hits")) & 
  scale_colour_gradientn(colors= c("lightgrey", "blue"),
                         name="Log Expression") & 
  theme(axis.text=element_text(size=16,color="black"))+
  theme(axis.title=element_text(size=16)) + 
  theme(legend.text = element_text(size=20)) +
  theme(legend.title = element_text(size=20)) &
  xlab("UMAP 1") &
  ylab("UMAP 2") 

plot1+plot2

ggsave("Figure_S2.pdf",
       width=15,height=12)

## identify differentially abundant genes between AD-associated clusters to 
## all other clusters
ad.high<-unique(overrep.df$cluster[overrep.df$AD_status=="AD"&
                                        overrep.df$padj<0.1])
ad.low<-unique(overrep.df$cluster[overrep.df$AD_status=="Control"&
                                       overrep.df$padj<0.1])

ad.low.unique<-ad.low[-which(ad.low%in%ad.high)]
ad.high.unique<-ad.high[-which(ad.high%in%ad.low)]
print(ad.low.unique)
print(ad.high.unique)

ad.markers <- FindMarkers(alz.ex, ident.1=ad.high,
                             ident.2=ad.low.unique,
                             group.by="cell_type",only.pos = FALSE, min.pct = 0.25,
                             test.use="MAST")

## identify if DEGs are enriched for nd genes

deg.dt<-data.table()
excitatory.neuron.cna.all<-ad.markers
pop.success<-nrow(ad.markers[
  ad.markers$avg_log2FC<1&
    ad.markers$p_val_adj<0.1])

subset.success<-ad.markers[
  ad.markers$avg_log2FC<1&
    ad.markers$p_val_adj<0.1]
subset.success<-subset.success[order(subset.success$p_val_adj,
                                     -abs(subset.success$avg_log2FC))]

nd.success<-nrow(ad.markers[
  ad.markers$avg_log2FC<1&
    ad.markers$p_val_adj<0.1&
    ad.markers$Gene%in%nd.mod$`human gene`])

nd.success<-nrow(subset.success[subset.success$Gene%in%
                                  nd.mod$`human gene`])

n.nd<-length(unique(nd.mod$`human gene`))
print(hypergeom_test(nd.success,pop.success,background,n.nd))

temp.pop<-data.table(`Cell Type`=groups[i],
                     group="All Protein-Coding Genes",
                     DEG_proportion=pop.success/background,
                     stringsAsFactors=F)
pop.ci<-binom.test(pop.success,
                   background,temp.pop$DEG_proportion)$conf.int

temp.pop$upper_CI<-pop.ci[2]
temp.pop$lower_CI<-pop.ci[1]

temp.nd<-data.table(
  `Cell Type`=groups[i],
  group="Neurodegeneration Modifiers",
  DEG_proportion=nd.success/n.nd,
  stringsAsFactors=F)

pop.ci<-binom.test(nd.success,
                   n.nd,temp.nd$DEG_proportion)$conf.int

temp.nd$upper_CI<-pop.ci[2]
temp.nd$lower_CI<-pop.ci[1]

deg.dt<-rbind(temp.pop,temp.nd,deg.dt)


ggplot(deg.dt,
       aes(y=DEG_proportion, x=`Cell Type`,
           fill=group)) + 
  geom_bar(stat="identity",position="dodge")+
  #geom_boxplot(outlier.size=0,coef=0) +
  # stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
  #              geom = "point", width = 0.5)+
  scale_colour_manual(values = c("#1f77b4","#ff7f0e"))+
  #stat_boxplot(geom="errorbar") + 
  geom_errorbar(aes(ymin=lower_CI,ymax=upper_CI),
                stat="identity",position="dodge")+
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16,color="black"),
        legend.text = element_text(size=24),
        legend.title = element_text(size=24))+
  xlab("Cell Type")+
  ylab("Proportion of Differentially Expressed Genes")
ggsave("Figure 2E.pdf",width=15)
