library(tximport)
library(DESeq2)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(apeglm)
library(pheatmap)
library(RColorBrewer)

tcpy.rna<-fread("auxiliary_files/genes.htseqcount.cufflinks.TCPY.tsv",sep='\t',stringsAsFactors=F)

## data available on request, contact Clemens Scherzer
tcpy.clin<-fread("clinical_metadata.tsv")

key.iso<-fread("auxiliary_files/isoform_map_hg19.txt",sep='\t',stringsAsFactors=F)

colnames(tcpy.rna)<-gsub("AD_|HC_|_TCPY_[0-9]+_rep[0-9]+|BC",
                         "",colnames(tcpy.rna))
temp.cols<-colnames(tcpy.rna)
txi<-tcpy.rna[,colnames(tcpy.rna)%in%
                tcpy.clin$SOURCE_SUBJECT_ID[
                  tcpy.clin$Sex%in%c("M","F")
                ]&!duplicated(colnames(tcpy.rna)),with=FALSE]
colnames(txi)<-temp.cols[colnames(tcpy.rna)%in%
                           tcpy.clin$SOURCE_SUBJECT_ID[
                             tcpy.clin$Sex%in%c("M","F")
                           ]&
                           !duplicated(colnames(tcpy.rna))]

txi<-txi[,!duplicated(colnames(txi)),with=FALSE]
rownames(txi)<-tcpy.rna$tracking_id

tcpy.clin<-tcpy.clin[tcpy.clin$SOURCE_SUBJECT_ID%in%colnames(txi),]
rownames(tcpy.clin)<-tcpy.clin$SOURCE_SUBJECT_ID
tcpy.clin$NPDX<-factor(tcpy.clin$NPDX)

tcpy.clin.nobc<-tcpy.clin[!tcpy.clin$SOURCE_SUBJECT_ID%in%tcpy.bc$clin.id]

dds <- DESeqDataSetFromMatrix(txi, tcpy.clin, ~NPDX)
#keep <- rowSums(counts(dds)) >= 10
#dds <- dds[keep,]

dds <- DESeq(dds)

print(resultsNames(dds))
perturbations<-c("HC")
for(i in 1:length(perturbations)){
  res <- results(dds,name=paste0("NPDX_",perturbations[i],"_vs_AD"))
  
  resLFC <- lfcShrink(dds, coef=paste0("NPDX_",perturbations[i],
                                       "_vs_AD"),type="normal")
  
  lfc.table<-as.data.frame(resLFC)
  lfc.table$`Gene stable ID version`<-rownames(txi)
  
  lfc.table.id<-merge(lfc.table,key.iso[,c("Gene stable ID version",
                                           "HGNC symbol")],
                      by="Gene stable ID version")
  
  write.table(lfc.table.id,
              file=paste0("TCPY_counts/TCPY_AD_DEG_normshrink.txt"),
              sep='\t',quote=F,row.names=F)
  
  
  lfc.labeled<-lfc.table.id[!is.na(lfc.table.id$padj)&
                              !is.na(lfc.table.id$log2FoldChange),]
  lfc.labeled<-lfc.labeled[!duplicated(lfc.labeled$`Gene stable ID version`),]
  lfc.labeled$log2FoldChange<-lfc.labeled$log2FoldChange*-1
  lfc.sig<-lfc.labeled[lfc.labeled$padj<0.1&abs(lfc.labeled$log2FoldChange)>1,]
  lfc.upreg<-lfc.labeled[lfc.labeled$padj<0.1&lfc.labeled$log2FoldChange>1,]
  lfc.downreg<-lfc.labeled[lfc.labeled$padj<0.1&lfc.labeled$log2FoldChange<(-1),]
  
  ggplot(data=lfc.labeled, aes(y=-log10(padj), 
                               x=log2FoldChange
                               #color=label
  )) + 
    geom_point(stat="identity",alpha=1) +
    geom_point(data=lfc.downreg,aes(y=-log10(padj),x=log2FoldChange),colour="blue")+
    geom_point(data=lfc.upreg,
               aes(y=-log10(padj),x=log2FoldChange),colour="red")+
    #labs(color = "RBP target\n") +
    geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
    geom_vline(aes(xintercept=1), colour="red", linetype="dashed")+
    geom_vline(aes(xintercept=-1), colour="red", linetype="dashed")+
    theme_classic()+
    geom_text_repel(data=lfc.downreg,aes(label=`HGNC symbol`),
                    size=6)+
    geom_text_repel(data=lfc.upreg,aes(label=`HGNC symbol`),
                    size=6)+
    theme(legend.position = "none")+
    theme(axis.text=element_text(size=16,color="black"),
          axis.title=element_text(size=18,color="black"))+
    xlab(paste0('log2 fold change between ',perturbations[i],' AD and control'))+
    ylab('negative log10 FDR-adjusted p-value')
  ggsave(paste0("plots/tcpy_degs_volcano_norm.pdf"))
  ggsave(paste0("plots/tcpy_degs_volcano_norm.png"))
  
}

## plot cell type marker genes according to Allen
rownames(dds)<-tcpy.rna$tracking_id
marker.genes<-c("GAD1","ADARB2","LAMP5","PAX6","VIP",
                "LHX6","SST","PVALB","SLC17A7","CUX2",
                "RORB","THEMIS","FEZF2","CTGF","AQP4",
                "PDGFRA","MOG","CDLN5","MUSTN1","FYB")
marker.keys<-unique(key.iso$`Gene stable ID version`[key.iso$`HGNC symbol`%in%
                                                       marker.genes&
                                                       key.iso$`Gene stable ID version`%in%
                                                       tcpy.rna$tracking_id])

names(marker.keys) <- marker.keys

pc.function<-function(x){
  print(x)
  y<-plotCounts(dds,x,"NPDX",returnData=TRUE)
  y$feature <- x
  return(y)
}
df <- lapply(marker.keys,pc.function)

df <- do.call(rbind, df)
df$gene_pat<-rownames(df)

df.gene.label<-merge(key.iso[,c("Gene stable ID version",
                                "Gene name")],
                     df,by.x="Gene stable ID version",
                     by.y="feature")

df.gene.label.nodup<-df.gene.label[
  !duplicated(df.gene.label$gene_pat)]

df.gene.label.nodup$cell_type<-"unknown"
df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("ADARB2",
                                       "VIP","GAD1",
                                       "LAMP5","PAX6","LHX6",
                                       "SST","PVALB"
  )]<-"Inhibitory neuron"
df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("RORB",
                                       "SLC17A7","CUX2",
                                       "THEMIS","FEZF2",
                                       "CTGF")]<-"Excitatory neuron"

df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("AQP4")]<-"Astrocyte"

df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("PDGFRA")]<-"Oligodendrocyte progenitor cell"

df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("MOG")]<-"Oligodendrocyte"

df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("CDLN5")]<-"Endothelial cell"

df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("FYB")]<-"Microglial cell"

df.gene.label.nodup$cell_type[
  df.gene.label.nodup$`Gene name`%in%c("MUSTN1")]<-"Pericyte"

df.gene.label.nodup$`Gene name`<-factor(
  df.gene.label.nodup$`Gene name`,
  levels=c("ADARB2",
           "VIP","GAD1",
           "LAMP5","PAX6",
           "LHX6",
           "SST","PVALB","RORB",
           "SLC17A7","CUX2",
           "THEMIS","FEZF2",
           "CTGF","AQP4","PDGFRA","MOG","FYB","MUSTN1"))

ggplot(df.gene.label.nodup, aes(
  x=`Gene name`, y=count,color=cell_type)) + 
  #geom_point(position=position_jitter(w=0.01,h=0)) + 
  #geom_violin(width=2)+
  geom_boxplot(width=1,outlier.shape = NA) +
  geom_point(size=1,position=position_jitter(width=0.1)
  )+
  #facet_wrap(~Gene,scale="free_y")+
  scale_color_manual(values=c("#aec7e8","#ffbb78","#2ca02c","#98df8a",
                              "#9467bd","#c5b0d5","#8c564b"))+
  theme_classic()+
  theme(axis.text.x=element_text(size=16,color="black",
                                 angle=45,hjust=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.title=element_text(size=16,color="black"),
        legend.text=element_text(size=16,color="black"),
        legend.title=element_text(size=16,color="black"))+
  guides(color=guide_legend("Cell type"))+
  xlab('Gene')+
  ylab('Count')
ggsave(paste0("TCPY_cell_type_marker_dist.png"),
       width=15)
ggsave(paste0("TCPY_cell_type_marker_dist.pdf"),
       width=15)