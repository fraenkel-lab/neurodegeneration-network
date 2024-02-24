library(tximport)
library(DESeq2)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(apeglm)
library(pheatmap)
library(RColorBrewer)

num.ids<-rep(c(seq(5926,5945),seq(5947,5949)),2)
id.num.lanes<-num.ids[order(num.ids)]
lane.ids<-rep(c(1,2),length(id.num.lanes)/2)

## output files from salmon
sample.ids<-paste0("D23-",id.num.lanes,"-",lane.ids,"_quant")
salmon.files<-file.path("salmon_perturbseq_outs/",
                        sample.ids,"quant.sf")

names(salmon.files)<-gsub("_quant","",sample.ids)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
tx2gene<-tx2gene[!is.na(tx2gene$GENEID),]

txi <- tximport(salmon.files, type = "salmon", tx2gene = tx2gene)
names(txi)

bmc.key<-fread("auxiliary_files/BMC_gene_perturb.txt")
sample.groups<-data.table(ID=gsub("-[0-9]_quant","",sample.ids))
bmc.key.groups<-merge(bmc.key,sample.groups,by.x="BMC_ID",by.y="ID")

sampleTable <- data.frame(condition = factor(bmc.key.groups$Gene_perturb),
                          lane = lane.ids,
                          replicate_number=bmc.key.groups$Rep_Number,
                          guide=factor(bmc.key.groups$Guide))
sampleTable$condition_guide<-factor(paste0(bmc.key.groups$Gene_perturb,
                                           "_",bmc.key.groups$Guide))
sampleTable$condition_guide_rep<-factor(paste0(bmc.key.groups$Gene_perturb,
                                               "_",bmc.key.groups$Guide,
                                               bmc.key.groups$Rep_Number))
rownames(sampleTable) <- colnames(txi$counts)

dds <- DESeqDataSetFromTximport(txi, sampleTable, ~condition+guide)

ddsColl <- collapseReplicates(dds,groupby=dds$condition_guide_rep,dds$lane)

dds <- DESeq(ddsColl)

print(resultsNames(dds))
perturbations<-c("CSNK2A1","NOTCH1","HNRNPA2B1")
entrezids<-c("1457","4851","3181")
fig.label<-c("Figure_S7A","Figure_S7B","Figure_5D")
for(i in 1:length(perturbations)){
  res <- results(dds,name=paste0("condition_",perturbations[i],"_vs_Control"))
  
  resLFC <- lfcShrink(dds, coef=paste0("condition_",perturbations[i],
                                       "_vs_Control"),type="apeglm")
  
  name.map<-fread("auxiliary_files/gene_name_mappings.txt")
  lfc.table<-as.data.frame(resLFC)
  lfc.table$entrezgene<-rownames(resLFC)
  
  lfc.table.id<-merge(lfc.table,name.map[,c("entrezgene","hgnc_symbol")],
                      by="entrezgene")
  
  write.table(lfc.table.id,
              file=paste0("shrink_LFC_results_",
                          perturbations[i],"_guide_tech_collapse.txt"),
              sep='\t',quote=F,row.names=F)
  
  
  lfc.labeled<-lfc.table.id[!is.na(lfc.table.id$padj),]
  lfc.sig<-lfc.labeled[lfc.labeled$padj<0.1&abs(lfc.labeled$log2FoldChange)>1,]
  lfc.upreg<-lfc.labeled[lfc.labeled$padj<0.1&lfc.labeled$log2FoldChange>1,]
  lfc.downreg<-lfc.labeled[lfc.labeled$padj<0.1&lfc.labeled$log2FoldChange<(-1),]
  
  ggplot(data=lfc.labeled, aes(y=-log10(padj),
                               x=log2FoldChange
                               #color=label
  )) +
    geom_point(stat="identity",alpha=0.1) +
    geom_point(data=lfc.downreg,aes(y=-log10(padj),x=log2FoldChange),colour="blue")+
    geom_point(data=lfc.upreg,
               aes(y=-log10(padj),x=log2FoldChange),colour="red")+
    #labs(color = "RBP target\n") +
    geom_hline(aes(yintercept=(-log10(0.1))), colour="red", linetype="dashed")+
    geom_vline(aes(xintercept=1), colour="red", linetype="dashed")+
    geom_vline(aes(xintercept=-1), colour="red", linetype="dashed")+
    theme_classic()+
    geom_text_repel(data=lfc.downreg,aes(label=hgnc_symbol),
                    size=6)+
    geom_text_repel(data=lfc.upreg,aes(label=hgnc_symbol),
                    size=6)+
    theme(legend.position = "none")+
    theme(axis.text=element_text(size=16,color="black"),
          axis.title=element_text(size=18,color="black"))+
    xlab(paste0('log2 fold change between ',perturbations[i],' CRISPRi and control'))+
    ylab('negative log10 FDR-adjusted p-value')
  ggsave(paste0(fig.label,".pdf"))
  
}

## data transformation for PCA and plotting
vsd <- vst(ddsColl, blind=FALSE)

cols<-colorRampPalette(c("blue","white","red"))(50)

vsd.sub<-vsd[,colnames(vsd)%in%sampleTable$condition_guide_rep[
  !sampleTable$condition%in%c("HNRNPA2B1")
]]
pcaData <- plotPCA(vsd.sub, intgroup=c("condition","guide",
                                       "replicate_number"), returnData=TRUE,
                   ntop=dim(vsd.sub)[1])
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition,shape=as.factor(guide))) +
  geom_point(size=3) +
  theme_classic()+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=18,color="black"),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))
ggsave("Figure_S6D.pdf")

## plot distributions of knockdowns, guides pooled
tpm.sub<-txi$abundance[rownames(txi$abundance)
                       %in%name.map$entrezgene[
                         name.map$hgnc_symbol%in%c("HNRNPA2B1",
                                                   "NOTCH1",
                                                   "CSNK2A1")
                       ],]
rownames(tpm.sub)<-c("CSNK2A1","HNRNPA2B1","NOTCH1")
figure.subpanel<-c("B","A","C")
for(i in 1:nrow(tpm.sub)){
  temp.dt.gene<-data.table(Value=tpm.sub[i,colnames(tpm.sub)%in%
                                           rownames(df[df$Condition%in%
                                                         rownames(tpm.sub)[i],])],
                           Genotype=rownames(tpm.sub)[i])
  temp.dt.control<-data.table(Value=tpm.sub[i,colnames(tpm.sub)%in%
                                              rownames(df[df$Condition%in%
                                                            c("Control"),])],
                              Genotype="Control")
  plot.df<-rbind(temp.dt.gene,temp.dt.control)
  
  ggplot(plot.df, aes(y=as.numeric(Value), x=as.character(Genotype))) + 
    
    geom_violin()+
    geom_boxplot(width=0.05, alpha=0) +
    geom_point(position=position_jitter(width=0.01)
    )+
    facet_wrap(~Gene,scale="free_y")+
    theme_classic()+
    theme(axis.text=element_text(size=16,color="black"),
          axis.title=element_text(size=24,color="black"))+
    xlab('Genotype')+
    ylab('TPM')
  ggsave(paste0("Figure_S6",figure.subpanel[i],".pdf"))
}

perturbations<-c("NOTCH1","HNRNPA2B1","CSNK2A1")

collected.pathways<-data.table(pathway="")

for(perturbation in perturbations){
  traj.df<-fread(paste0(
    "shrink_LFC_results_",perturbation,
    "_guide_tech_collapse.txt"))
  
  traj.sig<-traj.df[
    !is.na(traj.df$padj)]
  traj.sig<-traj.sig[order(traj.sig$log2FoldChange)]
  
  ranks<-traj.sig$log2FoldChange
  
  ranks[is.infinite(ranks)]<-.Machine$double.xmax
  names(ranks)<-traj.sig$hgnc_symbol
  ranks<-ranks[!duplicated(names(ranks))]
  ranks<-ranks[!is.na(names(ranks))]
  
  pathways<-gmtPathways("auxiliary_files/Reactome_2022.txt")
  names(pathways)<-gsub(" R-HSA-[0-9]+","",names(pathways))
  names(pathways)<-gsub(" \\(GO:[0-9]+\\)","",names(pathways))
  
  set.seed(42)
  fgseaRes <- fgsea(pathways, ranks,nproc=1,
                    minSize=5,maxSize=500)
  print(min(fgseaRes$padj,na.rm=TRUE))
  
  fgsea.sig<-fgseaRes[fgseaRes$padj<0.1,]
  fgsea.sig$group<-perturbation
  
  for(i in 1:nrow(fgsea.sig)){
    fgseaRes$leading_edge[i]<-paste(unlist(fgseaRes$leadingEdge[i]),
                                    collapse=",")
  }
  
  fgseaRes<-fgseaRes[,!which(colnames(fgseaRes)%in%c("leadingEdge")),
                     with=FALSE]
  
  write.table(fgseaRes[order(fgseaRes$padj),],
              paste0(perturbation,"_enriched_pathways_sorted.txt"),
              sep='\t',quote=F,row.names=F)
  
  collected.pathways<-merge(collected.pathways,fgsea.sig,by="pathway",
                            all.x=TRUE,all.y=TRUE)
  
  ## write sig gene tables
  traj.sig<-traj.sig[!is.na(traj.sig$log2FoldChange),]
  
  traj.sig.order<-traj.sig[order(traj.sig$padj,-abs(traj.sig$log2FoldChange)),]
  
  # write.table(traj.sig.order,paste0("RNA_seq_knockdown_results_summary/tables/",
  #                                   perturbation,"_sort_DEGs.txt"),
  #             sep='\t',quote=F,row.names=F)
}

group.1<-collected.pathways[,1:9]
colnames(group.1)<-gsub("\\.x","",colnames(group.1))

group.2<-collected.pathways[,c(1,10:17)]
colnames(group.2)<-gsub("\\.y","",colnames(group.2))

group.3<-collected.pathways[,c(1,18:25)]

pathway.plot.dt<-rbind(group.1)

pathway.plot.dt<-pathway.plot.dt[!is.na(pathway.plot.dt$group)
                                 #&pathway.plot.dt$padj<0.1
]
pathway.plot.order<-pathway.plot.dt[order(-pathway.plot.dt$NES),]

#get top 10 pathways by genotype
pathway.plot.order.top<-pathway.plot.order[,head(.SD,10),by="group"]
pathway.plot.order.top<-pathway.plot.order[pathway.plot.order$padj<0.1]

n.group.enrich<-as.data.frame(table(pathway.plot.order.top$pathway))

pathway.plot.order.top.combine<-pathway.plot.order[pathway.plot.order$pathway%in%
                                                     n.group.enrich$Var1[
                                                       n.group.enrich$Freq==1]]
pathway.plot.order.top.combine$total_pathway_size<-0
pathway.plot.order.top.combine$lead_size<-0

for(i in 1:nrow(pathway.plot.order.top.combine)){
  pathway.plot.order.top.combine$total_pathway_size[i]<-length(pathways[[
    pathway.plot.order.top.combine$pathway[i]]])
  pathway.plot.order.top.combine$lead_size[i]<-length(strsplit(unlist(
    pathway.plot.order.top.combine$leadingEdge[[i]]),","))
}

pathway.plot.order.top.combine$size_proportion<-(
  pathway.plot.order.top.combine$lead_size/
    pathway.plot.order.top.combine$total_pathway_size)

colnames(pathway.plot.order.top.combine)[
  ncol(pathway.plot.order.top.combine)]<-"Proportion of enriched genes in pathway"

# rename long pathway name
pathway.plot.order.top.combine$pathway[pathway.plot.order.top.combine$pathway%in%
       pathway.plot.order.top.combine$pathway[2]
       ]<-"Respiratory Electron Transport, ATP Synthesis"

#limit to the top 10 pathways for HNRNPA2B1
pathway.plot.order.top.combine<-pathway.plot.order.top.combine[
  pathway.plot.order.top.combine$group%in%c("HNRNPA2B1")
]

ggplot(pathway.plot.order.top.combine, aes(y=reorder(pathway,
                                                     NES), x=NES,
)) + 
  geom_point(aes(colour=-log10(padj),size=`Proportion of enriched genes in pathway`))+
  # scale_colour_gradient2(low = scales::muted("blue"),
  #                        mid = "white",
  #                        high = scales::muted("red"),
  # )+
  scale_colour_gradient2(low = scales::muted("red"),
                         #mid = "white",
                         high = scales::muted("blue")
  )+
  #geom_vline(xintercept=1,color="red",linetype="dashed")+
  facet_wrap(~group,scale="free_x",ncol=1)+
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))+
  xlab('Normalized enrichment score')+
  ylab('enriched pathways')
ggsave("Figure_5E.pdf",
       width=20,height=15)

### merge DEGs into one table
ck.degs<-fread("RNA_seq_knockdown_results_summary/tables/CSNK2A1_DEGs.txt")
notch.degs<-fread("RNA_seq_knockdown_results_summary/tables/NOTCH1_DEGs.txt")
hnrn.degs<-fread("RNA_seq_knockdown_results_summary/tables/HNRNPA2B1_DEGs.txt")

ck.notch.merge<-merge(ck.degs,notch.degs,by="entrezgene",all.x=TRUE,
                      all.y=TRUE)
colnames(ck.notch.merge)<-gsub("\\.x","_CSNK2A1",colnames(ck.notch.merge))
colnames(ck.notch.merge)<-gsub("\\.y","_NOTCH1",colnames(ck.notch.merge))

ck.notch.hnrn.merge<-merge(ck.notch.merge,
                           hnrn.degs,by="entrezgene",all.x=TRUE,
                           all.y=TRUE)

colnames(ck.notch.hnrn.merge)[14:19]<-paste0(colnames(
  ck.notch.hnrn.merge)[14:19],"_HNRNPA2B1")

deg.merge.out<-ck.notch.hnrn.merge[,-which(colnames(ck.notch.hnrn.merge)%in%
                                             c("hgnc_symbol_NOTCH1",
                                               "hgnc_symbol_HNRNPA2B1")),with=FALSE]
colnames(deg.merge.out)[7]<-"hgnc_symbol"

write.table(deg.merge.out,
            "Table_S8.txt",
            sep='\t',quote=F,row.names=F)

## get compiled pathway results
ck.enriched.pathways<-fread(
  "CSNK2A1_enriched_pathways.txt")
notch.enriched.pathways<-fread(
  "NOTCH1_enriched_pathways.txt")
hnrn.enriched.pathways<-fread(
  "HNRNPA2B1_enriched_pathways.txt")

ck.not.path.merge<-merge(ck.enriched.pathways,notch.enriched.pathways,
                         by="pathway")
ck.not.hn.path.merge<-merge(ck.not.path.merge,hnrn.enriched.pathways,
                            by="pathway")
colnames(ck.not.hn.path.merge)<-gsub("\\.x","_CSNK2A1",colnames(ck.not.hn.path.merge))
colnames(ck.not.hn.path.merge)<-gsub("\\.y","_NOTCH1",colnames(ck.not.hn.path.merge))

pathway.compiled<-ck.not.hn.path.merge[,-which(colnames(ck.not.hn.path.merge)%in%
                                                 c("size_NOTCH1","size")),with=FALSE]
colnames(pathway.compiled)[15:20]<-paste0(colnames(pathway.compiled)[15:20],
                                          "_HNRNPA2B1")

## fix the leading edge analysis
perturbations<-c("CSNK2A1","NOTCH1","HNRNPA2B1")
deg.list<-list(ck.degs,notch.degs,hnrn.degs)

pathway.compiled.copy<-pathway.compiled
pathway.compiled.copy$leading_edge_size_CSNK2A1<-0
pathway.compiled.copy$leading_edge_size_NOTCH1<-0
pathway.compiled.copy$leading_edge_size_HNRNPA2B1<-0

for(i in 1:length(perturbations)){
  traj.sig<-deg.list[[i]][!is.na(deg.list[[i]]$padj),]
  traj.sig<-traj.sig[order(-traj.sig$log2FoldChange)]
  
  ranks<-traj.sig$log2FoldChange
  ranks[is.infinite(ranks)]<-.Machine$double.xmax
  names(ranks)<-traj.sig$hgnc_symbol
  ranks<-ranks[!duplicated(names(ranks))]
  ranks<-ranks[!is.na(names(ranks))]
  
  for(j in 1:nrow(pathway.compiled.copy)){
    nd.df<-data.table(ID=path,Gene=pathways[[pathway.compiled.copy$pathway[j]]])
    nd.df<-nd.df[nd.df$Gene%in%names(ranks)]
    
    ranks[is.na(ranks)]<-0
    order.ranks<-ranks[order(-ranks)]
    em2 <- GSEA(order.ranks, TERM2GENE = nd.df,seed=42,pvalueCutoff=1,
                minGSSize=1,maxGSSize=1000000)
    
    lead.edge<-unlist(strsplit(em2@result$core_enrichment,"/"))
    lead.edge.size<-length(lead.edge)
    
    pathway.compiled.copy[j,paste0("leading_edge_",
                                   perturbations[i])]<-paste0(
                                     lead.edge,collapse=",")
    pathway.compiled.copy[j,paste0("leading_edge_size_",
                                   perturbations[i])]<-lead.edge.size
  }
}

pathway.compiled.copy.sub<-pathway.compiled.copy

pathway.compiled.copy$leading_edge_size_CSNK2A1<-lengths(
  strsplit(pathway.compiled.copy$leading_edge_CSNK2A1,","))
pathway.compiled.copy$leading_edge_size_NOTCH1<-lengths(
  strsplit(pathway.compiled.copy$leading_edge_NOTCH1,","))
pathway.compiled.copy$leading_edge_size_HNRNPA2B1<-lengths(
  strsplit(pathway.compiled.copy$leading_edge_HNRNPA2B1,","))

write.table(pathway.compiled.copy,
            "Table_S9.txt",
            sep='\t',quote=F,row.names=F)

## plot barplot of cell cycle, DNA damage for NOTCH/CSNK2A1
enriched.pathways<-fread("Table_S9.txt")

ck.sub<-enriched.pathways[enriched.pathways$pathway%in%
                            c("DNA Repair","Cell Cycle, Mitotic",
                              "Cell Cycle Checkpoints",
                              "Activation Of NMDA Receptors And Postsynaptic Events"),
                          c("pathway","NES_CSNK2A1"),with=FALSE]
ck.sub$group<-"CSNK2A1"
colnames(ck.sub)<-gsub("_CSNK2A1","",colnames(ck.sub))

notch.sub<-enriched.pathways[enriched.pathways$pathway%in%
                               c("DNA Repair","Cell Cycle, Mitotic",
                                 "Cell Cycle Checkpoints",
                                 "Activation Of NMDA Receptors And Postsynaptic Events"),
                             c("pathway","NES_NOTCH1"),with=FALSE]
notch.sub$group<-"NOTCH1"
colnames(notch.sub)<-gsub("_NOTCH1","",colnames(notch.sub))

pathway.plot.order.top.combine<-rbind(
  ck.sub,notch.sub
)

pathway.plot.order.top.combine$pathway<-factor(pathway.plot.order.top.combine$pathway,
                                               levels=c("Activation Of NMDA Receptors And Postsynaptic Events",
                                                        "DNA Repair",
                                                        "Cell Cycle, Mitotic",
                                                        "Cell Cycle Checkpoints"))
ggplot(pathway.plot.order.top.combine, aes(y=pathway, x=NES,
)) + 
  geom_col(aes(fill=NES))+
  scale_fill_gradient2(low = scales::muted("blue"),
                       high = scales::muted("red")
  )+
  facet_wrap(~group,nrow=1)+
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        strip.text=element_text(size=16))+
  xlab('Normalized enrichment score')+
  ylab('enriched pathways')
ggsave("Figure_6F.pdf")
