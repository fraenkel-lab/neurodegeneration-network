library(data.table)
library(fgsea)
library(igraph)
library(ggplot2)
library(see)
library(readxl)

nodes<-fread("auxiiliary_files/cytoscape_node_anno.txt")
edges<-fread("auxiliary_files/full_network_edges_abeta_tau_lipid.txt")

perturb.degs<-read_xlsx(
  "Supplementary\ Table\ 8.xlsx")

perturbations<-c("NOTCH1","HNRNPA2B1","CSNK2A1")

## generate ego graphs, then plot fold change relative
## to distance from perturbed gene

ad.network<-graph_from_edgelist(as.matrix(
  edges),directed = FALSE)

collected.pathways<-data.table(pathway="")

for(perturbation in perturbations){
  pathways<-list()
  
  dist.df<-stack(as.data.frame(distances(ad.network,
                                         perturbation, 
                                         unlist(ego(ad.network,10, 
                                                    perturbation)))))
  dist.df.merge<-merge(dist.df[dist.df$values>0,],
                       perturb.degs[,c("hgnc_symbol",
                                       paste0("log2FoldChange_",perturbation)),
                                    with=FALSE],
                       by.x="ind",by.y="hgnc_symbol")
  
  dist.df.merge$values<-ordered(dist.df.merge$values,
                                levels = unique(
                                  dist.df.merge$values)[order(unique(
                                    dist.df.merge$values),
                                    decreasing=TRUE)])
  dist.df.merge$absolute_FC<-abs(dist.df.merge[[paste0("log2FoldChange_",
                                                       perturbation)]])
  print(kruskal.test(as.formula(paste0(
    "absolute_FC ~ values")), 
    data = dist.df.merge))
  
  dist.df.merge$degree<-ordered(dist.df.merge$values,
                                levels = unique(
                                  dist.df.merge$values)[order(unique(
                                    dist.df.merge$values),
                                    decreasing=TRUE)])
  #dist.df.merge$degree<-as.numeric(dist.df.merge$degree)
  
  ggplot(dist.df.merge,aes(x=degree,
                           y=abs(.data[[paste0("log2FoldChange_",perturbation)]])))+
    geom_violindot(#dots_size=0.5,
      #position_dots = ggplot2::position_nudge(x = -0.0005, y = 0),
      binwidth=0.01)+
    #geom_boxplot(width=0.1,outlier.shape=NA)+
    #geom_jitter(width=0.1)+
    #geom_smooth(method="lm")+
    theme_classic()+
    theme(axis.text=element_text(size=16,color="black"),
          axis.title=element_text(size=16),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16))+
    ylab('absolute log2FoldChange KO vs control')+
    xlab(paste0('Degrees from ',perturbation))
  ggsave(paste0(perturbation,"_lfc_degree_plot.pdf"),
         width=15)
  
  for(i in 1:length(unique(nodes$louvain_clusters))){
    pathways[[i]]<-perturb.degs$hgnc_symbol[perturb.degs$hgnc_symbol%in%
                                              nodes$NAME[nodes$louvain_clusters%in%c(i)]]
    
    names(pathways)[i]<-paste0("louvain_cluster ",
                               unique(nodes$louvain_clusters)[i])
  }
  
  for(degree in unique(dist.df.merge$values)){
    pathways[[length(pathways)+1]]<-dist.df.merge$ind[
      dist.df.merge$values==degree]
    names(pathways)[length(pathways)]<-paste0("nodes ",degree,
                                              " away from ",perturbation)
    pathways[[length(pathways)+1]]<-dist.df.merge$ind[
      dist.df.merge$values<=degree]
    names(pathways)[length(pathways)]<-paste0("nodes at most ",degree,
                                              " away from ",perturbation)
  }
  
  perturb.order<-perturb.degs[order(
    perturb.degs[[paste0("log2FoldChange_",perturbation)]],
    decreasing=TRUE)]
  ranks<-perturb.degs[[paste0("log2FoldChange_",perturbation)]]
  
  ranks[is.infinite(ranks)]<-.Machine$double.xmax
  names(ranks)<-perturb.order$hgnc_symbol
  ranks<-ranks[!duplicated(names(ranks))]
  ranks<-ranks[!is.na(names(ranks))]
  
  # pathways<-gmtPathways("enrichr_libs/Reactome_2022.txt")
  # names(pathways)<-gsub(" R-HSA-[0-9]+","",names(pathways))
  # names(pathways)<-gsub(" \\(GO:[0-9]+\\)","",names(pathways))
  
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
  
  collected.pathways<-merge(collected.pathways,fgsea.sig,by="pathway",
                            all.x=TRUE,all.y=TRUE)
  
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
#pathway.plot.order<-pathway.plot.dt[order(pathway.plot.dt$padj),]

#pathway.plot.order.top<-pathway.plot.order[,head(.SD,10),by="group"]
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

# pathway.plot.order.top.combine$pathway[pathway.plot.order.top.combine$pathway%in%
#       pathway.plot.order.top.combine$pathway[2]]<-"Respiratory Electron Transport, ATP Synthesis"
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
ggsave("top_differential_pathways_stack.pdf",
       width=20,height=15)
ggsave("top_differential_pathways_stack.png",
       width=20,height=15)

#barplot
ggplot(pathway.plot.order.top.combine, aes(y=reorder(pathway,
                                                     NES), x=NES,
)) + 
  geom_col(aes(fill=NES))+
  # scale_colour_gradient2(low = scales::muted("blue"),
  #                        mid = "white",
  #                        high = scales::muted("red")
  # )+
  scale_fill_gradient2(high = scales::muted("red"),
                       mid = "white",
                       low = scales::muted("blue"))+
  #geom_vline(xintercept=1,color="red",linetype="dashed")+
  facet_wrap(~group)+
  theme_classic()+
  #coord_flip()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16))+
  xlab('Normalized enrichment score')+
  ylab('enriched pathways')
ggsave("shared_enriched_pathways_hnrnpa2b1_notch1.pdf",
       width=15)
ggsave("shared_enriched_pathways_hnrnpa2b1_notch1.png",
       width=15)

## instead, identify number of DEGs per cluster and 
## determine if there's hypergeometric enrichment

stat.df<-data.table()

for(perturbation in perturbations){
  deg<-perturb.degs[perturb.degs[[paste0("padj_",
                                         perturbation)]]<0.1]
  for(cluster in unique(nodes$louvain_clusters)){
    x<-length(which(nodes$louvain_clusters%in%
                      c(cluster)&
                      nodes$NAME%in%deg$hgnc_symbol))
    m<-nrow(deg)
    n<-nrow(perturb.degs)
    k<-length(which(nodes$louvain_clusters%in%c(cluster)))
    
    hyper_p<-hypergeom_test(x,m,n,k)
    
    temp.df<-data.table(perturb=perturbation,
                        louvain_cluster=cluster,
                        p_val=hyper_p)
    stat.df<-rbind(stat.df,temp.df)
  }
}

stat.df$padj<-p.adjust(stat.df$p_val,method="fdr")
