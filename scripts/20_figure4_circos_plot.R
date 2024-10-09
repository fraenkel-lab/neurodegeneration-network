library(data.table)
library(circlize)
library(ggplot2)

data.net.count<-fread("auxiliary_files/n_node_type_by_louvain_cluster.txt")
data.net.count$general_datatype[data.net.count$general_datatype%in%
                                  c("drosophila phoshoproteomics")]<-"drosophila phosphoproteomics"

pathways.dt<-data.table(louvain_clusters=seq(0,max(data.net.count$louvain_clusters)),
                        pathway=c("Condensation of prometaphase chromosomes",
                                  "Golgi complex activity","Double-stranded DNA break repair",
                                  "Gɑq signaling events","NOTCH signaling","Hedgehog signaling",
                                  "Postsynaptic activity","Postsynaptic activity",
                                  "Autophagy","Insulin signaling",
                                  "Mitotic regulation","Electron transport chain",
                                  "Mitochondrial biogenesis",
                                  "Clathrin-mediated endocytosis"))

data.net.count.anno<-merge(data.net.count,pathways.dt,by="louvain_clusters")

color.map<-fread("auxiliary_files/pie_chart_color_map.txt",header=FALSE)

edges<-fread("auxiliary_files/full_network_edges_abeta_tau_lipid.txt")
nodes<-fread("auxiliary_files/cytoscape_anno_nodes.txt")
nodes$general_datatype[nodes$general_datatype%in%
                         c("drosophila phoshoproteomics")]<-"drosophila phosphoproteomics"
n.edges.between.clusts<-data.table()

## chord df simple
chord.df<-data.table()
nodes$louvain_clusters[nodes$louvain_clusters==7]<-6
for(clust1 in unique(nodes$louvain_clusters)){
  for(clust2 in unique(nodes$louvain_clusters)){
    clust1.df<-nodes[nodes$louvain_clusters%in%c(clust1)]
    clust2.df<-nodes[nodes$louvain_clusters%in%c(clust2)]
    n_edges<-nrow(edges[(edges$source%in%clust1.df$NAME&
                           edges$target%in%clust2.df$NAME)|
                          (edges$target%in%clust1.df$NAME&
                             edges$source%in%clust2.df$NAME)])
    
    temp.df<-data.table(cluster1=clust1,
                        cluster2=clust2,
                        clust1_name=pathways.dt$pathway[
                          pathways.dt$louvain_clusters%in%c(clust1)],
                        clust2_name=pathways.dt$pathway[
                          pathways.dt$louvain_clusters%in%c(clust2)],
                        n_nodes_in_edge=n_edges)
    
    
    chord.df<-rbind(chord.df,temp.df)
  }
}

chord.df$min_clust<-do.call(pmin, chord.df[,c("cluster1","cluster2")])
chord.df$max_clust<-do.call(pmax, chord.df[,c("cluster1","cluster2")])

chord.df.nodup<-chord.df[!duplicated(chord.df[,c("min_clust",
                                                 "max_clust")]),
                         c("clust1_name","clust2_name","n_nodes_in_edge")]

all_species = unique(c(chord.df.nodup$clust1_name,
                       chord.df.nodup$clust2_name))

## get tab20 pallette, skirting the reds
color.pal<-c("#aec7e8","#ffbb78","#2ca02c","#98df8a",
             "#9467bd","#c5b0d5","#8c564b","#c49c94",
             "#e377c2","#f7b6d2","#7f7f7f","#c7c7c7",
             "#bcbd22")

color_species = structure(color.pal, 
                          names = all_species)

colnames(chord.df.nodup)<-c("from","to","value")
chordDiagram(chord.df.nodup, #annotationTrack = c("grid", "axis"),
             grid.col = color_species, 
             #order=unique(chord.df.nodup$clust1_name),
             directional = FALSE#,
             # preAllocateTracks = list(
             #   track.height = 0.04,
             #   track.margin = c(0.02, 0)
             # )
)

pdf("Fig4a.pdf")
chordDiagram(chord.df.nodup, grid.col = color_species, 
             annotationTrack = "grid", 
             preAllocateTracks = list(track.height = max(
               strwidth(unlist(dimnames(chord.df.nodup))))))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA) # here set bg.border to NA is important
dev.off()
