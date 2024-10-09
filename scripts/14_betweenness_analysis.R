library(data.table)
library(fgsea)
library(reactome.db)
library(igraph)

## FOR AD create network of nd modifiers and nodes associated with DNA damage
## assess betweenness and see which has no previous association with AD



pathways<-gmtPathways("auxiliary_files/Reactome_2022.txt")
names(pathways)<-gsub(" R-HSA-[0-9]+","",names(pathways))

damage.paths<-c("Base Excision Repair","Base-Excision Repair, AP Site Formation",
                "DNA Damage Bypass","DNA Damage Reversal",
                "DNA Damage Recognition In GG-NER","DNA Damage/Telomere Stress Induced Senescence",
                "DNA Double Strand Break Response","DNA Double-Strand Break Repair",
                "DNA Repair","Dual Incision In GG-NER","Dual Incision In TC-NER",
                "G1/S DNA Damage Checkpoints","G2/M DNA Damage Checkpoint",
                "Gap-filling DNA Repair Synthesis And Ligation In GG-NER",
                "Gap-filling DNA Repair Synthesis And Ligation In TC-NER",
                "HDR Thru Homologous Recombination (HRR)",
                "HDR Thru Homologous Recombination (HRR) Or Single Strand Annealing (SSA)",
                "HDR Thru MMEJ (alt-NHEJ)","HDR Thru Single Strand Annealing (SSA)",
                "Homology Directed Repair","Mismatch Repair",
                "Mismatch Repair (MMR) Directed By MSH2:MSH3 (MutSbeta)",
                "Mismatch Repair (MMR) Directed By MSH2:MSH6 (MutSalpha)",
                "Nonhomologous End-Joining (NHEJ)","Nucleotide Excision Repair",
                "Processing Of DNA Double-Strand Break Ends",
                "Recognition Of DNA Damage By PCNA-containing Replication Complex",
                "Recruitment And ATM-mediated Phosphorylation Of Repair And Signal Proteins At DNA Double Strand Breaks",
                "Sensing Of DNA Double Strand Breaks","TP53 Regulates Transcription Of DNA Repair Genes",
                "Transcription-Coupled Nucleotide Excision Repair (TC-NER)","Translesion Synthesis By POLH",
                "Translesion Synthesis By POLI","Translesion Synthesis By POLK",
                "Translesion Synthesis By REV1","Translesion Synthesis By Y Family DNA Polymerases Bypasses Lesions On DNA Template",
                "p53-Dependent G1 DNA Damage Response","Fanconi Anemia Pathway"
)
## get genes in DNA damage pathways
dt<-data.table()
for(path in damage.paths){
  temp.dt<-data.table(Gene=unlist(pathways[[path]]))
  dt<-rbind(dt,temp.dt)
}

nd.mod<-fread("Table1.txt")
edges<-fread("auxiliary_files/full_network_edges_abeta_tau_lipid.txt")
# edge.damage.nd.sub<-edges[edges$source%in%dt$Gene|
#                             edges$target%in%dt$Gene]
edge.damage.nd.sub<-edges[(edges$source%in%dt$Gene&
                             edges$target%in%dt$Gene)|
                            (edges$target%in%dt$Gene&
                               edges$source%in%nd.mod$`human gene`)|
                            (edges$source%in%dt$Gene&
                               edges$target%in%nd.mod$`human gene`)|
                            (edges$target%in%nd.mod$`human gene`&
                               edges$source%in%nd.mod$`human gene`)]

ad.assoc<-fread("openTargets/AD_association.tsv")

damage.genes.nosoc<-nodes[(nodes$NAME%in%ad.assoc$symbol[
  (ad.assoc$geneticAssociations%in%c("No data")&
     ad.assoc$somaticMutations%in%c("No data")&
     ad.assoc$drugs%in%c("No data")&
     ad.assoc$pathwaysSystemsBiology%in%c("No data")&
     ad.assoc$rnaExpression%in%c("No data")&
     ad.assoc$animalModels%in%c("No data"))])&
    nodes$NAME%in%nd.mod$`human gene`
]

dna.damage.genes.nosoc<-nodes[(nodes$NAME%in%ad.assoc$symbol[
  (ad.assoc$geneticAssociations%in%c("No data")&
     ad.assoc$somaticMutations%in%c("No data")&
     ad.assoc$drugs%in%c("No data")&
     ad.assoc$pathwaysSystemsBiology%in%c("No data")&
     ad.assoc$rnaExpression%in%c("No data")&
     ad.assoc$animalModels%in%c("No data"))])&
    nodes$NAME%in%dt$Gene]

dna.damage.genes.soc<-nodes[!(nodes$NAME%in%ad.assoc$symbol[
  (ad.assoc$geneticAssociations%in%c("No data")&
     ad.assoc$somaticMutations%in%c("No data")&
     ad.assoc$drugs%in%c("No data")&
     ad.assoc$pathwaysSystemsBiology%in%c("No data")&
     ad.assoc$rnaExpression%in%c("No data")&
     ad.assoc$animalModels%in%c("No data"))])&
    nodes$NAME%in%dt$Gene]

damage.graph<-graph_from_edgelist(as.matrix(
  edge.damage.nd.sub),directed = FALSE)
dna.damage.between<-betweenness(damage.graph)  
dna.damage.between.order<-dna.damage.between[order(-dna.damage.between)]

dna.damage.no.lit<-dna.damage.between.order[names(
  dna.damage.between.order)%in%damage.genes.nosoc$NAME]

dna.damage.lit.df<-data.table(Gene=names(dna.damage.between.order),
                              betweenness=unlist(dna.damage.between.order))

write.table(dna.damage.lit.df,"betweenness_damage_analysis.txt",
            sep='\t',quote=F)