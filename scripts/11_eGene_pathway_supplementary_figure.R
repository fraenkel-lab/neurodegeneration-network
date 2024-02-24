library("data.table")
library("GSVA")
library("GSVAdata")
library("stringr")
library("edgeR")
data("c2BroadSets")
data("genderGenesEntrez")

tcpy.rna<-fread("auxiliary_files/genes.htseqcount.cufflinks.TCPY.tsv",
                sep='\t',stringsAsFactors=F)

#get the gene names from the tracking IDs
key<-fread("auxiliary_files/isoform_map_hg19.txt",sep='\t',stringsAsFactors=F)

tcpy.merge<-merge(tcpy.rna,key,by.x="tracking_id",by.y="Gene stable ID version")

tcpy.merge.unique<-tcpy.merge[!duplicated(tcpy.merge$`Gene name`)]

#get numerical expression matrix
obs_df_t<-data.matrix(tcpy.merge.unique[,grep("AD|HC",
                                             colnames(tcpy.merge.unique)),
                                       with=FALSE])

rownames(obs_df_t)<-tcpy.merge.unique$`Gene name`
gene_hgnc<-tcpy.merge.unique[,c("Gene name"),with=FALSE]

gene_mappings<-fread('auxiliary_files/gene_name_mappings.txt',stringsAsFactors=F)

## merge to get the entrezgene ids
gene_hgnc_merge<-merge(gene_hgnc,gene_mappings,by.x='Gene name',by.y='hgnc_symbol')

gene_hgnc_merge<-gene_hgnc_merge[!is.na(gene_hgnc_merge$`Gene name`),]
gene_hgnc_merge<-gene_hgnc_merge[!duplicated(gene_hgnc_merge$`Gene name`)]
gene_hgnc_merge<-gene_hgnc_merge[!duplicated(gene_hgnc_merge$entrezgene)]

obs_df_t<-obs_df_t[which(rownames(obs_df_t)%in%gene_hgnc_merge$`Gene name`),]
gene_hgnc_merge<-gene_hgnc_merge[which(gene_hgnc_merge$`Gene name`%in%
                                         rownames(obs_df_t)),]

#use EdgeR's CPM transform
exp_mat<-cpm(obs_df_t)
exp_dt<-as.matrix(exp_mat)

colnames(exp_dt)<-colnames(obs_df_t)
rownames(exp_dt)<-gene_hgnc_merge$entrezgene

canonicalC2BroadSets <- c2BroadSets[c(grep("^REACTOME", names(c2BroadSets)))]                          

MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"), setName="MSY")

XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"), setName="XiE")

canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))

esrnaseq <- gsva(exp_dt, canonicalC2BroadSets, method="gsva", min.sz=5, 
                 max.sz=500,parallel.sz=1,
                 kcdf="Poisson")

esrnalim<-data.matrix(esrnaseq)
esrnalim.table<-data.table(esrnalim)
esrnalim.table$pathway<-rownames(esrnalim)
esrnalim.table$pathway<-gsub("^REACTOME_","",esrnalim.table$pathway)

exp_mat<-cpm(obs_df_t)
exp_dt<-as.matrix(exp_mat)

colnames(exp_dt)<-colnames(obs_df_t)
rownames(exp_dt)<-gene_hgnc_merge$entrezgene

eqt<-fread("Table_2.txt")
colnames(eqt)[3]<-"Gene"

exp_dt<-exp_dt[rownames(obs_df_t)%in%eqt$Gene,]
obs_df_t<-obs_df_t[rownames(obs_df_t)%in%eqt$Gene,]

path<-unique(esrnalim.table$pathway)

cor.table<-data.table()

#correlate GSVA gene signatures with eGene expression
for(pt in path){
  print(pt)
  
  for(i in 1:nrow(obs_df_t)){
    cor.out<-cor.test(unlist(esrnalim.table[esrnalim.table$pathway%in%
                                              c(pt),grep("AD|HC",
                                        colnames(esrnalim.table)),with=FALSE]),
                      exp_mat[i,grep("AD|HC",colnames(exp_mat))],
                      method="pearson")
    
    temp.dt<-data.table(gene=rownames(obs_df_t)[i],
                        pathway=pt,
                        coeff=cor.out$estimate,
                        pval=cor.out$p.val)
    
    cor.table<-rbind(cor.table,temp.dt)
  }
}

cor.table$padj<-p.adjust(cor.table$pval,method="fdr")

print(cor.table[cor.table$gene%in%eqt$Gene])
print(length(which(cor.table$padj<0.1)))

pathway.cors<-cor.table[!cor.table$pathway%in%
                             c("MSY","XiE","receptor_ligand_interactions")]
pathway.cors$pathway<-gsub("_"," ",pathway.cors$pathway)

pathway.cors$padj<-p.adjust(pathway.cors$pval,method="fdr")

pathway.cors.sig<-pathway.cors[pathway.cors$padj<0.1,]

pathway.order<-pathway.cors.sig[order(pathway.cors.sig$padj)]

#get top 5 correlated pathways by eGene
top.egene.paths<-pathway.order[,head(.SD,5),by="gene"]

mat.dt<-data.table(pathway=unique(top.egene.paths$pathway))

for(i in 1:length(unique(top.egene.paths$gene))){
  temp.gene<-unique(top.egene.paths$gene)[i]
  temp.dt<-data.table(pathway=pathway.cors.sig$pathway[
    pathway.cors.sig$gene%in%temp.gene],
    value=pathway.cors.sig$coeff[
      pathway.cors.sig$gene%in%
        temp.gene])
  colnames(temp.dt)[2]<-temp.gene
  
  mat.dt<-merge(mat.dt,temp.dt,by="pathway",
                all.x=TRUE)
}

#get all numerical values for heatmap
mat.heat<-data.matrix(mat.dt[,2:ncol(mat.dt)])
mat.heat[is.na(mat.heat)]<-0
mat.heat<-t(mat.heat)
rownames(mat.heat)<-unique(top.egene.paths$gene)
colnames(mat.heat)<-mat.dt$pathway

breaks<-c(seq(min(mat.heat,
                  na.rm=TRUE),0,length.out=ceiling(51/2)+1),
          seq(max(mat.heat,
                  na.rm=TRUE)/51,max(
                    mat.heat,na.rm=TRUE),
              length.out=floor(51/2)))

cols<-colorRampPalette(c("blue","white","red"))(52)

f1<-circlize::colorRamp2(breaks,cols)

breaks2<-c(seq(min(eqt.valid$BETA,na.rm=TRUE),0,length.out=ceiling(51/2)+1),
           seq(max(eqt.valid$BETA,
                   na.rm=TRUE)/51,max(
                     eqt.valid$BETA,na.rm=TRUE),
               length.out=floor(51/2)))
cols2<-colorRampPalette(c("blue","white","red"))(52)
f2<-circlize::colorRamp2(breaks2,cols2)

eqt.beta<-eqt[eqt$Egene3%in%rownames(mat.heat)]
eqt.beta<-eqt.beta[match(rownames(mat.heat),eqt.beta$Egene3)]
row_ha = rowAnnotation(`Beta Coefficient` = eqt.beta$BETA,
                       col=list(`Beta Coefficient`=f2))

pdf("Table_S3.pdf",
    width=20,height=8)

Heatmap(mat.heat,
        name="Pearson Correlation",
        show_row_names = TRUE,
        show_column_names = TRUE,
        col=f1,
        column_names_rot = 45,
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        right_annotation = row_ha,
        row_names_gp = gpar(fontsize = 16),
        column_names_gp = gpar(fontsize = 16))

dev.off()