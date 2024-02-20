## use Gprofiler to identify enriched pathways in the neurodegeneration
## screen hits

library(data.table)
library(gprofiler2)
library(ggplot2)

nd.mod<-fread("Table1.txt")

#use gProfiler to identify the enriched pathways among the hits from
#the neurodegeneration screen
gostres <- gost(query = nd.mod$`Drosophila gene`,
                organism = "dmelanogaster", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                measure_underrepresentation = FALSE, evcodes = TRUE,
                user_threshold = 0.05,correction="g_SCS",
                domain_scope = "annotated", 
                numeric_ns = "", as_short_link = FALSE)

enrich.k.r<-gostres$result

enrich.k.r<-enrich.k.r[order(enrich.k.r$p_value),]
enrich.k.r<-enrich.k.r[enrich.k.r$term_size<500&
                         enrich.k.r$intersection_size>3,]

enrich.k.r$overlapping_genes<-unlist(enrich.k.r$intersection)

write.table(enrich.k.r[,c("term_name","p_value",
                          "term_id","overlapping_genes")],
            "Table_S1.csv",
            sep=',',quote=F,row.names=F)

## plot top pathways for inspection
enrich.paths<-fread("Table_S1.csv",
                    sep=',',fill=TRUE)

enrich.paths$p_value<-as.numeric(enrich.paths$p_value)
ggplot(enrich.paths[1:20,], aes(y=reorder(term_name,
                                          -log10(p_value)), x=-log10(p_value))) + 
  geom_bar(stat="identity",position="dodge")+
  theme_classic()+
  theme(axis.text=element_text(size=16,color="black"),
        axis.title=element_text(size=16))+
  xlab(bquote("negative " ~ log[10] ~ " FDR-adjusted p-value"))+
  ylab('Enriched pathways')
