pcna.dat<-fread("auxiliary_files/PCNA_quantification.txt")
pcna.dat<-pcna.dat[1:6,]

plot.dt<-data.table(value=c(pcna.dat$control,pcna.dat$CkIIa_RNAi,
                            pcna.dat$CkIIb_RNAi),
                    Genotype=c(rep("Control",nrow(pcna.dat)),
                               rep("CkIIa RNAi",nrow(pcna.dat)),
                               rep("CkIIb RNAi",nrow(pcna.dat))))
anova.ck2a<-aov(value~factor(Genotype),data=plot.dt)
anova.summary<-summary(anova.ck2a)
tukey<-TukeyHSD(anova.ck2a)

plot.dt$Genotype_factor<-factor(plot.dt$Genotype,
                                levels=c("Control","CkIIa RNAi",
                                         "CkIIb RNAi"))
ggplot(plot.dt, aes(
  x=Genotype_factor, y=value)) + 
  #geom_point(position=position_jitter(w=0.01,h=0)) + 
  #geom_violin()+
  geom_boxplot(outlier.shape = NA) +
  geom_point(size=1,position=position_jitter(width=0.1,
                                             height=0)
  )+
  theme_classic()+
  theme(axis.text.x=element_text(size=16,color="black",
                                 angle=45,hjust=1),
        axis.text.y=element_text(size=16,color="black"),
        axis.title=element_text(size=16,color="black"),
        legend.text=element_text(size=16,color="black"),
        legend.title=element_text(size=16,color="black"))+
  guides(color=guide_legend("Cell type"))+
  xlab('Genotype')+
  ylab('PCNA (#)')
ggsave("Fig6h.pdf")
