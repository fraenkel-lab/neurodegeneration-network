## helper function for identifying overlaps, from IDR2D package
determine_anchor_overlap <- function(rep1_anchor, rep2_anchor, max_gap = -1L) {
  rep1_ranges <- GenomicRanges::makeGRangesFromDataFrame(rep1_anchor,
                                                         keep.extra.columns = TRUE)
  rep2_ranges <- GenomicRanges::makeGRangesFromDataFrame(rep2_anchor,
                                                         keep.extra.columns = TRUE)
  
  # adjust seq levels
  seq_levels_r1 <- GenomeInfoDb::seqlevels(rep1_ranges)
  seq_levels_r2 <- GenomeInfoDb::seqlevels(rep2_ranges)
  combined_seq_levels <- stringr::str_sort(union(seq_levels_r1,
                                                 seq_levels_r2))
  GenomeInfoDb::seqlevels(rep1_ranges) <- combined_seq_levels
  GenomeInfoDb::seqlevels(rep2_ranges) <- combined_seq_levels
  
  # get overlap between replicates, accept 1000 bp gap
  overlap_df <- data.frame(GenomicRanges::findOverlaps(rep1_ranges,
                                                       rep2_ranges,
                                                       maxgap = max_gap))
  colnames(overlap_df) <- c("rep1_idx", "rep2_idx")
  return(overlap_df)
}

eqt.valid<-fread("Table2.txt")

# reference regions for the SNPs that are in LD with the eQTL of interest
# in the TCPY eQTL analysis
phys.eqtls<-fread("auxiliary_files/ADLoci29_genes.txt")

# FIMO results
motif.dat<-fread("egene_FIMO/overlapped_tcpy_egene_motif.txt")
motif.dat$TF<-gsub("_HUMAN.*","",motif.dat$motif_id)

#identify motifs that overlap the variant of interest
motif.dat.overlap<-motif.dat[motif.dat$start<51&
                               motif.dat$stop>=51
]

#previously published eQTL data downloaded from 
#https://www.ebi.ac.uk/gwas/publications/30617256
jansen.eqtls<-fread("Alzheimers_Jansen_2018.txt.gz")
colnames(jansen.eqtls)[3]<-"start"
jansen.eqtls$end<-jansen.eqtls$start

## ChIP-seq regions downloaded from the UCSC genome browser
conserved.sites<-fread("encRegTfbsClustered.txt.gz")
colnames(conserved.sites)[c(2,3,4)]<-c("chr","start","end")
conserved.sites$chr<-gsub("chr","",conserved.sites$chr)

#get eQTL associated with HLA-DRB1
interest.eqt<-eqt.valid[eqt.valid$Egene3%in%c("HLA-DRB1")]
interest.eqt$start<-as.numeric(gsub("^[0-9A-Za-z]+:","",interest.eqt$`Chr:bp1`))
interest.eqt$end<-interest.eqt$start
interest.eqt$chr<-str_extract(interest.eqt$`Chr:bp1`,"^[0-9A-Za-z]+")

#allow for a 500KB distance for overlap on either direction, so an overall
#1MB overlap
eqt.overlaps<-determine_anchor_overlap(interest.eqt,jansen.eqtls[,c("chr","start","end"),
                                                                 with=FALSE],
                                       max_gap=5e5)
colnames(jansen.eqtls)[colnames(jansen.eqtls)%in%c("pvalue")]<-"P"
colnames(jansen.eqtls)[colnames(jansen.eqtls)%in%c("rsid")]<-"EQTL3"
jansen.eqtls$label<-""

lead.snp<-phys.eqtls[grep("HLA-DRB1",phys.eqtls$`Genes_InLD.R2>0.4`),]

jansen.eqtls$label[jansen.eqtls$EQTL3%in%interest.eqt$"EQTL3"|
                     jansen.eqtls$EQTL3%in%lead.snp$SNP]<-"red"


snp.overlap<-jansen.eqtls[eqt.overlaps$rep2_idx,]

snp.overlap$start<-snp.overlap$start/1e6

phys.eqtls$chr<-paste0("chr",phys.eqtls$Chr)
phys.eqtls$start<-phys.eqtls$bp
phys.eqtls$end<-phys.eqtls$bp

conserved.overlaps<-determine_anchor_overlap(phys.eqtls[,c("Chr","start","end")],
                                             conserved.sites,
                                             max_gap=-1)

conserved.sites$eGene<-""
conserved.sites$eGene[conserved.overlaps$rep2_idx]<-interest.eqt$Egene3[
  conserved.overlaps$rep1_idx
]

colnames(interest.eqt)[c(3,18,19)]<-c("eqt_start","eqt_end","eqt_chr")
conserved.sites$eqt_start<-""
conserved.sites$eqt_end<-""
conserved.sites$eqt_chr<-""

conserved.sites$eqt_start[conserved.overlaps$rep2_idx]<-interest.eqt$eqt_start[
  conserved.overlaps$rep1
]

conserved.sites$eqt_end[conserved.overlaps$rep2_idx]<-interest.eqt$eqt_end[
  conserved.overlaps$rep1
]

conserved.sites$eqt_chr[conserved.overlaps$rep2_idx]<-interest.eqt$eqt_chr[
  conserved.overlaps$rep1
]

overlapped.motifs<-conserved.sites[conserved.overlaps$rep2_idx,]

motif.dat.overlap$eGene<-gsub("_.*","",motif.dat.overlap$sequence_name)
overlapped.chip.motifs<-merge(overlapped.motifs,motif.dat.overlap,
                              by.x=c("V5","eGene"),by.y=c("TF","eGene"))

overlapped.chip.motifs$padj<-p.adjust(overlapped.chip.motifs$`p-value`,
                                      method="fdr")

colnames(overlapped.chip.motifs)[c(5,17,18)]<-c("start","overlap_start",
                                                "overlap_end")

overlapped.egene.motif.chip<-overlapped.chip.motifs[!duplicated(
  overlapped.chip.motifs[,
                         c("V5","eGene")])&
  overlapped.chip.motifs$eGene%in%c(
   "HLA-DRB1")
  ,]
overlapped.egene.motif.chip$overlap_start<-as.numeric(
  overlapped.egene.motif.chip$eqt_start)+
  as.numeric(overlapped.egene.motif.chip$overlap_start)-51

overlapped.egene.motif.chip$overlap_end<-as.numeric(
  overlapped.egene.motif.chip$eqt_start)+
  as.numeric(overlapped.egene.motif.chip$overlap_end)-51

motif.regions<-overlapped.egene.motif.chip[,c("V5",
                                              "overlap_start","overlap_end")]
colnames(motif.regions)<-c("tf","start","stop")
motif.regions$type<-"motif"
motif.regions$tf_factor<-c(1,2,3)

chip.regions<-overlapped.egene.motif.chip[,c("V5",
                                             "start","end")]
colnames(chip.regions)<-c("tf","start","stop")
chip.regions$type<-"peak"
chip.regions$tf_factor<-c(1,2,3)

overlapped.regions<-rbind(motif.regions,chip.regions)

# write tracks for visualization on the genome browser
write.table(overlapped.regions,"FIMO_chip_overlap_HLA_DRB.txt",
            sep='\t',quote=F,row.names=F)

no.motif.chip<-overlapped.motifs[overlapped.motifs$eGene%in%c("HLA-DRB1")&
  !overlapped.motifs$V5%in%
    overlapped.regions$tf]

chip.nomo.regions<-no.motif.chip[,c("V5","start","end")]
colnames(chip.nomo.regions)[c(1,3)]<-c("tf","stop")
chip.nomo.regions$type<-"peak"
chip.nomo.regions$tf_factor<-seq(4,nrow(chip.nomo.regions)+3)

all.overlapped.regions<-rbind(overlapped.regions,chip.nomo.regions)

upper.lim<-interest.eqt$eqt_start[interest.eqt$Egene3%in%c("HLA-DRB1")]+50
lower.lim<-interest.eqt$eqt_end[interest.eqt$Egene3%in%c("HLA-DRB1")]-50

all.overlapped.regions$start[all.overlapped.regions$start<lower.lim]<-lower.lim
all.overlapped.regions$stop[all.overlapped.regions$stop>upper.lim]<-upper.lim
write.table(all.overlapped.regions,"eqt_motif_search/motif_chip_seq_peaks_HLA_DRB.txt",
            sep='\t',quote=F,row.names=F)