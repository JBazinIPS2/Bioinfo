SUPPA analysis gene expression
============================================

  A=OE_ASCO
  B=Col
  
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/iso_all_sample.tpm.txt X35S.ASCO_1_S7_noribo_clean.kallisto,X35S.ASCO_2_S8_noribo_clean.kallisto,X35S.ASCO_3_S9_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.tpm out/$B.tpm  -i
  
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/psy_allsamples.psi X35S.ASCO_1_S7_noribo_clean.kallisto,X35S.ASCO_2_S8_noribo_clean.kallisto,X35S.ASCO_3_S9_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.psi out/$B.psi  -e
  
  mkdir out/${A}_vs_${B}
  
  python3 ~/.local/lib/python3.4/site-packages/SUPPA/suppa.py diffSplice -m empirical -gc -i AtRTD2_all_events.ioe \
  -p out/${A}.psi out/${B}.psi -e out/${A}.tpm out/${B}.tpm -o out/${A}_vs_${B}/comp
  

  A=RNAi_ASCO
  B=Col
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/iso_all_sample.tpm.txt RNAi.ASCO_1_S4_noribo_clean.kallisto,RNAi.ASCO_2_S5_noribo_clean.kallisto,RNAi.ASCO_3_S6_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.tpm out/$B.tpm  -i
  
  A=RNAi_ASCO
  B=Col
  Rscript ~/bin/SUPPA-master/scripts/split_file.R out/psy_allsamples.psi RNAi.ASCO_1_S4_noribo_clean.kallisto,RNAi.ASCO_2_S5_noribo_clean.kallisto,RNAi.ASCO_3_S6_noribo_clean.kallisto \
  WT_1_S1_noribo_clean.kallisto,WT_2_S2_noribo_clean.kallisto,WT_3_S3_noribo_clean.kallisto out/$A.psi out/$B.psi  -e
  
  mkdir out/${A}_vs_${B}
  
  python3 ~/.local/lib/python3.4/site-packages/SUPPA/suppa.py diffSplice -m empirical -gc -i AtRTD2_all_events.ioe \
  -p out/${A}.psi out/${B}.psi -e out/${A}.tpm out/${B}.tpm -o out/${A}_vs_${B}/comp
  
  
.. {r setup_R}


library(DESeq2)
library(IsoformSwitchAnalyzeR)
library(tidyr); library(dplyr)
library(sleuth)
library(gplots)
# source("https://bioconductor.org/biocLite.R")
# biocLite("tximport")
library(tximport)
library(AnnotationDbi)
library(GenomicFeatures)


# import function
makePaddedDataFrame <- function(l,...){
  # Make padded Data.Frame
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l, na.pad ,len=maxlen),...)
}

library('VennDiagram')
library("eulerr")

as_euler_diagram <- function(rep_list){
  # Retour an eulerr obejct from a list a genes
  #
  # rep_list : a list containing the object to compare
  unique_samples <-
    unlist(rep_list) %>%
    as.vector %>%
    sort %>%
    unique
  repartition_list <- list()
  for (current_column in names(rep_list)){
    if (length(rep_list[[current_column]]) > 0){
      repartition_list[[current_column]] = unique_samples %in% rep_list[[current_column]]
    }
  }
  repartition_matrix <- data.frame(repartition_list)
  euler(repartition_matrix)
}

na.pad <- function(x,len){
  x[1:len]
}

makePaddedDataFrame <- function(l,...){
  maxlen <- max(sapply(l,length))
  data.frame(lapply(l,na.pad,len=maxlen),...)
}

# .. ..
# 
# Common annotation lists
setwd("/K/DDEVE (8095)/REGARN (8097)/JeremieB/ASCO/JB_010419_ASCO_AtRTD2/")
### import and format data ###
## event type distribution

### import and format RNAseq data ###


### list of DEG ###○

deg=read.csv("../ASCO/ASCO_RNAseq.csv")
OE_deg=deg%>%filter(X35SASCO_WT_FDR<0.01 & abs(X35SASCO_WT_log2FC)>0.75)
RNAi_deg=deg%>%filter(RNAiASCO_WT_FDR<0.01 & abs(RNAiASCO_WT_log2FC)>0.75)


deg=list(OE_deg$gene, RNAi_deg$gene); names(deg)=c("OE_ASCOvs_Col_DEG", "RNAi_ASCO_vs_Col_DEG")
vennset<- overLapper(deg, type="vennsets")


### tpm counts 

txdb <- makeTxDbFromGFF("in/database/AtRTD2_19April2016.gtf")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

files=paste0("in/kallisto/",list.files("in/kallisto/"),"/abundance.tsv")
genecounts = tximport(files, type = "kallisto", tx2gene = tx2gene,
         countsFromAbundance="scaledTPM", ignoreTxVersion = TRUE)
colnames(genecounts$counts)=gsub("noribo_clean.kallisto","",list.files("in/kallisto/"))

tpm=genecounts$counts
rowMeans(tpm[,1:3])

tpm_mean=data.frame(row.names = row.names(tpm),
                    OE=rowMeans(tpm[,1:3]),
                    RNAi=rowMeans(tpm[,4:6]),
                    WT=rowMeans(tpm[,7:9]))%>%mutate()

ggplot(tpm_mean, aes(y=log(RNAi), x=log(OE)))+geom_point()

filter=tpm_mean[rowSums(tpm_mean>50)>1,]
deg=deg[deg$gene%in%rownames(filter),]
ggplot(deg, aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC))+geom_point(color="gray", alpha=0.2)+theme_bw()+
  geom_point(data=deg%>%filter(X35SASCO_WT_FDR<0.01&abs(X35SASCO_WT_log2FC)>0.75),
mapping = aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC), color="yellow3", alpha=0.7)+
  geom_point(data=deg%>%filter(RNAiASCO_WT_FDR<0.01&abs(RNAiASCO_WT_log2FC)>0.75), 
mapping = aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC), color="blue", alpha=0.3)

  
ggplot(deg, aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC))+geom_point(color="gray", alpha=0.2)+theme_bw()+
  geom_point(data=deg%>%filter(X35SASCO_WT_FDR<0.01&abs(X35SASCO_WT_log2FC)>0.75),
             mapping = aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC), color="red", alpha=0.7)

ggplot(deg, aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC))+geom_point(color="gray", alpha=0.2)+theme_bw()+
  geom_point(data=deg%>%filter(X35SASCO_WT_FDR<0.01&abs(X35SASCO_WT_log2FC)>0.75)%>%
               filter(RNAiASCO_WT_log2FC<0.01&abs(RNAiASCO_WT_log2FC)>0.75)%>%,
             mapping = aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC), color="red", alpha=0.7)

ggplot(deg, aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC))+geom_point(color="gray", alpha=0.2)+theme_bw()+
  geom_point(data=deg%>%filter(X35SASCO_WT_FDR<0.01&abs(X35SASCO_WT_log2FC)>0.75),
             mapping = aes(y=X35SASCO_WT_log2FC, x=RNAiASCO_WT_log2FC), color="red", alpha=0.7)








### DEG and RNAser

## event type distribution

prp8=read.csv("../data revision/prp8_IR.csv")
nsr=read.table("../data revision/database_nsrRNA_seq_withAPA.txt")%>%filter(introns_genes_nsrNAA==1)
dirs=list.dirs("out")
comp=dirs[grep("vs_Col",dirs)]

data=data.frame()
for(i in comp) {
  tmp=read.delim(paste0(i,"/comp.dpsi.temp.0"))%>%mutate(Event=Event_id)%>%separate(Event_id,c("geneID", "AS_event"), sep=";")%>%
    separate(AS_event, c("AS_type", "chromosome", "exon_start", "exon_end", "strand"), sep=":")%>%
    mutate(comparison=gsub("out/","",i))
  colnames(tmp)=c("geneID","AS_type","chromosome","exon_start","exon_end","strand","dPSI","p.val","Event","comparison") 
  tmp=tmp[tmp$p.val<0.05,]
  data=rbind(data,tmp)
}

write.csv(data, "out/diff_events_all_comparison.csv")
data_events=data%>%dplyr::group_by(comparison, AS_type)%>%dplyr::summarise(AS_n=n())

### mutant effect ####

ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  ylab("number of DAS events")+ggtitle("all_comparison")+ggsave("out/DAS_events_per_comp.pdf",width = 5, height = 5)

data_events=data_events%>%filter(!AS_type %in% c("AF","AL","MX"))
ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  ylab("number of DAS events")+ggtitle("all_comparison")+ggsave("out/DAS_events_per_comp_noAFALMX.pdf",width = 5, height = 5)

data_unique=data[!duplicated(data$geneID),]
data_events=data_unique%>%dplyr::group_by(comparison, AS_type)%>%dplyr::summarise(AS_n=n())

ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  ylab("number of DAS genes")+ggtitle("all_comparison")+ggsave("out/DAS_genes_per_comp.pdf",width = 5, height = 5)


#########################################################################################################
################## VENN DIAGRAMMS #######################################################################
#########################################################################################################
library(systemPipeR)


my_list=split(as.character(data$geneID), f = data$comparison)
my_list_events=split(as.character(data$Event), f = data$comparison)

### list of DAS ###○
das=my_list
names(das)=c("OE_ASCOvs_Col_DAS", "RNAi_ASCO_vs_Col_DAS")
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_35S_RNAi.pdf")
vennPlot(vennset, mymain="DAS_genes", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()


### list of DEG ###○

deg=read.csv("../ASCO/ASCO_RNAseq.csv")
OE_deg=deg%>%filter(X35SASCO_WT_FDR<0.01 & abs(X35SASCO_WT_log2FC)>0.75)
RNAi_deg=deg%>%filter(RNAiASCO_WT_FDR<0.01 & abs(RNAiASCO_WT_log2FC)>0.75)


deg=list(OE_deg$gene, RNAi_deg$gene); names(deg)=c("OE_ASCOvs_Col_DEG", "RNAi_ASCO_vs_Col_DEG")
vennset<- overLapper(deg, type="vennsets")

pdf("out/intersect_35S_RNAi_deg.pdf")
vennPlot(vennset, mymain="DEG", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()

pdf("out/intersect_35S_RNAi_das_deg.pdf")

all=das
all$OE_ASCOvs_Col_DEG=OE_deg$gene
all$RNAi_ASCO_vs_Col_DEG=RNAi_deg$gene

vennset<- overLapper(all,type="vennsets")

vennPlot(vennset, mymain="DEG", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()

mat=data.frame(row.names = deg$gene, RNAi=deg$RNAiASCO_WT_log2FC, OE=deg$X35SASCO_WT_log2FC)
mat=mat[c(OE_deg$gene,RNAi_deg$gene),]
pheatmap::pheatmap(mat,scale = "row")


### list of IR prp8, nsr, RNAi, 35S ###○
data=filter(data, AS_type=="RI")
my_list=split(as.character(data$geneID), f = data$comparison)
das=my_list
das$prp8=prp8$GeneID
das$nsr=nsr$id
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_prp8_nsr_35S_RNAi.pdf")
vennPlot(vennset, mymain="DAS_genes", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()

### list of IR prp8, nsr, RNAi, 35S ###○
data=filter(data, AS_type=="RI")%>%filter(comparison=="RNAi_ASCO_vs_Col")
my_list=split(as.character(data$geneID), f = data$comparison)
das=my_list
das$prp8=prp8$GeneID
das$nsr=nsr$id
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_prp8_nsr_RNAi.pdf")
vennPlot(vennset, mymain="DAS_genes", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()


### list of DAS events ###○
das=my_list_events
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_35S_RNAi_events.pdf")
vennPlot(vennset, mymain="DAS_events", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()

###############################
##### compare dPSI ###########
#############################

common=intersectmatrix(vennset)%>%as.data.frame()%>%
  add_rownames(var="event_id")%>%filter(OE_ASCO_vs_Col==1 & RNAi_ASCO_vs_Col==1)

OE=data%>%filter(comparison =="OE_ASCO_vs_Col")%>%dplyr::select(Event,dPSI)%>%dplyr::rename(dPSI_OE=dPSI)

RNAi=data%>%filter(comparison =="RNAi_ASCO_vs_Col")%>%dplyr::select(Event,dPSI)%>%dplyr::rename(dPSI_RNAi=dPSI)

data_plot=inner_join(OE,RNAi)

ggplot(data_plot, aes(x=dPSI_OE,y=dPSI_RNAi))+geom_point()+ggsave("out/dPSI_correlation.pdf")



################################
### higher stringency #########
################################"

data=data.frame()
for(i in comp) {
  tmp=read.delim(paste0(i,"/comp.dpsi.temp.0"))%>%mutate(Event=Event_id)%>%separate(Event_id,c("geneID", "AS_event"), sep=";")%>%
    separate(AS_event, c("AS_type", "chromosome", "exon_start", "exon_end", "strand"), sep=":")%>%
    mutate(comparison=gsub("out/","",i))
  colnames(tmp)=c("geneID","AS_type","chromosome","exon_start","exon_end","strand","dPSI","p.val","Event","comparison") 
  tmp=tmp[tmp$p.val<0.01 & abs(tmp$dPSI)>0.1,]
  data=rbind(data,tmp)
}

write.csv(data, "out/diff_events_all_comparison_high_conf.csv")
data_events=data%>%dplyr::group_by(comparison, AS_type)%>%dplyr::summarise(AS_n=n())


ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  
  ylab("number of DAS events")+ggtitle("all_comparison")+ggsave("out/DAS_events_per_comp_high_conf.pdf",width = 5, height = 5)


data_events=data_events%>%filter(!AS_type %in% c("AF","AL","MX"))
ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  ylab("number of DAS events")+ggtitle("all_comparison")+ggsave("out/DAS_events_per_comp_noAFALMX_high_conf.pdf",width = 5, height = 5)


data_unique=data[!duplicated(data$geneID),]
data_events=data_unique%>%dplyr::group_by(comparison, AS_type)%>%dplyr::summarise(AS_n=n())

ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  ylab("number of DAS genes")+ggtitle("all_comparison")+ggsave("out/DAS_genes_per_comp_high_conf.pdf",width = 5, height = 5)


data_events=data_unique%>%dplyr::group_by(comparison, AS_type)%>%dplyr::summarise(AS_n=n())
data_events=data_events%>%filter(!AS_type %in% c("AF","AL","MX"))
ggplot(data_events, aes(x=AS_type, y=AS_n, fill=comparison))+geom_col(position = "dodge")+
  ylab("number of DAS genes")+ggtitle("all_comparison")+ggsave("out/DAS_genes_per_compnoAFALMX_high_conf.pdf",width = 5, height = 5)



#########################################################################################################
################## VENN DIAGRAMMS #######################################################################
#########################################################################################################

library(systemPipeR)

my_list=split(as.character(data$geneID), f = data$comparison)
my_list_events=split(as.character(data$Event), f = data$comparison)

### list of DAS ###○
das=my_list
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_35S_RNAi_highconf.pdf")
vennPlot(vennset, mymain="DAS_genes", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()

### list of DAS events ###○
das=my_list_events
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_35S_RNAi_events_highconf.pdf")
vennPlot(vennset, mymain="DAS_events", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()


###############################
##### compare dPSI ###########
#############################

common=intersectmatrix(vennset)%>%as.data.frame()%>%
  add_rownames(var="event_id")%>%filter(OE_ASCO_vs_Col==1 & RNAi_ASCO_vs_Col==1)

OE=data%>%filter(comparison =="OE_ASCO_vs_Col")%>%select(Event,dPSI)%>%rename(dPSI_OE=dPSI)
RNAi=data%>%filter(comparison =="RNAi_ASCO_vs_Col")%>%select(Event,dPSI)%>%rename(dPSI_RNAi=dPSI)

data_plot=inner_join(OE,RNAi)

ggplot(data_plot, aes(x=dPSI_OE,y=dPSI_RNAi))+geom_point()+ggsave("out/dPSI_correlation_highconf.pdf")

























### list of DAS events ###○
das=my_list_events
names(das)=c("prp39_vs_col", "smd1_vs_col")
vennset<- overLapper(das, type="vennsets")

pdf("out/intersect_DAS_prp39_smd1_events_highconf.pdf")
vennPlot(vennset, mymain="DAS_genes", mysub="", colmode=2, ccol=c("red", "blue"))
dev.off()




