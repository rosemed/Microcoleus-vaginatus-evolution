#PCA
library(ggplot2)
pca_data=function(data,group){
  pca=prcomp(data)
  pca_result<-as.data.frame(pca$x)[,1:2]
  pca_result$group=group
  pca_result$group=as.factor(pca_result$group)
  variance_prop <- pca$sdev^2 / sum(pca$sdev^2)
  pc1=round(variance_prop[1]*100,2)
  pc2=round(variance_prop[2]*100,2)
  p=ggplot(data=pca_result,aes(x=PC1,y=PC2))+
    theme_bw()+
    geom_point(aes(color = group), shape = 19, size=2)+
    theme(panel.grid = element_blank())+
    geom_vline(xintercept = 0,lty="dashed", size = 1, color = 'grey50')+
    geom_hline(yintercept = 0,lty="dashed", size = 1, color = 'grey50')+
    stat_ellipse(data=pca_result,
                 geom = "polygon",level=0.95,
                 linetype = 2,size=0.5,
                 aes(fill=group),
                 alpha=0.4)+
    theme(legend.title = element_blank(),
          panel.border = element_blank(),panel.grid.major = element_blank())+
    theme(axis.title.x=element_text(size=12),
          axis.title.y=element_text(size=12,angle=90),
          axis.text.y=element_text(size=10),
          axis.text.x=element_text(size=10),
          panel.grid=element_blank())+
    scale_color_manual(values = color)+scale_fill_manual(values = color)+
    labs(x =paste0("PC1(",pc1,"%)"),y = paste0("PC2(",pc2,"%)"))
  return(p)
}
data=read.table('genome_select_codon_fraction.txt',header = T,row.names = 1,sep = "\t")
group=data[,1]
color=c("#e6550d","#756bb1","#31a354","#bcbddc",
        "#3182bd","#fdae6b","#a1d99b","#9ecae1")
pca_codon=pca_data(data[,-1],group)

#gRodon
library(gRodon)
library(Biostrings)
library(ggridges)
library(ggplot2)
library(ggsci)
mgr_sum=NULL
for (i in dir()){
genes <- readDNAStringSet(paste0(i,"/",i,".ffn"))
CDS_IDs <- readLines(paste0(i,"/",i,"_CDS_names.txt"))
gene_IDs <- gsub(" .*","",names(genes))
genes <- genes[gene_IDs %in% CDS_IDs]
highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
mgr=predictGrowth(genes, highly_expressed)
mgr$CUBHE
mgr_sum=rbind(mgr_sum,cbind(i,mgr$CUBHE,mgr$ConsistencyHE,mgr$CUB,mgr$dCUB,mgr$d,mgr$LowerCI,mgr$LowerCI))
}
colnames(mgr_sum)=c("strain","CUBHE","ConsistencyHE","CUB","dCUB","d","LowerCI","LowerCI")
mgr_sum=as.data.frame(mgr_sum)
mgr_sum[,2:8]=sapply(mgr_sum[,2:8], as.numeric)
mgr_sum=mgr_sum[order(mgr_sum[,6]),]
write.table(mgr_sum,'../132maximal_growth_rate.txt',row.names = F,sep = "\t",quote = F)
setwd("..")
mgr_sum=read.table('132maximal_growth_rate.txt',header = T,sep = "\t")
rownames(mgr_sum)=mgr_sum$strain
group=read.table('phylocluster.txt',header = T,row.names = 1)
mgr_sum=merge(mgr_sum, group, by = "row.names")
mgr_sum$group=factor(mgr_sum$group)

color=c("#e6550d","#756bb1","#31a354","#bcbddc",
        "#3182bd","#fdae6b","#a1d99b","#9ecae1")
p=ggplot(data=mgr_sum,aes(x=d,y=group,fill=group))+
  geom_density_ridges() +
  theme_ridges() +
  theme_bw()+scale_fill_manual(values = color)+
  scale_x_continuous(expand = expansion(mult = 0))
