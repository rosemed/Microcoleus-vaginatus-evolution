library(ggplot2)
library(ggpubr)
library(ggstatsplot)
library(dplyr)
library(reshape2)
library(ggsci)
library(tidyr)
library(pheatmap)

data=read.table("132_ANI_16S.txt",sep = "\t")
colnames(data)=c("value","type")
p=ggplot(data,aes(x=value,fill=type))+
  geom_density(aes(y=..scaled..),alpha=.2)+ 
  theme_bw()+ theme(legend.position = c(0.25, 0.5))+
  ylab("density")+xlab("Z value")

cluster=read.table('phylocluster.txt',header = T)
anno_col=anno_row=read.table('phylocluster.txt',header = T,row.names = 1)
anno_color=list(group=c(G1="#e6550d",G2="#bcbddc",G3="#31a354",G4="#3182bd",
                         G5="#fdae6b",G6="#a1d99b",G7="#756bb1",G8="#9ecae1"))
bk <- c(0,90,95,100)

mydata=read.table("132_ANI.txt.matrix",header=T,row.names=1,sep="\t",fill = T,check.names = F)
mydata=mydata[order(match(rownames(mydata), cluster$stran)),order(match(rownames(mydata), cluster$stran))]
m=matrix("",132,132)
for (i in 1:ncol(mydata)){
  for (j in i:ncol(mydata)){
    m[j,i]=mydata[j,i]
  }
}
m=apply(m,2,as.numeric)
m[is.na(m)]<-0

colnames(m)=rownames(mydata)
rownames(m)=rownames(mydata)
p=pheatmap(m,cellwidth = 1,cellheight = 1,
           scale = "none",
           color = c("#ffffff","#afc7e8","#f19685"),
           border_color=NA,cluster_rows = F,cluster_cols = F,
           legend=F,display_numbers=F,annotation_col=anno_col,annotation_colors = anno_color,
           annotation_names_col=F,breaks=bk,show_rownames = F,show_colnames = F,annotation_legend = F)


mydata=read.table("single16s.selfblast.out",header=F,sep="\t",fill = T,check.names = F)[,1:3]
mydata = spread(mydata, V2, V3)[,-1]
rownames(mydata)=colnames(mydata)
mydata=mydata[order(match(rownames(mydata), cluster$stran)),order(match(rownames(mydata), cluster$stran))]
m=matrix("",132,132)
for (i in 1:ncol(mydata)){
  for (j in i:ncol(mydata)){
    m[i,j]=mydata[i,j]
  }
}
m=apply(m,2,as.numeric)
m[is.na(m)]<-0
colnames(m)=rownames(mydata)
rownames(m)=rownames(mydata)

p=pheatmap(m,cellwidth = 1,cellheight = 1,
         scale = "none",
         color = c("#ffffff","#afc7e8","#f19685"),
         border_color=NA,cluster_rows = F,cluster_cols = F,
         legend=F,display_numbers=F,annotation_row=anno_col,annotation_colors = anno_color,
         annotation_names_row=F,breaks=bk,show_rownames = F,show_colnames = F,annotation_legend = F)
