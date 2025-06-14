library(ggplot2)
library(ggsci)
data=read.table('phylogroup_node.change',header = T,sep = "\t")
data=data[,c(1,which(colSums(data[,-1])>0)+1)]
tot=colSums(data[,-1])
merge_rows <- function(df) {
  col_sums <- colSums(df[,-1])
  
  low_percentage_rows <- apply(df[,-1], 1, function(row) all(row / col_sums < 0.1))
  
  low_percentage_df <- df[low_percentage_rows, -1]
  high_percentage_df <- df[!low_percentage_rows, -1]
  high_cate=df[!low_percentage_rows, 1]

  if (nrow(low_percentage_df) > 0) {
    merged_row <- colSums(low_percentage_df)
    result <- rbind(high_percentage_df, merged_row)
    result=data.frame(cate=c(high_cate,"other"),result)
  } else {
    result <- high_percentage_df
    result=data.frame(cate=df[,1],result)
  }
  
  return(result)
}
data <- merge_rows(data)
data$cate=factor(data$cate,levels = c("K","L","N","T","O","C","G","H","I","other","S","-"))
#color=c("#F19019","#A1D3B4","#C5AACA","#C4C4C3","#E3DDEE")
color=c("#E56F5E","#eea599","#A5BDF2","#b7cbd5","#81b3a9","#A1D3B4","#C5AACA","#E3DDEE","#FFC892","#D6CDBE","#a5a5a5","#DFE1E2")
for (i in 2:14) {
  data2=data[,c(1,i)]
  data2$node=colnames(data)[i]
  colnames(data2)[1:2]=c("category","count")
  data2$fraction=round(data2$count / tot[i-1]*100)
  data2$ymax <- cumsum(data2$fraction)
  data2$ymin <- c(0, head(data2$ymax, n=-1))
  data2$label=paste0(data2$fraction,"%")
  data2$labelPosition <- (data2$ymax + data2$ymin) / 2
  p=ggplot(data2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
    geom_rect() +
    geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
    scale_fill_manual(values=color) +
    coord_polar(theta="y") +
    xlim(c(2.8, 4)) +
    theme_void() +
    theme(legend.position = "none")+
    geom_text(x=2.8,aes(y=50,label="-"),size=3)
  assign(paste0("p",i-1),p)
}
