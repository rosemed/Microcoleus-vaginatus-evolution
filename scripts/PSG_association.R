library(tidyverse)
library(circlize)
file="132_BSM_anno.txt"
data=read.table(file,header = T,sep = "\t",quote="",comment.char = "")
xrange=table(data$group)
normal_sector_index=unique(data$group)
sector.width <- xrange[normal_sector_index] / sum(xrange[normal_sector_index])

circos.clear()
circos.par(start.degree = 90, points.overflow.warning = FALSE,
           gap.after = c(rep(2, n_distinct(data$group) - 1), 20),
           track.margin = c(0, 0.01),
           cell.padding = c(0, 0, 0, 0)
)
circos.initialize(factors = data$group, x=data$x_pos,sector.width = sector.width)

circos.track(
  factors = data$group, 
  ylim = c(0, 1), 
  track.height = 0.08, 
  bg.col = "grey90",
  panel.fun = function(x, y) {
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2]-0.7,
      CELL_META$sector.index, 
      facing = "bending.inside", 
      cex = 0.8, 
      adj = c(0.5, 0)
    )
  }
)

circos.trackPlotRegion(
  factors = data$group, 
  ylim=c(0,25), 
  track.height = 0.37, 
  panel.fun = function(x, y) {
    sector_data <- data[data$group == CELL_META$sector.index, ]
    
    for (i in 1:nrow(sector_data)) {
      x_pos <- sector_data$x_pos[i]
      y_pos <- sector_data$proportion[i]
      
      circos.lines(
        c(x_pos, x_pos), 
        c(0, y_pos), 
        col = sector_data$color[i], 
        lwd = 1
      )
      
      circos.points(
        x_pos, 
        y_pos, 
        pch = 16, 
        col = sector_data$color[i], 
        cex = 0.8
      )
    }
    

    if (CELL_META$sector.index == unique(data$group)[1]) {
      circos.yaxis(
        side = "left", 
        at = seq(0, 25, by = 5), 
        labels.cex = 0.6,
        tick.length = convert_length(2, "mm")
      )
      

      circos.text(
        CELL_META$cell.xlim[1] - convert_x(10, "mm"), 
        0.5, 
        "Selection proportion", 
        facing = "reverse.clockwise", 
        niceFacing = TRUE, 
        adj = c(0.5, 0.5), 
        cex = 0.6
      )
    }
  }
)
label_data=read.table('BSM_OG_label.txt',header=T,comment.char = "")
label_data=label_data[label_data$count >1,]
for (common_og in label_data$OG){
  common_og_pos=data[data$OG==common_og,c(6,8)]
  for (i in 1:(nrow(common_og_pos)-1)){
    for (j in (i+1):nrow(common_og_pos)){
      circos.link(common_og_pos$group[i],common_og_pos$x_pos[i],
                  common_og_pos$group[j],common_og_pos$x_pos[j],
                  col = "#4e6a74", lwd = 1)
    }
  }
}

df=NULL
for (common_og in label_data$OG){
  common_og_pos=data[data$OG==common_og,]
  for (i in 1:(nrow(common_og_pos)-1)){
    for (j in (i+1):nrow(common_og_pos)){
      df=rbind(df,cbind(common_og,common_og_pos$category[i],common_og_pos$description[i],
                        common_og_pos$group[i],common_og_pos$group[j]))
    }
  }
}
write.table(df,'common_og_anno.txt',row.names = F,quote = F,sep = "\t")

#pi value distribution
data=read.table('phylogroup_selection.location')
total_pi=NULL
total_sel_pi=NULL
for (group in unique(data$V2)){
  sub=data[data$V2==group,]
  pi=read.table(paste0(group,'_pi.windowed.pi'),header=T)
  total_pi=rbind(total_pi,data.frame(group=group,value=pi$PI))
  mid=round((sub$V6+sub$V5)/2)
  midpoint=round((pi$BIN_END+pi$BIN_START)/2)
  indices <- NULL
  for (i in seq_along(mid)) {
      diffs <- abs(midpoint - mid[i])
      closest_pos <- which.min(diffs)
      indices=c(indices,closest_pos)
  }
  indices=unique(indices)
  sub_pi=pi[indices,]
  total_sel_pi=rbind(total_sel_pi,data.frame(group=group,value=sub_pi$PI))
}
library(ggplot2)
library(ggridges)
p=ggplot() +
  geom_density_ridges(data = total_pi, aes(x = value, y = group, fill = "total"), 
                      scale = 2, alpha = 0.7) +
  geom_density_ridges(data = total_sel_pi, aes(x = value, y = group, fill = "selection"), 
                      scale = 2, alpha = 0.7) +
  theme_ridges() +
  theme_bw()+scale_fill_manual(values = c("total" = "#828d93", "selection" = "#edb176")) +
  scale_x_continuous(expand = expansion(mult = 0))
