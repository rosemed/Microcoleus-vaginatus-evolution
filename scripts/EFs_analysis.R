#extract environmental factors
library(geodata)
library(raster)
library(sp)
library(rgdal)

# get soil grid layers (tif files)
gph <- soil_world(var="phh2o", depth=5, path="soilgrid")
gni <- soil_world(var="nitrogen", depth=5, path="soilgrid")
gsoc <- soil_world(var="soc", depth=5, path="soilgrid")
gsand <- soil_world(var="sand", depth=5, path="soilgrid")
gsilt <- soil_world(var="silt", depth=5, path="soilgrid")
gclay <- soil_world(var="clay", depth=5, path="soilgrid")
ph <- raster(gph)
ni <- raster(gni)
soc <- raster(gsoc)
sand <- raster(gsand)
silt <- raster(gsilt)
clay <- raster(gclay)

# get WorldClim and evapo layers,
# download from https://www.worldclim.org/data/worldclim21.html and
# https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/2
data_path1<-dir("database/worldclim",pattern = "tif",full.names=TRUE)
worldclim=stack(data_path)
data_path2<-dir("database/Global Aridity and PET Database",pattern = "tif",full.names=TRUE)
evapo=raster(data_path2[1])

#change the filename as your own location file
Samples<-read.table("location.txt",row.names = 1,header =T ,sep = "\t") 
#choose worldclim as the example
EF_data<-extract(worldclim,Samples,method="bilinear")

#db-RDA
p_list = c("vegan", "ggplot2", "ggpubr", "ggrepel", "rdacca.hp", "psych", "reshape2")
for(p in p_list){if (!requireNamespace(p)){install.packages(p)}
  library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)}

set.seed(123)
dirname="db-RDA"
dir.create(dirname)
env_total=read.table('microcoleus_env0213.txt',header = T,row.names = 1)
# Perform Random Forest imputation for NA Values
imputed_data_list <- missForest(env_total)
env_total <- imputed_data_list$ximp
env=as.data.frame(rda(env_total,scale = T)$Ybar)

og=read.table('Orthogroups.GeneCount.txt',header = T,row.names = 1,sep = "\t",check.names = F)
og_dist=as.data.frame(as.matrix(vegdist(t(og), method="bray")))

pro_dist=read.table('SingleCopyGene.concatenated.phy.mldist',header = T,row.names = 1,sep = "\t",check.names = F)
gene_dist=read.table('SingleCopyGeneSeq.concatenated.phy.mldist',header = T,row.names = 1,sep = "\t",check.names = F)

ani=read.table('ANI.txt.matrix',header = T,row.names = 1,sep = "\t",fill = T,check.names = F)
ani_dist=(100-ani)/10

group <- read.table("phylocluster.txt", header = T,row.names = 1,sep="\t", comment.char="")
cap=NULL
for (dist_index in c("og","pro","gene","ani")){
if (dist_index=="ani"){dist_mat=get(paste0(dist_index,"_dist"))/10
}else{dist_mat=get(paste0(dist_index,"_dist"))}
dist_mat=dist_mat[rownames(dist_mat) %in% rownames(env),colnames(dist_mat) %in% rownames(env)]
dis=as.dist(dist_mat)
pcoa <- cmdscale(dis, k = nrow(env) - 1, eig = TRUE, add = TRUE)
pcoa_site <- pcoa$point
db_rda <- rda(pcoa_site, env, scale = FALSE)
env_fit <- envfit(db_rda, env, permutations  = 999)
r <- as.matrix(env_fit$vectors$r)
p <- as.matrix(env_fit$vectors$pvals)
env.p <- cbind(r,p)
colnames(env.p) <- c("r2","p-value")
KK <- as.data.frame(env.p)
KK$p.adj = p.adjust(KK$`p-value`, method = 'BH')
select_env=rownames(KK[KK$p.adj<0.01,])
env=env[,colnames(env) %in% select_env]
dist_mat=dist_mat[rownames(dist_mat) %in% rownames(env),colnames(dist_mat) %in% rownames(env)]
dis=as.dist(dist_mat)

pcoa <- cmdscale(dis, k = nrow(env) - 1, eig = TRUE, add = TRUE)
pcoa_site <- pcoa$point
db_rda <- rda(pcoa_site, env, scale = FALSE)
env_fit <- envfit(db_rda, env, permutations  = 999)
r <- as.matrix(env_fit$vectors$r)
p <- as.matrix(env_fit$vectors$pvals)
env.p <- cbind(r,p)
colnames(env.p) <- c("r2","p-value")
KK <- as.data.frame(env.p)
KK$p.adj = p.adjust(KK$`p-value`, method = 'BH')
write.table(KK,paste0(dirname,"/",dist_index,"_envfit.txt"),col.names = NA,sep="\t",quote = F)
cap.hp = rdacca.hp(dis, env, method = 'dbRDA', type = 'R2', scale = FALSE)
cap=rbind(cap,cbind(dist_index,cap.hp$Total_explained_variation))
write.table(cap.hp$Hier.part, paste0(dirname,"/",dist_index,"_env_effect_HP.txt"),col.names = NA,sep="\t",quote = F)

score = scores(db_rda)
CAP1 = score$sites[,1]
CAP2 = score$sites[,2]
seg = as.data.frame(db_rda$CCA$biplot)
CPA_data = as.data.frame(score$sites)
plotdata=merge(CPA_data, group, by = "row.names")
colnames(plotdata) = c('sample','CAP1','CAP2','Group')
CAP1_exp = round(db_rda$CCA$eig[1]/sum(db_rda$CCA$eig)*100,2)
CAP2_exp = round(db_rda$CCA$eig[2]/sum(db_rda$CCA$eig)*100,2)
p = ggplot(plotdata, aes(CAP1, CAP2)) +
  geom_point(aes(fill = Group, color = Group),size = 2) + 
  scale_fill_manual(values = c("#d3d3d3","#dfbf80","#f1a93b","#d96558","#b43970","#692f7c","#282a62"))+
  scale_color_manual(values = c("#d3d3d3","#dfbf80","#f1a93b","#d96558","#b43970","#692f7c","#282a62"))+
  xlab(paste('dbRDA1 ( ',CAP1_exp,'%',' )', sep = '')) + 
  ylab(paste('dbRDA2 ( ',CAP2_exp,'%',' )', sep = '')) +
  geom_segment(data = seg, aes(x = 0, y = 0, xend = seg[,1], yend = seg[,2]),
               colour = "#e26844", size = 1,
               arrow = arrow(angle = 30, length = unit(0.4, 'cm'))) +
  geom_text_repel(data = seg, segment.colour = 'black',
                  aes(x = seg[,1], y = seg[,2], 
                      label = rownames(seg)),size = 3) +
  geom_vline(aes(xintercept = 0), linetype = 'dotted') +
  geom_hline(aes(yintercept = 0), linetype = 'dotted') +
  theme_bw()+
  theme(text = element_text(family = 'sans', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
        ,legend.position = 'right'
  )+ 
  theme(axis.text = element_text(colour = 'black',size = 12))+
  stat_ellipse(level=0.95,
    linetype = 2,size=0.7,
    aes(color=Group),
    alpha=0.8)+
  coord_equal(ratio=0.9)
ggsave(paste0(dirname,"/",dist_index,"_dbRDA.pdf"), p, width=149 * 1.5, height=80 * 1.5, unit='mm')
}
write.table(cap,paste0(dirname,'/env_effect_total.txt'),row.names = F,quote=F,sep="\t")
