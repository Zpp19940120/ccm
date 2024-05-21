# ccm
cancer cluster
setwd("E:\\32.岳主任单细胞\\1.组间比较CD8T")      #设置工作目录
scRNA_raw=readRDS('CD8+-T-cells_Seurat.rds')
library(Matrix)
library(stringr)
library(Seurat)
#install.packages("Rcpp")
library(SingleR)
library(dplyr)
library(reshape2)
library(stringr)
library(RColorBrewer)
library(plyr)
library(tidyverse)
library(patchwork)
library(ggplot2)
#devtools::install_github("eddelbuettel/harmony",force = TRUE)
library(harmony)
library(devtools)
library(tidydr)
library(ggunchull)
library(ggpubr)
library(paletteer)
library(SCP)
library(BiocParallel)
library(ggsci)
library(biomaRt)
#sessionInfo()
#packageVersion("Seurat")
#devtools::install_github("zhanghao-njmu/SCP")

color1=c('#313c63','#b42e20','#ebc03e','#377b4c',
         '#7bc7cd','#5d84a4','#4B4B5D',"#EC7232")
color2=rev(brewer.pal(n=11,name="Spectral"))
color3<-c("#af2934","#ffe327","#2f4e87","#b0b9b8","#f0eedf",
          "#aed4e9","#f4a69a","#3ba889","#4593c3","#f18e0c",
          "#262a35","#c5942e","#a2a7ab")
# 将颜色向量转化为调色板函数
color_palette <- colorRampPalette(color3)
# 生成自定义数量个颜色
color4 <- color_palette(25)
##UMAP
CellDimPlot(
  srt = scRNA, group.by = c("seurat_clusters"),pt.size = 0.00001,
  palcolor = color1,raster = F,
  reduction = "umap", theme_use = "theme_blank"
)
ggsave(filename = "1.cluster.pdf",height =6 ,width = 8)
##
meta=read.csv('cluster_Anno.csv',row.names = 1)
sc_meta=scRNA@meta.data
meta=meta[rownames(sc_meta),,drop=F]
sc_meta=cbind(sc_meta,meta)
colnames(sc_meta)[ncol(sc_meta)]='celltype'
scRNA@meta.data=sc_meta
##
CellDimPlot(
  srt = scRNA, group.by = c("celltype"),pt.size = 0.00001,
  palcolor = color1,raster = F,
  reduction = "umap", theme_use = "theme_blank"
)
ggsave(filename = "1.celltype.pdf",height =6 ,width = 8)
##提取dln
id=read.table('ann.txt')[,2]
DLN=subset(scRNA,orig.ident%in%c(id))
##添加分组
ann=read.table("ann.txt",sep="\t",check.names=F)
DLN@meta.data$LBJ = "NA"
for(i in 1:nrow(ann)){
  DLN@meta.data[which(DLN@meta.data$orig.ident == ann$V2[i]),'LBJ'] <- ann$V1[i]}
unique(as.vector(DLN@meta.data$LBJ))
#
CellDimPlot(
  srt = DLN, group.by = c("celltype"),pt.size = 0.00001,
  palcolor = color1,raster = F,
  reduction = "umap", theme_use = "theme_blank"
)
ggsave(filename = "1.celltype-DLN.pdf",height =6 ,width = 8)
#####
###重新整合
expression_matrix <- scRNA_raw@assays$RNA$counts
scRNA = CreateSeuratObject(expression_matrix,min.cells = 3, min.features = 200)
scRNAsub=scRNA
meta=scRNAsub@meta.data
meta1=scRNA_raw@meta.data[,'orig.ident',drop=F]
meta=cbind(meta1,meta[,c(2,3)])
scRNAsub@meta.data=meta
scRNAsub <- NormalizeData(scRNAsub, normalization.method = "LogNormalize", scale.factor = 1e4) 
scRNAsub <- FindVariableFeatures(scRNAsub, selection.method = "vst", nfeatures = 2000)
scale.genes <- row.names(scRNAsub)
scRNAsub <- ScaleData(scRNAsub, features = scale.genes)
scRNAsub <- RunPCA(scRNAsub, features = VariableFeatures(scRNAsub))
#整合
system.time({scRNA_harmony <- RunHarmony(scRNAsub, group.by.vars = "orig.ident")})
ElbowPlot(scRNA_harmony, ndims=50, reduction="harmony") 
###
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)
scRNA1=scRNA_harmony
scRNA1 <- FindClusters(scRNA1, resolution = 0.5)
scRNA1 <- RunUMAP(scRNA1, reduction = "harmony", dims = 1:20)#zuodaozhele
scRNA1 <- RunTSNE(scRNA1, reduction = "harmony",dims = 1:20)
CellDimPlot(
  srt = scRNA1, group.by = c("seurat_clusters"),#palcolor = color1,
  raster = F,pt.size = 0.00000001,
  reduction = "umap", theme_use = "theme_blank"
)
ggsave(filename = "re_CD8_T_cluster.pdf",height =3 ,width = 4)
#save(scRNA1,file = "scRNANKT.rdata")
scRNA.markers <- FindAllMarkers(scRNA1, 
                                only.pos = TRUE,  #特异性高表达marker
                                min.pct = 0.25, 
                                logfc.threshold = 0.25
)
save(scRNA.markers,file ="scRNA.re_CD8T_markers.Rdata" )

#挑选每个细胞亚群中特意高表达的10个基因
top10 <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC) 
#整理成表格，只显示基因名字
top10_table=unstack(top10, gene ~ cluster)
names(top10_table)=gsub("X","cluster",names(top10_table))
write.csv(file="reCD8_T_top20_marker_genes.csv",top10_table,row.names=F)
##
##气泡图
Tex_genes <- c("TIGIT","LAYN","CXCL13","PDCD1","TIM3")
Naive_genes <- c("SELL","CCR7","LEF1","TCF7") 
Teff_genes=c( 'GZMK', 'CD27','CXCR5')
Tisg_genes=c('ISG15','IFI6','MX1')
memory_genes <- c("CXCR4","ZFP36L2","GPR183","IL7R")
tissue_resident_genes <- c("XCL2","XCL1","CXCR6","ID2",
                           "HOPX","CD52","ZNF683")
Terminal_Effector_Memory_genes <- c("KLRG1","CX3CR1",
                                    "ASCL2","TBX21")
Tc17_genes <- c("IL17A","IL23R","IL26","ZBTB16","RORC",
                "RORA","TMIGD2","KLRB1","SLC4A10")
Cycling_genes <- c("MKI67","STMN1")
γδT_genes<-c('TRDC','TRDV2','TRGV9')
CD4_genes<-c('CD4')
CD8_genes <- c("CD8A")
NK_genes <- c('GNLY',"NKG7","TYROBP","KLRD1","FCGR3A","NCAM1")	




# PvC Clusters Enriched in Tumor Lesions (富含于肿瘤病灶的外周血管细胞簇)
#tumor_pericytes_genes <- c("RGS5", "COL1A1", "COL3A1", "COL6A3")
#immunomodulatory_cluster_genes <- c("CCL19", "CCL21")
genes_to_check =list(
  Naive=Naive_genes,
  Tisg=Tisg_genes,
  CD4=CD4_genes,
  γδT=γδT_genes,
  Tc17=Tc17_genes,
  NK=NK_genes,
  Teff=Teff_genes,
  Temra_genes=Terminal_Effector_Memory_genes,
  memory=memory_genes,
  Tex=Tex_genes,
  tissue_resident=tissue_resident_genes,
  CD8=CD8_genes
  

)

cols=rev(brewer.pal(n=11,name="Spectral"))[6:11]

genes_to_check = lapply(genes_to_check, str_to_upper)

DotPlot(scRNA1 , 
        features = genes_to_check,
        #scale = T,
        assay='RNA',group.by = "seurat_clusters",#dot.scale = 6
)+
  scale_color_gradientn(colors = color2)+
  
  theme(legend.position = "right",legend.box = "vertical",
        axis.text.x= element_text(face = "bold.italic",colour = "#441718",size = 10,angle = 45,hjust=1),
        axis.text.y= element_text(face = "bold.italic",colour = "#441718",size = 10),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        axis.title.x =  element_blank(),
        axis.title.y = element_blank(),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        strip.text.x = element_text(face = "bold.italic",colour = "white",size = 9),
        strip.background = element_rect(fill = "#313c63"),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic",size = 10),
        panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F8F4EB"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 10)
        
  )

ggsave('re_CD8_dotplot.pdf',height = 4,width = 16)
scRNA=scRNA1
##注释
celltype=read.table("celltype_ann.txt",sep="\t",header=T,check.names=F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$id[i]),'celltype'] <- celltype$celltype[i]}
unique(as.vector(scRNA@meta.data$celltype))

##挑选想要得细胞
Idents(scRNA)="celltype"
#unique(scRNA$orig.ident)
##cluster_tsne_umap
scRNA=subset(scRNA,celltype!='CD4T')
CellDimPlot(
  srt = scRNA, group.by = c("celltype"),palcolor = color2[-6],
  raster = F,pt.size = 0.000000000001,
  reduction = "umap", theme_use = "theme_blank"
)
ggsave(filename = "re_CD8_celltype.pdf",height =4 ,width = 5)
##
##提取dln
id=read.table('ann.txt')[,2]
DLN=subset(scRNA,orig.ident%in%c(id))
##添加分组
ann=read.table("ann.txt",sep="\t",check.names=F)
DLN@meta.data$LBJ = "NA"
for(i in 1:nrow(ann)){
  DLN@meta.data[which(DLN@meta.data$orig.ident == ann$V2[i]),'LBJ'] <- ann$V1[i]}
unique(as.vector(DLN@meta.data$LBJ))
###比例柱状图
Idents(DLN)='celltype'
table(Idents(DLN), DLN$LBJ)#各样本不同细胞群细胞数
Cellratio <- prop.table(table(Idents(DLN), DLN$LBJ), margin = 2)#计算各组样本不同细胞群比例
Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio) <- c("celltype","sample","ratio") #对列名重命名
#Cellratio$sample <- factor(Cellratio$sample,levels = c())  #调整画图的x轴坐标顺序
library(ggalluvial)
#Cellratio$celltype=factor(Cellratio$celltype,levels = c('CD8.Tm','CD8.Tk','CD8.Trm','CD8.Tex','CD8.Temra','CD8.MAIT','CD8.Tisg'))
ggplot(Cellratio,aes(x=sample,y=ratio,fill=celltype,stratum=celltype,alluvium=celltype))+
  geom_col(width = 0.6,color=NA)+
  #geom_flow(width=0.4,alpha=0.2,knot.pos=0)+ # knot.pos参数可以使连线变直
  scale_fill_manual(values=color3)+
  theme_classic()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
#theme(legend.position = 'none') #取消图例展示
ggsave("re_celltype_ratio.pdf",width = 4,height = 4)

##定义比较
my_comparisons <- list(c("DLN","n-DLN"))
meta=DLN@meta.data
fq <- prop.table(table(meta$celltype, as.character(meta$orig.ident)), margin=2) *100#计算细胞比例，celltype如果你细胞类型命名是其他，请修改这里
df <- reshape2::melt(fq, value.name = "freq", varnames = c("celltype", "orig.ident"))  #宽数据转长数据
ei <- unique(meta[,colnames(meta) %in% c("orig.ident","LBJ")])
df <- merge(df, ei, by = "orig.ident",)
df$LBJ <- factor(df$LBJ,levels = c("DLN","n-DLN"))

custom_fill_colors = c('#313c63','#b42e20','#ebc03e','#377b4c',
                       '#7bc7cd','#5d84a4','#4B4B5D',"#EC7232",
                       pal_lancet(palette = c("lanonc"))(9),
                       pal_npg("nrc")(9),
                       pal_nejm("default", alpha = 0.6)(8),
                       pal_d3("category20")(20))
unique(df$celltype)

ggplot(data = df,aes(x = LBJ, #分组列名
                     y = freq, #连续变量列名
                     fill = LBJ))+ #按分组填充颜色
  scale_fill_manual(values = c("#106C61","#C31820")) + #用自定义颜色填充
  geom_violin(alpha = 0.4, position = position_dodge(width = .75), 
              size = 0.8, color="black") +
  geom_boxplot(notch = FALSE, outlier.size = -1, 
               color="black", lwd=0.8, alpha = 1) +
  geom_point(shape = 21, size=3, # 点的性状和大小
             position = position_jitterdodge(), # 让点散开
             color="black", alpha = 1) +
  #ylim(0,45)+
  #theme_classic() + 
  ylab("Proportion (%)") +
  xlab("") +
  #labs(title = "TOP2A+ neoplasm")+
  facet_wrap("celltype", scales = "free_y", ncol =3 )+
  theme(axis.text= element_text(face = "bold.italic",colour = "#441718",size = 12,hjust = 0.5),
        strip.background = element_rect(fill = "#074481"),
        strip.text = element_text(face = "bold.italic",colour = "white",size = 8),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        #axis.ticks = element_blank(),
        
        #axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
        title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        #legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F8F4EB"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        #legend.title = element_text(face ="bold.italic",size = 13)
        legend.position = "none"
        #axis.title = element_blank() ,
  )+
  # 如果不要组间比较就注释掉下面这行
  stat_compare_means(method = "wilcox.test"#,label.y = max(df$freq)-5
  )
ggsave(filename = "re_celltype_prop.pdf",height = 9,width = 9)
save(DLN,file = "DLN.rdata")
save(scRNA,file = "CD8.rdata")
##面积图
##数据错误调整
cli=DLN@meta.data
library(ggplot2)
cli=cli[,c("LBJ","celltype")]
A <- prop.table(table(cli$celltype, cli$LBJ), margin = 2)
A <- as.data.frame(A)
colnames(A) <- c("celltype", "LBJ", "Freq")
cluster_cols <- c('#313c63','#b42e20','#ebc03e','#377b4c',
                  '#7bc7cd','#5d84a4','#4B4B5D',"#EC7232",
                  pal_lancet(palette = c("lanonc"))(9),
                  pal_npg("nrc")(9),
                  pal_nejm("default", alpha = 0.6)(8),
                  pal_d3("category20")(20))[-7]
unique(A$LBJ)
ggplot(A,aes(x = LBJ,y =Freq,
             group=celltype))+
  stat_summary(geom = 'line',fun='mean',cex=1,col='white')+#先要有折线
  geom_area(data = A,aes(fill=celltype))+#折线下面积，填充用celltype
  scale_fill_manual(values=color3)+
  ggtitle("Cell percentage")+
  labs(x=NULL,y=NULL)+
  theme(axis.text.x = element_text(face = "bold.italic",colour = "#441718",size = 16,angle = 45,hjust = 1),
        axis.text.y = element_text(face = "bold.italic",colour = "#441718",size = 16),
        #panel.border = element_rect(colour = "black", fill=NA, size=0.2),
        axis.ticks = element_blank(),
        #axis.title = element_text(face = "bold.italic",colour = "#441718",size = 16),#Font face ("plain", "italic", "bold", "bold.italic")
        #axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_blank(),
        plot.title = element_text(face = "bold.italic",colour = "#441718",size = 16,hjust = 0.5),
        #title = element_text(face = "bold.italic",colour = "#441718",size = 13),
        legend.text = element_text(face ="bold.italic"),
        panel.border = element_rect(fill=NA,color="black",size=1.5,linetype="solid"),
        panel.background = element_rect(fill = "#F1F6FC"),
        panel.grid.major = element_line(color = "#CFD3D6", size =.5,linetype ="dotdash" ),
        legend.title = element_text(face ="bold.italic",size = 13)
        #legend.position = "none"
        #axis.title = element_blank() ,
  )+
  geom_vline(aes(xintercept ="DLN"),
             linetype="dashed", size=1.2, colour="white")+
  geom_vline(aes(xintercept ="n-DLN"),
             linetype="dashed", size=1.2, colour="white")
ggsave("cell_bar_CD8.pdf",height = 7,width = 6)






