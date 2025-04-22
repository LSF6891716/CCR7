data <- readRDS("immune.combined.rds")
library(Seurat)
library(SeuratData)
library(viridis)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(devtools)
library(celldex)
library(nichenetr)
library(SingleR)
library(dplyr)
library(SCENIC)
library(reshape2)
library(Matrix)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ClusterGVis) 
library(profvis)
library(htmltools)
library(scRNAtoolVis)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(patchwork)
########Macrophages#####
Macrophages<-subset(data,type=="Macrophages")
DimPlot(Macrophages, label = T, cols= brewer.pal(8, "Set2"), 
        pt.size = 1.5,
        repel = T, reduction = "umap",group.by ="type")
sce=Macrophages
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce)) 
sce <- JackStraw(sce, num.replicate = 100)
sce <- ScoreJackStraw(sce, dims = 1:20)
sce <- FindNeighbors(sce, dims = 1:10)
sce <- FindClusters(sce, resolution = 0.4)
sce <- RunUMAP(sce, dims = 1:10)
sce <- RunTSNE(sce, dims = 1:10)
# sce$stim[sce$stim=="case"]<-"KO"
# sce$stim[sce$stim=="control"]<-"WT"
cell_type_cols <- c(brewer.pal(9, "Set1"),
                    "#FF34B3","#BC8F8F","#20B2AA",
                    "#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF",
                    "#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C",
                    "#8B008B","#FF3030", "#7CFC00",
                    "#000000","#708090") 
DimPlot(sce, label = T, cols= cell_type_cols, #???????????? 
        pt.size = 1.5,#?????????????????? 
        repel = T, reduction = "umap")
DimPlot(sce, cols = cell_type_cols ,reduction = "umap",group.by ="seurat_clusters",
        cols.highlight = "#AD002A",label =T,split.by = "stim",pt.size = 1.5,#?????????????????? 
        repel = T)
############Spp1+Macrophages##########
count1<-sce1@assays[["RNA"]]@counts
count1<-count1[which(count@Dimnames[[1]]=="Spp1"),]
count1<-as.data.frame(count1)
count1$custer[count1$count1>0]<-"Spp1+"
count1$custer[count1$count1==0]<-"Spp1-"
count1$GROUP<-sce@meta.data[["stim"]]
count1$
table(count1$custer,count1$GROUP)
sampleinf<-sce$seurat_clusters
stim<-sce$stim
cNMF<-sce$cNMF.cluster
cNMF<-as.data.frame(sce$cNMF.cluster)
library(ggplot2)
sample1<-as.data.frame(cbind(sampleinf,cNMF))
sample2<-as.data.frame(cbind(stim,cNMF))
colnames(sample1)<-c("Group","cNMF.program")
colnames(sample2)<-c("Group","cNMF.program")
sample<-rbind(sample1,sample2)
sample$Group[sample$Group=="case"]<-"KO"
sample$Group[sample$Group=="control"]<-"WT"
sample$Type<-NA
sample$Type[sample$Group=="KO"|sample$Group=="WT"]<-"Group"
sample$Type[sample$Group=="B30"|sample$Group=="B76"|sample$Group=="B77"]<-"KO sample"
sample$Type[sample$Group=="B2"|sample$Group=="B9"|sample$Group=="B47"]<-"WT sample"
table(sample$Type)
# sample$sampleinf<-as.factor(sample$sampleinf)
sample$Group<-factor(sample$Group,levels=c("B30","B76","B77","B2","B9","B47","KO","WT"))
sample$cNMF.program<-as.factor(sample$cNMF.program)
ggplot() + 
  geom_bar(data = sample,mapping = aes(x = Group, fill = cNMF.program), 
           position = "fill")+ 
  scale_fill_manual(values = c("#00A087FF","#E64B35FF","#4DBBD5FF","#8491B4FF","#D9D9D9"))+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))+
  facet_grid(. ~ Type,scales = "free",space= "free")+
  theme(strip.background.x = element_rect(fill =c("#C7E9C0FF")))+theme(text=element_text(size=16,  family="serif"))
#7*9\
#####DEGs#######
df <- read.csv('DEG.csv',header = T)
head(df)
pvalue = 0.05
log2FC = 1
df$group <- case_when(
  df$log2FoldChange > log2FC & df$pvalue < pvalue ~ "up",
  df$log2FoldChange < -log2FC & df$pvalue < pvalue ~ "down",
  TRUE ~ 'none'
)
head(df)
df$'-log10(pvalue)' <- -log10(df$pvalue) 
df$group <- factor(df$group, levels = c("up","down","none"))
p <- ggplot(data = df,
            aes(x = log2FoldChange, y = -log10(pvalue), color = group)) + 
  geom_point(size = 2.2) 
p
p1 <- p +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2))+
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 40), breaks = seq(0, 40, by = 10))
p1
mycol <- c("#EB4232","#2DB2EB","#d8d8d8")
mytheme <- theme_classic() +
  theme(
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
p2 <- p1 +
  scale_colour_manual(name = "", values = alpha(mycol, 0.7)) +
  mytheme
p2
p3 <- p2 +
  geom_hline(yintercept = c(-log10(pvalue)),size = 0.7,color = "black",lty = "dashed") +
  geom_vline(xintercept = c(-log2FC, log2FC),size = 0.7,color = "black",lty = "dashed") 
p3
f <- function(x){
  inputx <- seq(0.0001, x, by = 0.0001)
  y <- 1/(inputx) + (-log10(pvalue))
  dff <- rbind(data.frame(x = inputx + log2FC, y = y),
               data.frame(x = -(inputx + log2FC), y = y))
  return(dff)
}
dff_curve <- f(5)
head(dff_curve)
df$curve_y <- case_when(
  df$log2FoldChange > 0 ~ 1/(df$log2FoldChange-log2FC) + (-log10(pvalue)),
  df$log2FoldChange <= 0 ~ 1/(-df$log2FoldChange-log2FC) + (-log10(pvalue))
)
df$group2 <- case_when(
  df$`-log10(pvalue)` > df$curve_y & df$log2FoldChange >= log2FC ~ 'up',
  df$`-log10(pvalue)` > df$curve_y & df$log2FoldChange <= -log2FC ~ 'down',
  TRUE ~ 'none'
)

df$group2 <- factor(df$group2, levels = c("up","down","none"))
mycol2 <- c("#F8B606","#4A1985","#d8d8d8")
p4 <- ggplot(data = df,
             aes(x = log2FoldChange, y = -log10(pvalue), color = group2)) + 
  geom_point(size = 2.2) +
  scale_x_continuous(limits = c(-4, 4), breaks = seq(-4, 4, by = 2)) +
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 40), breaks = seq(0, 40, by = 10)) +
  scale_colour_manual(name = "", values = alpha(mycol2, 0.7)) +
  geom_line(data = dff_curve,
            aes(x = x, y = y), 
            color = "black",lty = "dashed", size = 0.7) +
  mytheme
p4
p3 + p4
top30 <- filter(df, group2 != "none") %>%
  distinct(Symbol, .keep_all = T) %>%
  top_n(30, abs(log2FoldChange))
for_label <- df %>% 
  filter(Symbol == c("PTGS2", "ALOX5AP"))

p5 <- p3 +
  geom_text_repel(data = for_label,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)
                  
p5
p6 <- p4 +
  geom_text_repel(data = for_label,
                  aes(x = log2FoldChange, y = -log10(pvalue), label = Symbol),
                  force = 80, color = 'black', size = 3.2,
                  point.padding = 0.5, hjust = 0.5,
                  arrow = arrow(length = unit(0.02, "npc"),
                                type = "open", ends = "last"),
                  segment.color="black",
                  segment.size = 0.3,
                  nudge_x = 0,
                  nudge_y = 1)
p6
df <- read.csv('diff_Macrophages.csv',header = T,check.names = F)
head(df)
pvalue = 0.05
log2FC = 1
df$group <- case_when(
  df$`log2(fc)` > log2FC & df$PValue < pvalue ~ 'up',
  df$`log2(fc)` < -log2FC & df$PValue < pvalue ~ 'down',
  TRUE ~ 'none'
)
table(df$group)
df$group <- factor(df$group, levels = c('up','none','down'))
mytheme <- theme_classic() +
  theme(plot.title = element_text(size = 17),
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.margin = margin(15,5.5,5.5,5.5))
mycol <- c("#EB4232","grey90","#2DB2EB")
p2 <- ggplot() +
  geom_point(data = df,
             aes(x = `log2(fc)`, y = -log10(PValue), color = group), size = 1.6) +
  scale_colour_manual(name = '', values = alpha(mycol, 0.7)) +
  scale_x_continuous(limits = c(-2, 2),
                     breaks = seq(-12, 12, by = 4)) + #x??????
  scale_y_continuous(expand = expansion(add = c(2, 0)),
                     limits = c(0, 30),
                     breaks = seq(0, 30, by = 10)) + #y??????
  geom_hline(yintercept = c(-log10(pvalue)), size = 0.5, color = "black", lty = 'dashed') + #ˮƽ??ֵ??
  geom_vline(xintercept = c(-log2FC, log2FC), size = 0.5, color = "black", lty = 'dashed') + #??ֱ??ֵ??
  labs(title = 'Volcano_Plot') +
  mytheme
p2
p
up <- filter(df, group == 'up') %>% distinct(Symbol, .keep_all = T) %>%
  top_n(10, -log10(PValue))

down <- filter(df, group == 'down') %>% distinct(Symbol, .keep_all = T) %>%
  top_n(10, -log10(PValue))
sig <- rbind(up,down)
p22 <- p2 +
  geom_point(data = sig,
             aes(x = `log2(fc)`, y = -log10(PValue), color = group),
             size = 3.5, alpha = 0.2) +
  geom_text_repel(data = up,
                  aes(x = `log2(fc)`, y = -log10(PValue), label = Symbol),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.5,
                  nudge_x = 8 - up$`log2(fc)`,
                  direction = "y",
                  hjust = 0) +
  geom_text_repel(data = down,
                  aes(x = `log2(fc)`, y = -log10(PValue), label = Symbol),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.5,
                  nudge_x = -6 - down$`log2(fc)`,
                  direction = "y",
                  hjust = 1)
p22
p22 <- p2 +
  geom_point(data = sig,
             aes(x = `log2(fc)`, y = -log10(PValue), color = group),
             size = 3.5, alpha = 0.2) +
  geom_text_repel(data = up,
                  aes(x = `log2(fc)`, y = -log10(PValue), label = Symbol),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.5,
                  nudge_x = 8 - up$`log2(fc)`,
                  direction = "y",
                  hjust = 0) +
  geom_text_repel(data = down,
                  aes(x = `log2(fc)`, y = -log10(PValue), label = Symbol),
                  seed = 233,
                  size = 3.5,
                  color = 'black',
                  min.segment.length = 0,
                  force = 2,
                  force_pull = 2,
                  box.padding = 0.1,
                  max.overlaps = Inf,
                  segment.linetype = 3,
                  segment.color = 'black',
                  segment.alpha = 0.5,
                  nudge_x = -6 - down$`log2(fc)`,
                  direction = "y",
                  hjust = 1)
p22

##########cNMF###########
step1=function(dir_input="count_data",dir_output="res1",k=3:10,iteration=200,feature=2000){
  files=dir(dir_input)
  for (i in files) {
    
    one.file.path=paste(dir_input,"/",i,sep = "")
    prefix=strsplit(i,"\\.")[[1]][1]
    k.choices=paste(k,collapse = " ")
    junk.file=paste(dir_output,"/",prefix,"/","cnmf_tmp/",prefix,".spectra.k_*.iter_*.df.npz",sep = "")
    
    system(paste("python -W ignore cnmf.py prepare --output-dir ",dir_output," --name ",prefix," -c ",one.file.path," -k ",k.choices," --n-iter ",iteration," --total-workers 1 --numgenes ",feature,sep=""))
    system(paste("python -W ignore cnmf.py factorize --output-dir ",dir_output," --name ",prefix," --worker-index 0",sep = ""))
    system(paste("python -W ignore cnmf.py combine --output-dir ",dir_output," --name ",prefix,sep = ""))
    system(paste("rm -f ",junk.file,sep = ""))
    system(paste("MPLBACKEND='Agg' python -W ignore cnmf.py k_selection_plot --output-dir ",dir_output," --name ",prefix,sep = ""))
    
  }
  
  dir.create( paste(dir_output,"/k_selection",sep = ""))
  file.create(paste(dir_output,"/k_selection.txt",sep = ""))
  path3=paste(dir_output,"/k_selection.txt",sep = "")
  
  for (i in str_replace(files,"\\..*$","")) {
    path1=paste(dir_output,"/",i,"/",i,".k_selection.png",sep = "")
    path2=paste(dir_output,"/k_selection/",i,".k_selection.png",sep = "")
    system(paste("cp ",path1," ",path2,sep = ""))
    
    cat(i,",\n",sep = "",file=path3,append = T)
  }
  
}

step2=function(dir_input="res1",dir_output="res2",dir_count="count_data",usage_filter=0.03,top_gene=50,cor_min=0,cor_max=0.6,color=NULL,cluster_method="complete",scale_min = -2,scale_max = 2,cluster_method2="complete"){
  dir.create(dir_output)
  dirs=setdiff(dir(dir_input),c("k_selection","k_selection.txt"))
  ref.file=read.table(paste(dir_input,"/k_selection.txt",sep = ""),header = F,sep = ",",stringsAsFactors = F)
  colnames(ref.file)=c("sample","k")
  rownames(ref.file)=ref.file$sample
  
  for (i in dirs) {
    system(paste("MPLBACKEND='Agg' python -W ignore cnmf.py consensus --output-dir ",dir_input," --name ",i," --components ",ref.file[i,"k"]," --local-density-threshold 0.02 --show-clustering",sep = ""))
    
    path1=paste(dir_input,"/",i,"/",i,".usages.k_",ref.file[i,"k"],".dt_0_02.consensus.txt",sep = "")
    path2=paste(dir_input,"/",i,"/",i,".gene_spectra_score.k_",ref.file[i,"k"],".dt_0_02.txt",sep = "")
    
    if( file.exists(path1) & file.exists(path2) ){
      system(paste("cp ",path1," ",path2," ",dir_output,"/",sep = ""))
    }
  }
  ###################################################################################
  for (i in dirs) {
    usage.file=dir(dir_output,pattern = paste(i,".usages",sep = ""))
    usage.df=read.table(paste(dir_output,"/",usage.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    colnames(usage.df)=paste(i,1:dim(usage.df)[2],sep = "_")
    
    #normalize
    usage.df=usage.df / rowSums(usage.df)
    write.table(usage.df,file = paste(dir_output,"/",i,"_program.usage.norm.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
    
    #QC1
    tmpdf1=gather(usage.df,"program","ratio")
    tmpdf1%>%ggplot(aes(x=program,y=ratio))+geom_boxplot(outlier.shape = NA)+geom_jitter(color="red",alpha=0.4,width = 0.2)+
      labs(title = i)+
      theme(
        axis.text.x.bottom = element_text(angle = 45,hjust = 1),
        plot.title = element_text(hjust = 0.5,size=20)
      )
    ggsave(paste(dir_output,"/",i,"_program.usage.norm.QC.png",sep = ""),device = "png",width = 20,height = 16,units = c("cm"))
    
    #score
    score.file=dir(dir_output,pattern = paste(i,".gene_spectra_score",sep = ""))
    score.df=read.table(paste(dir_output,"/",score.file,sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    score.df=as.data.frame(t(score.df))
    colnames(score.df)=paste(i,1:dim(score.df)[2],sep = "_")
    
    topn.df=as.data.frame(matrix(nrow = top_gene,ncol = ncol(score.df)))
    colnames(topn.df)=colnames(score.df)
    
    for (k in colnames(score.df)) {
      tmpv=score.df[,k]
      names(tmpv)=rownames(score.df)
      topn.df[,k]=names(rev(tail(sort(tmpv),top_gene)))
    }
    
    #save
    write.table(topn.df, file = paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
    score.df$gene=rownames(score.df)
    write.table(score.df,file = paste(dir_output,"/",i,"_program.Zscore.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
  }
  ###################################################################################
  check.usage=data.frame()
  for (i in dirs) {
    usage.file=paste(dir_output,"/",i,"_program.usage.norm.txt",sep = "")
    usage.df=read.table(usage.file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    check.usage=rbind(check.usage,as.data.frame(colMeans(usage.df)))
  }
  colnames(check.usage)=c("mean_ratio")
  
  check.usage$sample_programs=rownames(check.usage)
  check.usage=check.usage%>%arrange(mean_ratio)
  check.usage$sample_programs=factor(check.usage$sample_programs,levels = check.usage$sample_programs)
  
  linex=sum(check.usage$mean_ratio < usage_filter)
  check.usage%>%ggplot(aes(x=sample_programs,y=mean_ratio))+geom_point()+
    geom_hline(yintercept = usage_filter,color="red")+
    geom_vline(xintercept = linex+0.5,color="red")+
    theme(
      axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)
    )
  ggsave(paste(dir_output,"/","check.usage.png",sep = ""),width = 30,height = 16,device = "png",units = "cm")
  maybe.bg=as.character(check.usage$sample_programs[check.usage$mean_ratio < usage_filter])
  ###################################################################################
  library(pheatmap)
  library(RColorBrewer)
  library(scales)
  all.score.df=data.frame()
  all.score.topn.df=data.frame()
  for (i in dirs) {
    score.file=paste(dir_output,"/",i,"_program.Zscore.txt",sep = "")
    score.df=read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
    if (i==dirs[1]) {all.score.df=score.df}
    if (i!=dirs[1]) {
      all.score.df=all.score.df%>%inner_join(score.df,by="gene")
    }
    score.topn.file=paste(dir_output,"/",i,"_program.Zscore.top",top_gene,"gene.txt",sep = "")
    score.topn.df=read.table(score.topn.file,header = T,sep = "\t",stringsAsFactors = F)
    if (i==dirs[1]) {all.score.topn.df=score.topn.df}
    if (i!=dirs[1]) {
      all.score.topn.df=cbind(all.score.topn.df,score.topn.df)
    }
  }
  rownames(all.score.df)=all.score.df$gene
  all.score.df$gene=NULL
  all.score.df=all.score.df[rowSums(is.na(all.score.df)) == 0,] #可能有空值，需要去掉
  all.score.rm.df=all.score.df[,setdiff(colnames(all.score.df),maybe.bg)] #在质控这一步检测出来的噪声
  all.score.rm.df.cor=cor(all.score.rm.df,method = "pearson")
  all.score.rm.df.cor[all.score.rm.df.cor < cor_min]=cor_min
  all.score.rm.df.cor[all.score.rm.df.cor > cor_max]=cor_max
  colanno=as.data.frame(colnames(all.score.rm.df.cor))
  colnames(colanno)="colnames"
  colanno$sample=str_replace(colanno$colnames,"_.*","")
  rownames(colanno)=colanno$colnames
  colanno$colnames=NULL
  rowanno=as.data.frame(rownames(all.score.rm.df.cor))
  colnames(rowanno)="rownames"
  rowanno$sample=str_replace(rowanno$rownames,"_.*","")
  rownames(rowanno)=rowanno$rownames
  rowanno$rownames=NULL
  if (is.null(color)){
    color_v=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(dirs))
    names(color_v)=dirs
  }else{
    color_v=color
  }
  ann_colors = list(sample = color_v)
  tmpp=pheatmap(all.score.rm.df.cor,cluster_rows = T,cluster_cols = T,
                clustering_method = cluster_method, 
                show_colnames = F,
                treeheight_row=30,treeheight_col=0,
                border_color=NA,
                annotation_row = rowanno,annotation_col = colanno,
                annotation_names_row = F,annotation_names_col = F,
                annotation_colors = ann_colors,
                color = colorRampPalette(c("white","yellow", "red","#67001F"))(50),
                fontsize_row=12,
                width = 11.5,height = 9,
                filename = paste(dir_output,"/","program_pearson_cor.",cluster_method,".heatmap.pdf",sep = "")
  )
  write.table(all.score.rm.df.cor,file = paste(dir_output,"/","cor_heatmap_data.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
  all.score.topn.rm.df=all.score.topn.df[,setdiff(colnames(all.score.topn.df),maybe.bg)]#在质控这一步检测出来的噪声
  write.table(all.score.topn.rm.df,file = paste(dir_output,"/","program_topngene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
  ###################################################################################
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(xlsx)
  hsets <- read.gmt("hallmark_cancersea.gmt")
  enrich.result=data.frame()
  pathway_v=c()
  program_v=c()
  program_topn=read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
  for (i in 1:dim(program_topn)[2]) {
    tmp <- enricher(program_topn[,i], TERM2GENE = hsets)
    if (is.null(tmp)) {
      next
    }
    tmp1=head(tmp@result)
    tmp1$program=colnames(program_topn)[i]
    rownames(tmp1)=NULL
    enrich.result=rbind(enrich.result,tmp1)
    program_v=append(program_v,colnames(program_topn)[i])
    pathway_v=append(pathway_v,paste(tmp1$Description,collapse = ","))
  }
  write.xlsx(enrich.result,file = paste(dir_output,"/","program_topngene_enrichment.xlsx",sep = ""),row.names = F)
  enrich.df=data.frame(program=program_v,pathway=pathway_v)
  enrich.df$program=factor(enrich.df$program,levels = tmpp$tree_row$labels[tmpp$tree_row$order])
  enrich.df=enrich.df%>%arrange(program)
  write.csv(enrich.df,file = paste(dir_output,"/","program_topngene_enrichment_order.csv",sep = ""),row.names = F,quote = F)
  
  ###################################################################################
  for (i in dirs) {
    one_matrix=read.table(paste(dir_count,"/",i,".count.txt",sep = ""),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
    RowSum=rowSums(one_matrix)
    one_matrix=log1p((one_matrix / RowSum) * 10000)
    one_matrix=as.data.frame(t(one_matrix))
    plot_gene=read.table(paste(dir_output,"/","program_topngene.txt",sep = ""),header = T,sep = "\t",stringsAsFactors = F)
    plot_gene=plot_gene[,str_detect(colnames(plot_gene),i)]
    plot_gene=gather(plot_gene,key = "program",value = "gene")
    plot_gene=plot_gene%>%arrange(program)
    one_matrix=one_matrix[plot_gene$gene,]
    one_matrix=t(scale(t(one_matrix)))
    one_matrix[one_matrix < scale_min] = scale_min
    one_matrix[one_matrix > scale_max] = scale_max
    if(length(unique(plot_gene$gene)) < length(plot_gene$gene)){
      plot_gene$gene=rownames(one_matrix)
    }
    pn=length(unique(plot_gene$program))
    tmpp2=pheatmap(one_matrix,cluster_rows = F,cluster_cols = T,
                   color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
                   treeheight_col=0,
                   clustering_method = cluster_method2,
                   show_colnames=F,
                   border_color=NA,
                   gaps_row=as.numeric(cumsum(table(plot_gene$program))[- pn])
    )
    write.table(one_matrix,paste(dir_output,"/",i,"_data_heatmap.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)
    cell_sort=tmpp2$tree_col$labels[tmpp2$tree_col$order]
    gene_sort=plot_gene$gene
    matrix_new=as.data.frame(one_matrix)
    matrix_new$gene=rownames(matrix_new)
    matrix_new=matrix_new%>%reshape2::melt(id="gene")
    colnames(matrix_new)[c(2,3)]=c("cell","exp")
    matrix_new$gene=factor(matrix_new$gene,levels = rev(gene_sort))
    matrix_new$cell=factor(matrix_new$cell,levels = cell_sort)
    plot1=matrix_new%>%ggplot(aes(x=cell,y=gene,fill=exp))+geom_tile()+
      geom_hline(yintercept = as.numeric(cumsum(table(plot_gene$program))[- pn])+0.5,color="black",linetype=5)+
      labs(title = paste(i,": ",length(cell_sort)," cells; ",pn," programs",sep = ""))+
      scale_x_discrete("",expand = c(0,0))+
      scale_y_discrete("",expand = c(0,0))+
      scale_fill_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100))+
      theme_bw()+
      theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x.bottom = element_blank(),axis.text.y.left = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20),
        legend.position = "left"
      )
    dim1=5*pn
    dim2=(top_gene / 10) *2
    gene.text=t(matrix(plot_gene$gene,nrow = dim2,ncol = dim1))
    rownames(gene.text)= seq(-0.5,-(dim1-0.5),-1)
    colnames(gene.text)=seq(0.5,(dim2-0.5),1)
    gene.text=as.data.frame(gene.text)
    gene.text$dim1=rownames(gene.text)
    gene.text=reshape2::melt(gene.text,id="dim1")
    colnames(gene.text)[2:3]=c("dim2","gene")
    gene.text$dim1=as.numeric(as.character(gene.text$dim1))
    gene.text$dim2=as.numeric(as.character(gene.text$dim2))
    plot2=gene.text%>%ggplot(aes(x=dim2,y=dim1))+geom_text(aes(label=gene))+
      geom_hline(yintercept = seq(-dim1,0,5)[-c(1,length(seq(-dim1,0,5)))],color="black",linetype=5)+
      scale_x_continuous("",expand = c(0,0),limits = c(0,10))+
      scale_y_continuous("",expand = c(0,0),limits = c(-dim1,0))+
      labs(title = paste(unique(plot_gene$program),collapse = "; "))+
      theme_bw()+
      theme(
        panel.grid = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5,size = 20)
      )
    library(patchwork)
    plot3=plot1+plot2+plot_layout(widths = c(1,2))
    ggsave(filename = paste(dir_output,"/",i,"_program_gene.heatmap.pdf",sep = ""),plot = plot3,width = 46,height = 16,units = "cm")
  }
}
library(reticulate)
use_condaenv(condaenv = "cnmf_env", required = T,conda = "/home/hsy/miniconda3/bin/conda")
py_config() 
source("1.R")
step1(dir_input = "count_data",dir_output = "res1",k=3:5,iteration = 50) #这里为了演示方便，取值都比较小
source("2.R")
step2(dir_input = "res1",dir_output = "res2",dir_count = "count_data",usage_filter = 0.03,top_gene = 30,cor_min = 0,cor_max = 0.6)
library(tidyverse)
mat.files=dir("./res2/",pattern = "dt_0_02.txt$")
all.mat=data.frame()
for (fi in mat.files) {
  tmp.mat=read.table(paste0("./res2/",fi),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  tmp.mat=as.data.frame(t(tmp.mat))
  sampleid=str_replace(fi,"\\..*$","")
  colnames(tmp.mat)=paste(sampleid,colnames(tmp.mat),sep = "_")
  tmp.mat$gene=rownames(tmp.mat)
  if (sampleid == "HNSCC22") {
    all.mat=tmp.mat
  }else{
    all.mat=all.mat %>% full_join(tmp.mat,by="gene") #元素的并集进行合并
  }
}
signature.programs=c("HNSCC6_3","HNSCC22_4","HNSCC5_3")
signature.loading=all.mat[,c("gene",signature.programs)]
used.gene=c()
for (pi in signature.programs) {
  tmp.df=signature.loading[,c("gene",pi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)
signature.loading=signature.loading[signature.loading$gene %in% used.gene,]
rownames(signature.loading)=signature.loading$gene
signature.loading$gene=NULL
signature.loading[is.na(signature.loading)]<-0
signature.loading$total_loading=rowSums(signature.loading)
signature.loading$average_loading=signature.loading$total_loading / length(signature.programs)
signature.loading=signature.loading%>%arrange(desc(average_loading))
head(rownames(signature.loading),30)

########Cellchat single smaple###
library(CellChat)
library(patchwork)
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
View(cellchat_control_Inflammation)
cellchat<-cellchat_control_Inflammation
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
View(cellchat_case_pEMT)
cellchat<-cellchat_case_pEMT
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
View(cellchat_control_pEMT)
cellchat<-cellchat_control_pEMT
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
load("D:/RRRR/cellchat/cellchat_case_cell_cycle1.Rdata")
View(cellchat_case_cell_cycle)
View(cellchat_control_cell_cycle)
object.list <- list(WT = cellchat_control_cell_cycle, KO = cellchat_case_cell_cycle)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg1 + gg2
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
library(ComplexHeatmap)
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1])
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 5,  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 5,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 5,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "KO"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
pathways.show <- c("SPP1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_cell(object.list[[i]], signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network - ", names(object.list)[i]))
}
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> Plot the aggregated cell-cell communication network at the signaling pathway level
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Inflam.FIB to DC cells
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1:8), targets.use = 2, lab.cex = 0.5, title.name = paste0("Signaling from Inflam.FIB - ", names(object.list)[i]))
}
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "KO")) # set factor level
plotGeneExpression(cellchat, signaling = "CXCL", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
load("D:/RRRR/cellchat/cellchat_case_pEMT1.Rdata")
object.list <- list(WT = cellchat_control_pEMT, KO = cellchat_case_pEMT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
pathways.show <- c("SPP1")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
# Chord diagram
pathways.show <- c("SPP1")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "KO")) # set factor level
plotGeneExpression(cellchat, signaling = "SPP1", split.by = "datasets", colors.ggplot = T)
load("D:/RRRR/cellchat/cellchat_case_cell_cycle1.Rdata")
object.list <- list(WT = cellchat_control_cell_cycle, KO = cellchat_case_cell_cycle)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
patchwork::wrap_plots(plots = gg)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "functional")
#> Compute signaling network similarity for datasets 1 2
cellchat <- netEmbedding(cellchat, type = "functional")
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
cellchat <- computeNetSimilarityPairwise(cellchat, type = "structural")
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
#>   genomic data. Bioinformatics 2016.
#>
#> The new InteractiveComplexHeatmap package can directly export static
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#>
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1])
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 5,  comparison = c(1, 2), angle.x = 45)
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "KO"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in KO
net.up <- subsetCommunication(cellchat, net = net, datasets = "KO",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in WT, i.e.,downregulated in KO
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:8), targets.use = 2, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:8), targets.use = 2, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:8), targets.use = 2, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:8), targets.use = 2, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
#> If you use it in published research, please cite:
#> Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional
#>   genomic data. Bioinformatics 2016.
#>
#> The new InteractiveComplexHeatmap package can directly export static
#> complex heatmaps into an interactive Shiny app with zero effort. Have a try!
#>
#> This message can be suppressed by:
#>   suppressPackageStartupMessages(library(ComplexHeatmap))
#> ========================================
i = 1
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i])
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1])
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 2,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "KO")) # set factor level
plotGeneExpression(cellchat, signaling = "Wnt10a", split.by = "datasets", colors.ggplot = T)
library(CellChat)
library(patchwork)
pathways.show <- c("Wnt10a")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show, names(object.list)[i]))
}
i = 1
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "KO"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in KO
net.up <- subsetCommunication(cellchat, net = net, datasets = "KO",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in WT, i.e.,downregulated in KO
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = c(1:8), targets.use = 2, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = c(1:8), targets.use = 2, comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = c(1:8), targets.use = 5, slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
netVisual_chord_gene(object.list[[1]], sources.use = c(1:8), targets.use = 5, slot.name = 'net', net = net.down, lab.cex = 0.8, small.gap = 3.5, title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
pathways.show <- c("Wnt10a")
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
# Chord diagram
pathways.show <- c("Wnt10a")
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "chord", signaling.name = paste(pathways.show, names(object.list)[i]))
}
# compare all the interactions sending from fibroblast to inflamatory immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1:8), targets.use = 5,  title.name = paste0("Signaling received by Inflammation.Tumors - ", names(object.list)[i]), legend.pos.x = 10)
}
# show all the significant signaling pathways from fibroblast to immune cells
par(mfrow = c(1, 2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = c(1,2,3,4), targets.use = c(5:11),slot.name = "netP", title.name = paste0("Signaling pathways sending from fibroblast - ", names(object.list)[i]), legend.pos.x = 10)
}
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "KO")) # set factor level
plotGeneExpression(cellchat, signaling = "Wnt10a", split.by = "datasets", colors.ggplot = T)
plotGeneExpression(cellchat, signaling = "Wnt10a", split.by = "datasets", colors.ggplot = T)
load("D:/RRRR/cellchat/cellchat_control_Epi_dif.Rdata")
library(CellChat)
library(patchwork)
object.list <- list(WT = cellchat_control_Epi_dif, KO = cellchat_case_Epi_dif)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
library(ComplexHeatmap)
# combining all the identified signaling pathways from different datasets
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
# ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "all", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "OrRd")
# ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "all", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "OrRd")
# draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45)
gg1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in KO", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg1 + gg2
# define a positive dataset, i.e., the dataset with positive fold change against the other dataset
pos.dataset = "KO"
# define a char name used for storing the results of differential expression analysis
features.name = pos.dataset
# perform differential expression analysis
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets", pos.dataset = pos.dataset, features.name = features.name, only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.1, thresh.p = 1)
#> Use the joint cell labels from the merged CellChat object
# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name)
# extract the ligand-receptor pairs with upregulated ligands in KO
net.up <- subsetCommunication(cellchat, net = net, datasets = "KO",ligand.logFC = 0.2, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated recetptors in WT, i.e.,downregulated in KO
net.down <- subsetCommunication(cellchat, net = net, datasets = "WT",ligand.logFC = -0.1, receptor.logFC = -0.1)
gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = F]
gg1 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.up, sources.use = 3, targets.use = c(1:8), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Up-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
pairLR.use.down = net.down[, "interaction_name", drop = F]
gg2 <- netVisual_bubble(cellchat, pairLR.use = pairLR.use.down, sources.use = 3, targets.use = c(1:8), comparison = c(1, 2),  angle.x = 90, remove.isolate = T,title.name = paste0("Down-regulated signaling in ", names(object.list)[2]))
#> Comparing communications on a merged object
gg1 + gg2
########Cellchat case vs control###
cellchat<-cellchat_control_pEMT
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("SPP1"))
cellchat <- computeCommunProbPathway(cellchat)
#devtools::install_github("sqjin/CellChat")
library(CellChat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
pathways.show <- c("SPP1") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

pairLR.SPP1 <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.SPP1[2,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
View(cellchat_control_pEMT)
library(CellChat)
library(patchwork)
cellchat<-cellchat_control_pEMT
cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5))
df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
