#libaries
library(gplots)
library(gProfileR)
library(ggplot2)
library(ggpubr)

#reading expression data
data<-read.delim("Statistical_scores_with_coordinates.txt", header=T, sep="\t")
#Chromosome, Start, End, Strand, Gene.stable.ID.version, Gene, Day1_log2FC, Day1_padj,
#Day2_log2FC, Day2_padj, Day4_log2FC, Day4_padj, Day8_log2FC, Day8_padj, XEN_log2FC, XEN_padj, epi_log2FC, epi_padj

# Getting Differentially Expressed Genes
getdeg<-function(data,logfc,pvalue){which(abs(data[,2])>=logfc & data[,3]<=pvalue)->i; return(i)}
tgdegnames<-data[0,1]

# all degs
degs<-0;
for(i in c(7,9,11,13)){
  x<-i
  y<-i+1
  degs<-c(degs,getdeg(data[,c(6,x,y)], 1, 0.05))
}
degs<-unique(degs)
degs<-degs[2:length(degs)]
degnames<-data[degs,6]
degdata<-data[degs,]

totaldegs<-degs
totaldegnames<-degnames

# Creation of separate DEG lists
x<-7; y<-8
d1degs<-getdeg(data[,c(6,x,y)], 1, 0.05)
x<-9; y<-10
d2degs<-getdeg(data[,c(6,x,y)], 1, 0.05)
x<-11; y<-12
d4degs<-getdeg(data[,c(6,x,y)], 1, 0.05)
x<-13; y<-14
d8degs<-getdeg(data[,c(6,x,y)], 1, 0.05)


# cap log2FC values to [-3, 3]
for(i in c(7,9,11,13)){
  degdata[which(degdata[,i]>3),i]<-3
  degdata[which(degdata[,i]<=-3),i]<--3
}

# heatmaps
hcols<-colorRampPalette(c("steelblue4","white","firebrick4"))(256)
pcols<-c("dark red", "dark orange", "gold","limegreen", "dark green", "dodgerblue", "dark blue", "darkorchid4")
treatments<-c("d1", "d2", "d4", "d8")
m<-heatmap.2(as.matrix(degdata[,seq(from=7, to=13, by=2)]), col=hcols, labRow="", tracecol="black", vline=0, mar=c(8,5), hclustfun = function(x) hclust(x,method = 'ward.D2'), labCol = treatments, Colv = NA)
mt<-as.hclust(m$rowDendrogram);
cutree(mt, k=6)->totalcluster;
table(totalcluster)

heatmap.2(as.matrix(degdata[,seq(from=7, to=13, by=2)]), col=hcols, labRow="", tracecol="black", vline=0, mar=c(5,2), hclustfun = function(x) hclust(x,method = 'ward.D2'), density.info="none", key.xlab = "log2(FC)", labCol = treatments, RowSideColors = pcols[totalcluster], Colv=NA, trace="none", cexCol=2)

# Set Intersections
library(UpSetR)
setmat<-matrix(0, nrow=length(data[,1]), ncol=4)
setmat[d1degs,1]<-1
setmat[d2degs,2]<-1
setmat[d4degs,3]<-1
setmat[d8degs,4]<-1
setdatamat<-data.frame(Genes=data[,1], "d1"=setmat[,1], "d2"=setmat[,2], "d4"=setmat[,3], "d8"=setmat[,4])
upset(setdatamat, nsets=4, order.by = "freq", main.bar.color = "gray13", mainbar.y.label = "Common Genes", sets.bar.color = c("Dark Red", "Dark blue", "Dark green", "Dark orange"), matrix.dot.alpha = 0.5, matrix.color = "gray13", keep.order = T)

# Bean Plots
library(beanplot)
par(mfrow=c(2,3))
par(mar=c(3,5,3,1))
for (i in 1:6){
  j<-as.numeric(names(which(totalcluster==i)))
  beanplot(data[j,c(13,11,9,7)], names=rev(treatments), what=c(0,1,1,0), col=pcols[i], horizontal=T, main=paste("Cluster", i), las=1, ylim=c(-3,3));abline(v=0, lty=3, lwd=1, col="black")
}

m<-matrix(nrow=6, ncol=4)
for (i in 1:6){
  j<-as.numeric(names(which(totalcluster==i)))
  m[i,]<-colMeans(data[j,c(7,9,11,13)])
}
as.data.frame(m)
rownames(m)<-paste("Cluster", 1:6)
colnames(m)<-treatments
#

clusternames<-c("Up", "IntermediateUp", "IntermediateDown", "LateDown", "Down", "LateUp")

par(mfrow=c(2,3))
# rearrange order of clusters in plot
par(mar=c(3,5,3,1))
for (i in c(1,2,6,5,3,4)){
  j<-as.numeric(names(which(totalcluster==i)))
  beanplot(data[j,c(13,11,9,7)], names=rev(treatments), what=c(0,1,1,0), col=pcols[i], horizontal=T, main=clusternames[i], las=1, ylim=c(-3,3));abline(v=0, lty=3, lwd=1, col="black")
}

heatmap.2(m, labRow = clusternames, labCol = treatments, hclustfun = function(x) hclust(x,method = 'ward.D2'), trace="none",  cellnote = round(m, digits=2), notecex = 1.5, notecol = "black", density.info = "none", key.xlab = "log2(FC)", mar=c(5,15), col=hcols, Colv=NA, RowSideColors = pcols[1:6], cexCol=1.7, cexRow=1.7)


# Functional Enrichment Analysis
funcenrplot<-function(setofgenes, species, top, header){
  gprofiler(query=as.character(setofgenes), organism=species, significant=T, src_filter=c("KEGG","GO:BP","GO:MF","GO:CC","HP"))[,c(12,10,4,3)]->out
  out<-out[order(out$p.value, decreasing=F),]
  out<-out[!duplicated(out$term.name),]
  out$term.name<-factor(out$term.name, levels = out$term.name[order(out$p.value, decreasing=T)])
  #
  head(out,top)->out
  #  
  funcplot<-ggplot(out, aes(x=term.name, y=-log10(p.value), label=round(-log10(p.value)),5)) + 
    geom_point(stat='identity', aes(col=-log10(p.value)), size=6)  +
    geom_text(color="black", size=2) +
    scale_color_gradient(low="steelblue4", high="firebrick4")+
    labs(title=header, 
         subtitle="Functional Enrichments") + 
    ylim(0, max(-log10(out$p.value))+1) +
    coord_flip()
  return(funcplot)
}
for(i in 1:6){
  genenames<-totaldegnames[which(totalcluster==i)]
  fep<-funcenrplot(genenames, "mmusculus", 20, paste("Cluster", i))
  assign(paste0("FuncPlot", i, sep=""),fep)
}


# Functional Enrichment Analysis plots
# Up Regulated Clusters
# Late Up (Cluster 6)
genenames<-totaldegnames[which(totalcluster==6)]
LateUpClusterFuncPlot<-funcenrplot(genenames, "hsapiens", 20, "Late Up (Cluster 6)")
# Intermediate Up (Cluster 2)
genenames<-totaldegnames[which(totalcluster==2)]
IntermediateUpClusterFuncPlot<-funcenrplot(genenames, "hsapiens", 20, "Intermediate Up (Cluster 2)")
# Up (Cluster 1)
genenames<-totaldegnames[which(totalcluster==1)]
UpClusterFuncPlot<-funcenrplot(genenames, "hsapiens", 20, "Up (Cluster 1)")

# Down Regulated Clusters
# Late Down (Cluster 4)
genenames<-totaldegnames[which(totalcluster==4)]
LateDownClusterFuncPlot<-funcenrplot(genenames, "hsapiens", 20, "Late Down (Cluster 4)")
# Intermediate Down (Cluster 3)
genenames<-totaldegnames[which(totalcluster==3)]
IntermediateDownClusterFuncPlot<-funcenrplot(genenames, "hsapiens", 20, "Intermediate Down (Cluster 3)")
# Down (Cluster 5)
genenames<-totaldegnames[which(totalcluster==5)]
DownClusterFuncPlot<-funcenrplot(genenames, "hsapiens", 20, "Down (Cluster 5)")

library(gridExtra)
activation<-ggarrange(UpClusterFuncPlot, IntermediateUpClusterFuncPlot, LateUpClusterFuncPlot, nrow=3)
activation<-annotate_figure(activation, top=text_grob("Over-Expressed Clusters", face="bold", size=18))
activation
#
repression<-ggarrange(DownClusterFuncPlot, IntermediateDownClusterFuncPlot, LateDownClusterFuncPlot, nrow=3)
repression<-annotate_figure(repression, top=text_grob("Under-Expressed Clusters", face="bold", size=18))
repression


# Creation of DFDs
library(strucchange)
# order data
data<-data[order(data$Chromosome, data$Start),]

tempValues <- c(7,9,11,13)
tempFilenames <- c("d1_DFD.bed", "d2_DFD.bed", "d4_DFD.bed", "d8_DFD.bed")
for(j in 1:4){
  if (file.exists(tempFilenames[j])) {
    print("file exists")
  }
  else {
    #Chromosome, Start, End, Strand, Gene.stable.ID.version, Gene, Day1_log2FC, Day1_padj,
    #Day2_log2FC, Day2_padj, Day4_log2FC, Day4_padj, Day8_log2FC, Day8_padj, XEN_log2FC, XEN_padj, epi_log2FC, epi_padj
    value<-tempValues[j] # number of column from which to draw data
    
    dfd<-0
    for(ch in levels(data$Chromosome)){
      print(ch)
      grep(paste(ch,"$",sep=""), data[,1])->i;data1<-data[i,c(1,2,3,6,value, value+1)]
      
      # aggregate datapoints with filter
      win=50; # window size for running average
      
      #so that the filter doesnt block when few data
      if (length(i) <= win) {
        next
      }
      
      lfc<-filter(data1[,5], rep(1/win,win))[ceiling(win/2):(length(data1[,1])-floor(win/2))]
      # create time series object
      tsfc<-ts(lfc)
      # split in domains (adjust h parameter)
      breakpoints(tsfc ~ 1, h=0.05,  dynamic = FALSE, rescale = FALSE, type="fluctuation")->bp;confint(bp)->ci;plot(tsfc);lines(ci)
      domains<-matrix(0, nrow=length(bp$breakpoints)+1, ncol=5)
      domains[1,1]<-1;
      kinit<-1;
      for(i in 1:length(bp$breakpoints)){
        kfin<-bp$breakpoints[i]-1
        domains[i,2]<-data1$End[bp$breakpoints[i]-1]
        domains[i+1,1]<-data1$Start[bp$breakpoints[i]]
        domains[i,3]<-mean(data1[kinit:kfin,5], na.rm=T)
        domains[i,4]<-length(as.vector(na.omit(data1$Gene[kinit:kfin])))
        domains[i,5]<-paste(as.vector(na.omit(data1$Gene[kinit:kfin])), collapse=",")
        kinit<-bp$breakpoints[i]
      }
      domains[length(domains[,2]),2]<-data1$End[length(data1$End)]
      domains[length(domains[,2]),3]<-mean(data1[kinit:length(data1[,5]),5])
      domains[length(domains[,2]),4]<-length(as.vector(na.omit(data1$Gene[kinit:length(data1[,5])])))
      domains[length(domains[,2]),5]<-paste(as.vector(na.omit(data1$Gene[kinit:length(data1[,5])])), collapse=",")
      domains<-as.data.frame(domains)
      domains<-data.frame(ch, domains)
      colnames(domains)<-c("chrom", "start", "end", "mean_lfc", "gene_number", "gene_names")
      dfd<-rbind(dfd,domains)
    }
    dfd[2:length(dfd[,1]),]->dfd
    as.numeric(as.character(dfd$start))->dfd$start
    as.numeric(as.character(dfd$end))->dfd$end
    as.numeric(as.character(dfd$mean_lfc))->dfd$mean_lfc
    as.numeric(as.character(dfd$gene_number))->dfd$gene_number
    str(as.data.frame(dfd))
    
    write.table(dfd, tempFilenames[j], append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = F, quote = F)
  }
}

# Domainograms
par(mfrow=c(2,2))
# read chromosomal sizes from mm9_R.genome (in Supplementary Data)
hsgenome<-read.table("mm9_R.genome") 
# load domainogram function (in Supplementary Data)
source("mousegenomeplot.R")
dfd1<-read.delim("d1_DFD.bed", header=F, sep="\t")
dfd2<-read.delim("d2_DFD.bed", header=F, sep="\t")
dfd4<-read.delim("d4_DFD.bed", header=F, sep="\t")
dfd8<-read.delim("d8_DFD.bed", header=F, sep="\t")
which(abs(dfd1[,4])>=0.1)->xx; dfd1[xx,]->dfd1
which(abs(dfd2[,4])>=0.1)->xx; dfd2[xx,]->dfd2
which(abs(dfd4[,4])>=0.1)->xx; dfd4[xx,]->dfd4
which(abs(dfd8[,4])>=0.1)->xx; dfd8[xx,]->dfd8
mouse.genome.plot(hsgenome, dfd1, "Day 1")
mouse.genome.plot(hsgenome, dfd2, "Day 2")
mouse.genome.plot(hsgenome, dfd4, "Day 4")
mouse.genome.plot(hsgenome, dfd8, "Day 8")


# find overlap percents
if (file.exists("DFD_Percentages_overlaps.tsv")) {
  print("file exists")
} else {
  dataPicker <- function(index) {
    if (index == 1) {
      return(dfd1)
    }
    else if (index == 2) {
      return(dfd2)
    }
    else if (index == 3) {
      return(dfd4)
    }
    else if (index == 4) {
      return(dfd8)
    }
  }
  
  overlapMatrix<-matrix(data = c(0), nrow=4, ncol=4)
  for (k in 1:4) {
    for (l in 1:4) {
      firstBED <- dataPicker(k)
      secondBED <- dataPicker(l)
      
      for (i in 1:nrow(firstBED)) {
        for (j in 1:nrow(secondBED)) {
          if (firstBED[i, 1] == secondBED[j, 1]) {
            
            if (firstBED[i, 2] <= secondBED[j, 2] && secondBED[j, 2] <= firstBED[i, 3] && firstBED[i, 3] <= secondBED[j, 3]) {
              overlapMatrix[k, l] <- overlapMatrix[k, l] + firstBED[i, 3] - secondBED[j, 2] + 1
            }
            else if (secondBED[j, 2] <= firstBED[i, 2] && firstBED[i, 2] <= secondBED[j, 3] && secondBED[j, 3] <= firstBED[i, 3]) {
              overlapMatrix[k, l] <- overlapMatrix[k, l] + secondBED[j, 3] - firstBED[i, 2] + 1
            }
            else if (secondBED[j, 2] <= firstBED[i, 2] && firstBED[i, 2] <= firstBED[i, 3] && firstBED[i, 3] <= secondBED[j, 3]) {
              overlapMatrix[k, l] <- overlapMatrix[k, l] + firstBED[i, 3] - firstBED[i, 2] + 1
            }
            else if (firstBED[i, 2] <= secondBED[j, 2] && secondBED[j, 2] <= secondBED[j, 3] && secondBED[j, 3] <= firstBED[i, 3]) {
              overlapMatrix[k, l] <- overlapMatrix[k, l] + secondBED[j, 3] - secondBED[j, 2] + 1
            }
          }
        }
      }
      
      overlapMatrix[k, l] <- overlapMatrix[k, l]/sum(hsgenome[, 3])
    }
  }
  
  percOverlaps <- matrix(nrow = 4*4, ncol = 4)
  for (i in 1:4) {
    for (j in 1:4) {
      percOverlaps[(i-1)*4+j, 1] = as.character(treatments[i])
      percOverlaps[(i-1)*4+j, 2] = as.character(treatments[j])
      percOverlaps[(i-1)*4+j, 3] = overlapMatrix[i, j]
      percOverlaps[(i-1)*4+j, 4] = overlapMatrix[i, j]
    }
  }
  
  write.table(percOverlaps, "DFD_Percentages_overlaps.tsv", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = F, quote = F)
}


#Heatmap of Genome Coverage and Overlaps
library(viridis)
dfdGenomeOverlaps<-read.delim("DFD_Percentages_overlaps.tsv", header=F, sep="\t")
ovmat<-matrix(dfdGenomeOverlaps[,4], nrow=4, ncol=4)
heatmap.2(ovmat, trace="none", hclustfun = function(x) hclust(x,method = 'ward.D2'), cellnote = round(ovmat, digits=2), notecex = 1.2, notecol = "white", scale="none", col="viridis", labRow = dfdGenomeOverlaps[1:4,2], labCol = dfdGenomeOverlaps[1:4,2], density.info = "none", key.xlab = "Genome Coverage", key.title = "", Colv = FALSE, Rowv = FALSE)


#Expansion of DFD sizes
library("ggpubr")
dfd1<-read.delim("d1_DFD.bed", header=F, sep="\t")
dfd2<-read.delim("d2_DFD.bed", header=F, sep="\t")
dfd4<-read.delim("d4_DFD.bed", header=F, sep="\t")
dfd8<-read.delim("d8_DFD.bed", header=F, sep="\t")
x<-length(dfd1[,1])
y<-length(dfd2[,1])
z<-length(dfd4[,1])
w<-length(dfd8[,1])
df<-data.frame(rbind(dfd1,dfd2,dfd4,dfd8))
colnames(df)<-c("chrom","start","end","meanScore","noGenes", "Genes")
timepoint<-vector("character", length=length(df[,1]))
timepoint[1:x]<-"d1"
timepoint[(x+1):(x+y)]<-"d2"
timepoint[(x+y+1):(x+y+z)]<-"d4"
timepoint[(x+y+z+1):(x+y+z+w)]<-"d8"
df<-cbind(df, timepoint)
# Size of Domains
DFDSize<-df$end-df$start
df<-cbind(df, DFDSize)
#
# reorder factor levels
df$timepoint = factor(df$timepoint,levels(df$timepoint)[c(1,2,3,4)])
# filter out
df[which(abs(df$meanScore)>=0.1),]->df1
#
# Comparison of DFD Size
my_comparisons <- list( c("d1", "d2"), c("d2", "d4"), c("d4", "d8"))
pDFDSize<-ggboxplot(df1, x = "timepoint", y = "DFDSize", fill = "timepoint", title="Size distribution of DFDs", palette = c("firebrick4", "dark orange", "olivedrab", "steelblue4")
) + stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y =1.3* max(df1$DFDSize))                             #
pDFDSize



#Scores of DFDs
# Comparison of Mean Scores
pMeanScore<-ggboxplot(df1, x = "timepoint", y = "meanScore", fill = "timepoint", title="DFD Deregulation Score", palette = c("firebrick4", "dark orange", "olivedrab", "steelblue4")
) + stat_compare_means(comparisons = my_comparisons) + stat_compare_means(label.y = 2*max(df1$meanScore))                             
pMeanScore



##Topological-Functional Bipartite Networks map the progression of prolonged TNF stimulation
library(dplyr)
########
# d1
########
if (file.exists("network_d1.tsv")) {
  netd1<-read.delim("network_d1.tsv", header=T, sep="\t")
} else {
  d<-read.delim("d1_DFD.bed", header=F, sep="\t")
  limit<-0.1
  dup<-d[which(d[,4]>=limit),]
  ddown<-d[which(d[,4]<=-limit),]
  #
  dfd<-dup; 
  m<-matrix(ncol=3, nrow=1);
  line<-matrix(ncol=3, nrow=1);
  for (i in 1:length(dfd[,1])) {
    q<-unlist(strsplit(as.character(dfd[i,6]), ","))
    out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)];
    if (length(out[,1])>0) {
      for (j in 1:length(out[,1])) {
        line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep="");
        line[1,2]<-out[j,1];
        line[1,3]<-out[j,2];
        m<-rbind(m,line)
      }
    }
  };
  m<-m[2:length(m[,1]),]
  #
  mup<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  mup[,4]<-1
  #
  dfd<-ddown; 
  m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mdown<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  -1->mdown[,4]
  #
  mm<-rbind(mup,mdown)
  # Combine with func analysis
  degnames1<-data[getdeg(data[,c(6,7,8)], 1, 0.05),6]
  allfunc1<-gprofiler(as.character(degnames1), organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,4,3)]
  allfunc1<-allfunc1[which(allfunc1[,2]>=30 & allfunc1[,2]<=1000),]
  which(mm$Function %in% allfunc1$term.name)->x
  mm[x,]->mmfilt
  netd1<-distinct(mmfilt)
  
  write.table(netd1, "network_d1.tsv", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
}
################
# d2
################
if (file.exists("network_d2.tsv")) {
  netd2<-read.delim("network_d2.tsv", header=T, sep="\t")
} else {
  d<-read.delim("d2_DFD.bed", header=F, sep="\t")
  limit<-0.1
  dup<-d[which(d[,4]>=limit),]
  ddown<-d[which(d[,4]<=-limit),]
  #
  dfd<-dup; m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mup<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  mup[,4]<-1
  #
  dfd<-ddown; m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mdown<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  -1->mdown[,4]
  #
  mm<-rbind(mup,mdown)
  # Combine with func analysis
  degnames2<-data[getdeg(data[,c(6,9,10)], 1, 0.05),6]
  allfunc2<-gprofiler(as.character(degnames2), organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,4,3)]
  allfunc2<-allfunc2[which(allfunc2[,2]>=30 & allfunc2[,2]<=1000),]
  which(mm$Function %in% allfunc2$term.name)->x
  mm[x,]->mmfilt
  netd2<-distinct(mmfilt)
  
  write.table(netd2, "network_d2.tsv", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
}
################
# d4
################
if (file.exists("network_d4.tsv")) {
  netd4<-read.delim("network_d4.tsv", header=T, sep="\t")
} else {
  d<-read.delim("d4_DFD.bed", header=F, sep="\t")
  limit<-0.1
  dup<-d[which(d[,4]>=limit),]
  ddown<-d[which(d[,4]<=-limit),]
  #
  dfd<-dup; m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mup<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  mup[,4]<-1
  #
  dfd<-ddown; m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mdown<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  -1->mdown[,4]
  #
  mm<-rbind(mup,mdown)
  # Combine with func analysis
  degnames4<-data[getdeg(data[,c(6,11,12)], 1, 0.05),6]
  allfunc4<-gprofiler(as.character(degnames4), organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,4,3)]
  allfunc4<-allfunc4[which(allfunc4[,2]>=30 & allfunc4[,2]<=1000),]
  which(mm$Function %in% allfunc4$term.name)->x
  mm[x,]->mmfilt
  netd4<-distinct(mmfilt)
  
  write.table(netd4, "network_d4.tsv", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
}
################
# d8
################
if (file.exists("network_d8.tsv")) {
  netd8<-read.delim("network_d8.tsv", header=T, sep="\t")
} else {
  d<-read.delim("d8_DFD.bed", header=F, sep="\t")
  limit<-0.1
  dup<-d[which(d[,4]>=limit),]
  ddown<-d[which(d[,4]<=-limit),]
  #
  dfd<-dup; m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mup<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  mup[,4]<-1
  #
  dfd<-ddown; m<-matrix(ncol=3, nrow=1);line<-matrix(ncol=3, nrow=1);for(i in 1:length(dfd[,1])){q<-unlist(strsplit(as.character(dfd[i,6]), ",")); out<-gprofiler(q, organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,3)]; if(length(out[,1])>0) {for(j in 1:length(out[,1])){line[1,1]<-paste(dfd[i,1],":",dfd[i,2],"-",dfd[i,3], sep=""); line[1,2]<-out[j,1]; line[1,3]<-out[j,2]; m<-rbind(m,line)}}};m<-m[2:length(m[,1]),]
  #
  mdown<-data.frame("DFD"=m[,1], "Function"=m[,2], "Enrichment"= m[,3], "mode"=vector(mode="numeric", length=length(m)))
  -1->mdown[,4]
  #
  mm<-rbind(mup,mdown)
  # Combine with func analysis
  degnames8<-data[getdeg(data[,c(6,13,14)], 1, 0.05),6]
  allfunc8<-gprofiler(as.character(degnames8), organism="mmusculus", src_filter=c("GO:BP","GO:MF", "GO:CC", "kegg", "tf", "hp"))[,c(12,4,3)]
  allfunc8<-allfunc8[which(allfunc8[,2]>=30 & allfunc8[,2]<=1000),]
  which(mm$Function %in% allfunc8$term.name)->x
  mm[x,]->mmfilt
  netd8<-distinct(mmfilt)
  
  write.table(netd8, "network_d8.tsv", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
}
#
netd1$Enrichment<-as.numeric(as.character(netd1$Enrichment))
netd2$Enrichment<-as.numeric(as.character(netd2$Enrichment))
netd4$Enrichment<-as.numeric(as.character(netd4$Enrichment))
netd8$Enrichment<-as.numeric(as.character(netd8$Enrichment))

# Plotting networks
library(igraph)
binet<-netd1 # or netd2 or netd4 or netd8
g<-graph_from_edgelist(as.matrix(binet[,1:2]), directed = FALSE)
#
g1<-simplify(g, remove.multiple = TRUE, remove.loops = TRUE)
#
igraph::degree(g1)->dg1
#
x<-grep("chr", names(g1[,1]))
type<-vector("character", length=length(g1[,1]))
type[x]<-"region"
type[-x]<-"function"
V(g1)$type<-type
#
dg1<-log2(dg1+1)
cl1<-cluster_edge_betweenness(g1)
plot(cl1, g1, layout=layout_with_fr, vertex.size = 3*dg1, edge.arrow.size = 0.8, edge.color = "dark grey", vertex.label.color= "black",vertex.color=membership(cl1), vertex.frame.color="dark grey", label.dist=10, norm=T, main=paste("d1: Modularity:", round(modularity(cl1), digits=3)), cex.main=5, asp=0)

plot(cl1, g1, layout=layout_with_fr, vertex.label=NA, vertex.size = 3*dg1, edge.arrow.size = 0.8, edge.color = "dark grey", vertex.label.color= "black",vertex.color=membership(cl1), vertex.frame.color="dark grey", label.dist=10, norm=T, main=paste("d1: Modularity:", round(modularity(cl1), digits=3)), cex.main=5, asp=0)



#Numbers of functions in Bipartite Networks
library(viridis)
# re-factor all networks
netd1$Function<-factor(netd1$Function)
netd2$Function<-factor(netd2$Function)
netd4$Function<-factor(netd4$Function)
netd8$Function<-factor(netd8$Function)
#
netd1$DFD<-factor(netd1$DFD)
netd2$DFD<-factor(netd2$DFD)
netd4$DFD<-factor(netd4$DFD)
netd8$DFD<-factor(netd8$DFD)
#
union(union(union(netd1[,2], netd2[,2]),netd4[,2]),netd8[,2])->comfunc
mf<-matrix(0, nrow = length(comfunc), ncol=4)
#
nodd1<-aggregate(netd1$Enrichment, list(netd1$Function), length)
nodd2<-aggregate(netd2$Enrichment, list(netd2$Function), length)
nodd4<-aggregate(netd4$Enrichment, list(netd4$Function), length)
nodd8<-aggregate(netd8$Enrichment, list(netd8$Function), length)
#
mf[which(comfunc %in% nodd1$Group.1),1]<-nodd1$x
mf[which(comfunc %in% nodd2$Group.1),2]<-nodd2$x
mf[which(comfunc %in% nodd4$Group.1),3]<-nodd4$x
mf[which(comfunc %in% nodd8$Group.1),4]<-nodd8$x
#
newmf <- mf[rowSums(mf) > 2, ]
newcomfunc <- comfunc[rowSums(mf) > 2]
heatmap.2(newmf, labRow = newcomfunc, labCol = c("d1", "d2", "d4", "d8"), Colv=NA, mar=c(5,20), col=viridis, trace="none", density.info = "none", hclustfun = function(x) hclust(x,method = 'ward.D2'), key.title = "Enrichment", key.xlab = "Number of DFDs")

#table of most common enrichments
mfComFunc <- data.frame(comfunc, mf)
write.table(mfComFunc, "CommonEnrichments.tsv", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T, quote = F)
