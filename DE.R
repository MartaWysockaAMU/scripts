#https://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day3/rnaSeq_DE.pdf

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")
BiocManager::install("Glimma")
BiocManager::install("edgeR", force=TRUE)

library(limma)
library(Glimma)
library(edgeR)
library(RColorBrewer)

setwd("/media/marta/DATA/RIBO/DE")
file = '../fc_output/mRNA_all_genes_counts.csv'
e = read.table(file,header=TRUE)
samplenames <- colnames(e)[7:15]
rnames = e[,1]
exp = e[,7:15]
rownames(exp) <- rnames
dgList <- DGEList(counts=exp, genes=rnames)
countsPerMillion <- cpm(dgList)
logCountsPerMilion <- cpm(dgList, log=TRUE)
countCheck <- countsPerMillion > 1

dgList <- DGEList(counts=countsPerMillion, genes=rnames)
groups <- c('control','control','control','vankomycin','vankomycin','vankomycin','cefoxitin','cefoxitin','cefoxitin')
dgList$samples$group <- groups

L <- mean(dgList$samples$lib.size) * 1e-6
M <- median(dgList$samples$lib.size) * 1e-6
c(L,M)

logCountsPerMilion.cutoff<- log2(10/M +2/L)
summary(countsPerMillion)
summary(logCountsPerMilion)

nsamples<-ncol(dgList)
col<-brewer.pal(nsamples,"Paired")

pdf("DE.pdf") 
par(mfrow=c(1,2))
plot(density(logCountsPerMilion[,1]),col=col[1],lwd=2,ylim=c(0,0.26),las=2,main="",xlab="")
title(main="A.Raw data",xlab="Log-cpm")
abline(v=logCountsPerMilion.cutoff,lty=3)
for(i in 2:nsamples){
  den<-density(logCountsPerMilion[,i])
  lines(den$x,den$y,col=col[i],lwd=2)
}
legend("topright", colnames(dgList),text.col=col,bty="n")

keep.exprs <- filterByExpr(dgList, group=groups)
x <- dgList[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)
xlcmp <- cpm(x, log=TRUE)
plot(density(xlcmp[,1]),col=col[1],lwd=2,ylim=c(0,0.26),las=2,main="",xlab="")
title(main="B. Filtered data",xlab="Log-cpm")
abline(v=logCountsPerMilion.cutoff,lty=3)
for(i in 2:nsamples){
  den<-density(xlcmp[,i])
  lines(den$x, den$y,col=col[i],lwd=2)
}
legend("topright", colnames(dgList),text.col=col,bty="n")

#######################################
boxplot(xlcmp,las=2,col=col,main="")
title(main="A. Unnormalised data",ylab="Log-cpm")

norm_dgList<-calcNormFactors(x)
norm_dgList$samples$norm.factors

norm_cpm<-cpm(norm_dgList)
norm_lcpm<-cpm(norm_dgList,log=TRUE)
boxplot(norm_lcpm,las=2,col=col,main="")
title(main="B. Normalised data",ylab="Log-cpm")

par(mfrow=c(1,1))
col.group<-col
levels(col.group)<-brewer.pal(nlevels(col.group),"Set1")
col.group<-as.character(col.group)
plotMDS(norm_lcpm,labels=groups,col=col.group)
legend("right",colnames(norm_lcpm),text.col='black',bty="n", fill=col, border = 'black')
title(main="Sample groups")

library(Glimma)
glMDSPlot(norm_lcpm, labels=paste(colnames(logCountsPerMilion)), groups=dgList$samples, launch=TRUE)

design<-model.matrix(~0+groups)
colnames(design)<-gsub("groups","",colnames(design))
design
contr.matrix<-makeContrasts(
  KvsC= control-cefoxitin,
  KvsV= control-vankomycin,
  CvsV= cefoxitin-vankomycin,
  levels=colnames(design))
contr.matrix

par(mfrow=c(1,2))
v <-voom(x, design,plot=TRUE)
v
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

summary(decideTests(efit))

tfit <- treat(vfit, lfc=1)
dt <- decideTests(tfit)
summary(dt)

par(mfrow=c(1,1))
vennDiagram(dt[ ,1:3], circle.col=c("turquoise", "salmon", "yellow"))

par(mfrow=c(1,3))
plotMD(tfit, column=1, status=dt[,1], main=colnames(tfit)[1], xlim=c(-8,13))
plotMD(tfit, column=1, status=dt[,2], main=colnames(tfit)[2], xlim=c(-8,13))
plotMD(tfit, column=1, status=dt[,3], main=colnames(tfit)[3], xlim=c(-8,13))

par(mfrow=c(1,1))
glMDPlot(tfit, coef=1, status=dt, main=colnames(tfit)[1], side.main="Gene", counts=norm_lcpm, groups=groups, launch=TRUE, html="MD_plot_KvsC")
glMDPlot(tfit, coef=2, status=dt, main=colnames(tfit)[2], side.main="Gene", counts=norm_lcpm, groups=groups, launch=TRUE, html="MD_plot_KvsV")
glMDPlot(tfit, coef=3, status=dt, main=colnames(tfit)[3], side.main="Gene", counts=norm_lcpm, groups=groups, launch=TRUE, html="MD_plot_CvsV")

KvsC<- topTreat(tfit, coef=1, n=Inf)
KvsV<- topTreat(tfit, coef=2, n=Inf)
CvsV<- topTreat(tfit, coef=3, n=Inf)


ga1=data.frame(genes=CvsV$genes, row.names=row.names(CvsV))
glXYPlot(CvsV$logFC,-log10(CvsV$P.Value),
         xlab="logFC",
         ylab="-log(p-val)",
         status=as.numeric(CvsV$P.Value <= 0.05),
         anno=ga1, side.main = "genes", html ="CvsV")

ga2=data.frame(genes=KvsC$genes, row.names=row.names(KvsC))
glXYPlot(KvsC$logFC,-log10(KvsC$P.Value),
         xlab="logFC",
         ylab="-log(p-val)",
         status=as.numeric(KvsC$P.Value <= 0.05),
         anno=ga2, side.main = "genes", html ="KvsC")

ga3=data.frame(genes=KvsV$genes, row.names=row.names(KvsV))
glXYPlot(KvsV$logFC,-log10(KvsV$P.Value),
         xlab="logFC",
         ylab="-log(p-val)",
         status=as.numeric(KvsV$P.Value <= 0.05),
         anno=ga3, side.main = "genes", html ="KvsV")

library(gplots)

par(mfrow=c(1,1))
mycol <- colorpanel(1000,"blue","white","red")

KvsC.topgenes <- KvsC$genes[1:100]
i <- which(v$genes$genes %in% KvsC.topgenes)
heatmap.2(norm_lcpm[i,], scale="row",
          labRow=v$genes$genes[i], labCol=groups,
          col=mycol, trace="none", density.info="none",
          margin=c(10,10), lhei=c(2,10), dendrogram="column", main = "KvsC top 100")

KvsV.topgenes <- KvsV$genes[1:100]
i <- which(v$genes$genes %in% KvsV.topgenes)
heatmap.2(norm_lcpm[i,], scale="row",
          labRow=v$genes$genes[i], labCol=groups,
          col=mycol, trace="none", density.info="none",
          margin=c(10,10), lhei=c(2,10), dendrogram="column", main = "KvsV top 100")

CvsV.topgenes <- CvsV$genes[1:100]
i <- which(v$genes$genes %in% CvsV.topgenes)
mycol <- colorpanel(1000,"blue","white","red")
heatmap.2(norm_lcpm[i,], scale="row",
          labRow=v$genes$genes[i], labCol=groups,
          col=mycol, trace="none", density.info="none",
          margin=c(10,10), lhei=c(2,10), dendrogram="column", main = "CvsV top 100")

dev.off()


KvsC_adjPval = KvsC[KvsC["adj.P.Val"] <= 0.05,]
KvsC_adjPval_down = rownames(KvsC_adjPval[KvsC_adjPval["logFC"] < 0,])
write(KvsC_adjPval_down, file="KvsC_down.csv")

KvsV_adjPval = KvsV[KvsV["adj.P.Val"] <= 0.05,]
KvsV_adjPval_down = rownames(KvsV_adjPval[KvsV_adjPval["logFC"] < 0,])
KvsV_adjPval_up = rownames(KvsV_adjPval[KvsV_adjPval["logFC"] > 0,])
write(KvsV_adjPval_down, file="KvsV_down.csv")
write(KvsV_adjPval_up, file="KvsV_up.csv")

CvsV_adjPval = CvsV[CvsV["adj.P.Val"] <= 0.05,]
CvsV_adjPval_down = rownames(CvsV_adjPval[CvsV_adjPval["logFC"] < 0,])
CvsV_adjPval_up = rownames(CvsV_adjPval[CvsV_adjPval["logFC"] > 0,])
write(CvsV_adjPval_down, file="CvsV_down.csv")
write(CvsV_adjPval_up, file="CvsV_up.csv")
