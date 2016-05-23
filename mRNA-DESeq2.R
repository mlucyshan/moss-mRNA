library(methods)
library("gplots")
library("RColorBrewer")
library("DESeq2")
library("Rsamtools")
library("GenomicFeatures")
library("GenomicAlignments")

sampleTable<-read.table("/Data05/mengge/moss_mRNA_rai1/tophat_output/moss_ABA_sampleFile.txt", header=TRUE)
bamFiles<-file.path("/Data05/mengge/moss_mRNA_rai1/tophat_output", sampleTable$FileName)
hse<-makeTxDbFromGFF("/Data05/mengge/annotation/Ppatens/annotation/Ppatens_318_mRNA.gtf", format="gtf")
exonsByGene<-exonsBy(hse, by="gene")
se<-summarizeOverlaps(exonsByGene, BamFileList(bamFiles), mode="Union", singleEnd=FALSE, ignore.strand=FALSE, fragments=TRUE)
colData(se)<-DataFrame(sampleTable)
colnames(se)<-sampleTable$SampleName
dds<-DESeqDataSet(se, design =~ time)
dds<-dds[rowSums(counts(dds))>1,]
dds<-DESeq(dds,fitType="local")
res<-results(dds)
pdf("/Data05/mengge/moss_mRNA_rai1/whole_gene_count/DESeq_plots_mRNA.pdf",title="DESeq plots (mRNA)")
plotMA(res, ylim=c(-3,3))
plotDispEsts(dds, ylim=c(1e-6,1e1))
hist(res$pvalue, breaks=20, col="grey", main="histogram of pvalue")
hist(res$padj, breaks=20, col="grey", main="histogram of pvalue (adjusted)")
qs <-c(0, quantile(res$baseMean[res$baseMean > 0], 0:7/7))
bins<-cut(res$baseMean, qs)
levels(bins)<-paste("~", round(0.5*qs[-1] + 0.5*qs[-length(qs)]))
ratios<-tapply(res$pvalue, bins, function(p) mean(p<0.01, na.rm=TRUE))
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")

# 30min/00min#
res_30vs00<-results(dds, contrast=c("time", "30min", "00min"))
res_30vs00Sig<-subset(res_30vs00, res_30vs00$padj<0.1)
#60 vs 00#
res_60vs00<-results(dds, contrast=c("time", "60min", "00min"))
res_60vs00Sig<-subset(res_60vs00, res_60vs00$padj<0.1)
#60 min vs 30 min#
res_60vs30<-results(dds, contrast=c("time", "60min", "30min"))
res_60vs30Sig<-subset(res_60vs30, res_60vs30$padj<0.1)

sigGenes<-c(rownames(res_30vs00Sig), rownames(res_60vs00Sig), rownames(res_60vs30Sig))
sigGenes<-unique(sigGenes)
#construct matrix with fold change for genes that were significant in at least one comparison#
sigGenesTable<-as.matrix(cbind(res_30vs00[sigGenes, ]$log2FoldChange, res_60vs00[sigGenes, ]$log2FoldChange, res_60vs30[sigGenes, ]$log2FoldChange))
rownames(sigGenesTable)<-sigGenes
colnames(sigGenesTable)<-c("min30vsmin00", "min60vsmin00", "min60vsmin30")
sampleDists<-dist(t(assay(dds)))
sampleDistMatrix<-as.matrix(sampleDists)
rownames(sampleDistMatrix)<-colnames(dds)
colnames(sampleDistMatrix)<-NULL
colours=colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
heatmap.2(sampleDistMatrix, trace="none", col=colours)
rds<-rlog(dds) #note that PCA will only work with rlog or varianceStabilizingTransformation values#
plotPCA(rds, intgroup="time")

#PPS Clustering & Gene Expression#
PPS.cluster_1<-scan("/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_cluster_1.Ppatens_318_gene.txt", what="c")
PPS.cluster_2<-scan("/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_cluster_2.Ppatens_318_gene.txt", what="c")
PPS.cluster_3<-scan("/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_cluster_3.Ppatens_318_gene.txt", what="c")
clusters <- list(PPS.cluster_1, PPS.cluster_2, PPS.cluster_3)
count_value<-counts(dds)
count_value<-assay(se)
count_value_avg<-cbind((count_value[,1]+count_value[,4])/2, (count_value[,2]+count_value[,5])/2, (count_value[,3]+count_value[,6])/2 )
colnames(count_value_avg)<- c("min00", "min30", "min60")
ncluster<-3
cols <- brewer.pal(ncluster, "Set1")
for (i in 1:ncluster){
        tmp<-count_value_avg[unlist(clusters[i]),]
        plot(colMeans(tmp), col="black", lwd=3, type="b", ylim=c(min(tmp)*0.5, max(tmp)*1.1), main=sprintf("Expression of Genes in PPS cluster %d",
i), xlab="time point", ylab="Raw Counts") 
        for (j in 1:nrow(tmp)){
                 lines(tmp[j,], col=cols[i], lwd=0.6, type="b")
        }
        lines(colMeans(tmp), col="black", lwd=3, type="b")
}

count_value_avg_2<-sweep(count_value_avg, 1, count_value_avg[,1], "-")
for (i in 1:ncluster){
        tmp<-count_value_avg_2[unlist(clusters[i]),]
        plot(colMeans(tmp), col="black", lwd=3, type="b", ylim=c(min(tmp)*1.1, max(tmp)*1.1), main=sprintf("Expression of Genes in PPS cluster %d",
i), xlab="time point", ylab="Raw Counts (zeroed to 00min)")
        for (j in 1:nrow(tmp)){
                 lines(tmp[j,], col=cols[i], lwd=0.6, type="b")
        }
        lines(colMeans(tmp), col="black", lwd=3, type="b")
}
#Table of normalized counts#
final.Ppatens_318_mRNA_read_counts<-counts(dds, normalized=TRUE)
write.table(final.Ppatens_318_mRNA_read_counts, file="/Data05/mengge/moss_mRNA_rai1/mRNA_read_counts_norm.txt", quote=FALSE, sep="\t")
final.DE_Ppatens_318_mRNA_read_counts<-counts(dds, normalized=TRUE)[sigGenes,]
write.table(final.DE_Ppatens_318_mRNA_read_counts, file="/Data05/mengge/moss_mRNA_rai1/DESeq_sighits_mRNA_read_counts_norm.txt", quote=FALSE, sep="\t")
for (i in 1:ncluster) {
        final <- counts(dds, normalized=T)[intersect(unlist(clusters[i]),rownames(dds)),]
        write.table(final, file=sprintf("/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_cluster_%d.mRNA_Ppatens_318_mRNA_read_counts.txt", i), quote=F, sep="\t")
}
dev.off()
#Kmeans clustering on the DE mRNAs#
select<-sigGenes
count_value_centered<-(count_value_avg-count_value_avg[,1])[select,]
colnames(count_value_centered)<-c("min00", "min30", "min60")
pdf("/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_mRNA_diagnostics.pdf")
for (ncluster in 2:9){
      cols<-brewer.pal(ncluster, "Set1")
      kclus<-kmeans(as.matrix(count_value_centered[,c(2,3)]), centers=ncluster, iter.max=10000)
      plot(count_value_centered[,2], count_value_centered[,3], xlab="Average expression value (30min)", ylab="Average expression value (60min)", main=sprintf("k-means clustering (k=%d)", ncluster), type="p", col=cols[kclus$cluster])
      sig <-count_value_centered[,c(2,3)]
      for (i in 1:ncluster){
            tmp<-sig[which(kclus$cluster==i),]
            plot(c(0,kclus$centers[i,]), col="black", lwd=3, type="b", ylim=c(min(sig)*1.1, max(sig)*1.1), main=sprintf("cluster %d (%d genes)",i,kclus$size[i]), xlab="time point", ylab="Average expression" )
            for (j in 1:kclus$size[i]){
                   lines(c(0,tmp[j,]), col=cols[i], lwd=0.6, type="b")
            }
            lines(c(0,kclus$centers[i,]), col="black", lwd=3, type="b")
      }
}
ncluster<-3
cols<-brewer.pal(ncluster,"Set1")
kclus<-kmeans(as.matrix(count_value_centered[,c(2,3)]), centers=ncluster, iter.max=10000)
plot(count_value_centered[,2], count_value_centered[,3], xlab="Average expression value (30min)", ylab="Average expression value (60min)", main=sprintf("k-means clustering (k=%d)", ncluster), type="p", col=cols[kclus$cluster])
sig <-count_value_centered[,c(2,3)]
for (i in 1:ncluster){
            tmp<-sig[which(kclus$cluster==i),]
            plot(c(0,kclus$centers[i,]), col="black", lwd=3, type="b", ylim=c(min(sig)*1.1, max(sig)*1.1), main=sprintf("cluster %d (%d genes)",i,kclus$size[i]), xlab="time point", ylab="Average expression" )
            for (j in 1:kclus$size[i]){
                   lines(c(0,tmp[j,]), col=cols[i], lwd=0.6, type="b")
            }
            lines(c(0,kclus$centers[i,]), col="black", lwd=3, type="b")
}
dev.off()


final<-as.data.frame(count_value_centered[select,])
final$kclus_membership<-as.character(kclus$cluster)
final<-final[order(final$kclus_membership),]
for (i in 1:ncluster){
      tmp<-final[which(final$kclus_membership==i),]
      write.table(tmp, file=sprintf("/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_DE_cluster%d.txt", i), quote=F, sep="\t")
}

write.table(final, file="/Data05/mengge/moss_mRNA_rai1/kmeans_clustering_DE_genes_clusters.txt", sep="\t", row.names=T, quote=F)
