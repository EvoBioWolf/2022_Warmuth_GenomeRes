### differential expression analysis of mapped RNASeq reads ###

# install packages
BiocManager::install("apeglm")

library(Rsubread)
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)


#notes:
# - use STAR mapper and report a random valid alignment for multi-mappers; report up to 20 alignments per read
# - reporting randomly one position (in feature Counts: countMultiMappingReads=TRUE, -B=1)

# generate read count matrix - random 
cnt.r <- featureCounts(files=c("1.random.bam", "2.random.bam", ...), nthreads=4, useMetaFeatures=TRUE, annot.ext="myGTF.gtf", isGTFAnnotationFile=TRUE, allowMultiOverlap=FALSE, countMultiMappingReads=TRUE, strandSpecific=0) 
cts.r <- cnt.r$counts

# generate read count matrix - unique
cnt.u <- featureCounts(files=c("1.unique.bam", "2.unique.bam",...), nthreads=4, useMetaFeatures=TRUE, annot.ext="myGTF.gtf", isGTFAnnotationFile=TRUE, allowMultiOverlap=FALSE, countMultiMappingReads=FALSE, strandSpecific=0) 
cts.u<-cnt.u$counts

# load coldata for design
coldata<-read.csv('coldata.csv',header=T, row.names=1, stringsAsFactors=T)

# tissue + sex + species
dds.r<- DESeqDataSetFromMatrix(countData = cts.r,colData = coldata, design = ~ tissue + species + sex) 

# re-level for heatmap (DE)
dds.r$sex <- relevel(dds.r$sex, "male")

#normalised Counts
dds.r.sizeFacts<- estimateSizeFactors(dds.r)
normCounts.r <-counts(dds.r.sizeFacts, normalized=TRUE)

#run DESeq2
ddsDeseq.r <- DESeq(dds.r)

# get results - sex
ddsRes.r.sex <- results(ddsDeseq.r, contrast=c("sex","female","male"), lfcThreshold = 1)

# get results - species
ddsRes.r.spp <- results(ddsDeseq.r, contrast=c("species","cornix","corone"), lfcThreshold = 1)

# get top DE features - sex
DEfeats.sex.r.unshrunken_padj0.1 <- subset(ddsRes.r.sex, padj < 0.1)

#normalised counts for DETEs
normCounts.DETEs.sex.r<- normCounts.r.str[rownames(normCounts.r.str) %in% rownames(DEfeats.sex.r.unshrunken_padj0.1),]

# counts per tissue 
normCounts.DETEs.sex.r.gonads<-normCounts.DETEs.sex.r[, grepl( "gonads" , colnames(normCounts.DETEs.sex.r))]
normCounts.DETEs.sex.r.spleen<-normCounts.DETEs.sex.r[, grepl( "spleen" , colnames(normCounts.DETEs.sex.r))]
normCounts.DETEs.sex.r.liver<-normCounts.DETEs.sex.r[, grepl( "liver" , colnames(normCounts.DETEs.sex.r))]

	  
# estimate dispersion trend using vst instead of rlog (faster and makes fewer assumptions about the data)
vst.r<-vst(dds.r, blind=FALSE)

#calculate poisson dists
poisd <- PoissonDistance(t(counts(dds.r)))
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf("my.pdf") 
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( vst$tissue, vst$sex, sep="-" )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows=poisd$dd,
         clustering_distance_cols=poisd$dd,
         col=colors)
dev.off()



theme_set(
  theme_bw(base_size = 16)
)

# pca of all feats
pdf("myPCA.pdf") 
pcaData <- plotPCA(vst.r, intgroup=c("sex","tissue","species"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar")) 

p <-ggplot(pcaData, aes(x = PC1, y = PC2, color = sex, shape = tissue, fill = factor(ifelse(species == "cornix", sex, "white")))) +
	scale_color_manual(name = "sex", values=c(wes_palette("Zissou1", n = 5)[5] ,wes_palette("Zissou1", n = 5)[1])) +
	scale_shape_manual(name = "tissue", values = c(21,22,23,24)) +
    scale_fill_manual(name = "population", values = c("#0073C2FF" ,"#CD534CFF", "white"))+
	xlim(-100,100) 

p + 
  # add size for sex
  geom_point(aes(size = species),stroke = 1) +
  # defining size with 2 marginally different values
  scale_size_manual(name = "population", values = c(3, 3.01)) +
  # Remove fill legend and replace the fill legend using the newly created size
  guides(fill = "none", size = guide_legend(override.aes = list(shape = c(1, 16))))+
  theme(panel.border = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_rect(fill = 'white'),
  axis.line = element_line(colour = "black"),
  legend.background = element_rect(fill="white",
                                  size=0.2, linetype="solid", 
                                  colour ="black")) 
								  
dev.off()


