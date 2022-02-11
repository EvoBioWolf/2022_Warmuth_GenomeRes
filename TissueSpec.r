### calculate the tissue specificity index tau ###

# load libaries
library(DESeq2)
 
# set up function that calculates tau
tau<-function(x){
  if(any(is.na(x))) stop('NA\'s need to be 0.')
  if(any(x<0)) stop('Negative input values not permitted. Maybe data is log transformed?')
  t<-sum(1-x/max(x))/(length(x)-1)
}

# get rlog normalized count data from DESeq object (dds)
rld <- assay(rlog(dds, blind=FALSE))
rld.assay<-assay(rld)
rld.assay<-data.frame(rld.assay)

#counts per tissue across all samples
rldLiver<-apply(rld.assay[,grep("liver",names(rld.assay))],1,mean)
rldSpleen<-apply(rld.assay[,grep("spleen",names(rld.assay))],1,mean)
rldGonads<-apply(rld.assay[,grep("gonads",names(rld.assay))],1,mean)
rldByTissue<-cbind(rldLiver,rldSpleen,rldGonads)

# counts per tissue for females
rldLiver_f<-apply(rld.assay[,grep("f_liver",names(rld.assay))],1,mean)
rldSpleen_f<-apply(rld.assay[,grep("f_spleen",names(rld.assay))],1,mean)
rldGonads_f<-apply(rld.assay[,grep("f_gonads",names(rld.assay))],1,mean)
rldByTissue_f<-cbind(rldLiver_f,rldSpleen_f,rldGonads_f)

#set neg log values to 0
rldByTissue_f[rldByTissue_f<0] <- 0

#calculate tau
TissuesTaus_f <- apply(rldByTissue_f, 1, tau)
TissuesTaus_f<-as.data.frame(TissuesTaus_f)
rownames<-rownames(TissuesTaus_f)
TissuesTaus_f<-cbind(rownames,TissuesTaus_f)
colnames(TissuesTaus_f)<-c("gene_ID","Tau")

# counts per tissue for males
rldLiver_m<-apply(rld.assay[,grep("m_liver",names(rld.assay))],1,mean)
rldSpleen_m<-apply(rld.assay[,grep("m_spleen",names(rld.assay))],1,mean)
rldGonads_m<-apply(rld.assay[,grep("m_gonads",names(rld.assay))],1,mean)
rldByTissue_m<-cbind(rldLiver_m,rldSpleen_m,rldGonads_m)

#set neg log values to 0
rldByTissue_mf[rldByTissue_m<0] <- 0

#calculate tau
TissuesTaus_m <- apply(rldByTissue_m, 1, tau)
TissuesTaus_m<-as.data.frame(TissuesTaus_m)
rownames<-rownames(TissuesTaus_m)
TissuesTaus_m<-cbind(rownames,TissuesTaus_m)
colnames(TissuesTaus_m)<-c("gene_ID","Tau")

