##### Symbiont #####
######## T1 ########

M1B1<-read.table("../201602_RSEM/Symbiont_CDS/Control-M1B1_symbiont_cds.genes.results",
                 sep="\t",head=T)
M1D1<-read.table("../201602_RSEM/Symbiont_CDS/Control-M1D1_symbiont_cds.genes.results",
                 sep="\t",head=T)
M2A1pb<-read.table("../201602_RSEM/Symbiont_CDS/Control-M2A1pb_symbiont_cds.genes.results",
                   sep="\t",head=T)
M1G1b<-read.table("../201602_RSEM/Symbiont_CDS/Treat-M1G1b_symbiont_cds.genes.results",
                  sep="\t",head=T)
M2F1<-read.table("../201602_RSEM/Symbiont_CDS/Treat-M2F1_symbiont_cds.genes.results",
                 sep="\t",head=T)
M2J1<-read.table("../201602_RSEM/Symbiont_CDS/Treat-M2J1_symbiont_cds.genes.results",
                 sep="\t",head=T)

header_T1<-c("M1B1","M1D1","M2A1pb","M1G1b","M2F1","M2J1")

FPKM_T1<-cbind(M1B1$FPKM, 
            M1D1$FPKM,
            M2A1pb$FPKM,            
            M1G1b$FPKM,
            M2F1$FPKM,
            M2J1$FPKM, deparse.level=2)
FPKM_T1[1:5,]

colnames(FPKM_T1)<-header_T1
rownames(FPKM_T1)<-M1B1$gene_id

FPKM_T1<-FPKM_T1[(rowSums(FPKM_T1)>0),]
dim(FPKM_T1)

transFPKM_T1<-t(FPKM_T1)
expMin_T1<-apply(transFPKM_T1,1,function(x){min(x[x >0])}) #determine the minimal value for each sample
log_fpkm_T1<-log10(transFPKM_T1+expMin_T1)
log_fpkm_T1[,1:5]
log_fpkm_T1<-data.frame(t(log_fpkm_T1))
log_fpkm_T1[1:5,]

fpkm_long_T1<-stack(log_fpkm_T1) # transform the data.frame to have all FPKMs in one column
head(fpkm_long_T1)

names(fpkm_long_T1)=c("FPKM","Sample") # change colnames
library(ggplot2)
p_T1<-ggplot(fpkm_long_T1,aes(x=Sample,y=FPKM))+geom_boxplot()
p_T1<-p_T1+theme(axis.text.x=element_text(angle=90)) #turn x-axis labels to 90 degrees
p_T1+ylab(expression(log[10](FPKM_T1))) #label the y-axis properly

#Transcript abundances with TPM
TPM_T1<-cbind(M1B1$TPM,
           M1D1$TPM,
           M2A1pb$TPM,
           M1G1b$TPM,
           M2F1$TPM,
           M2J1$TPM, deparse.level=2)
TPM_T1[1:5,]
colnames(TPM_T1)<-header_T1
rownames(TPM_T1)<-M1B1$gene_id
expr_genes_T1<-apply(TPM_T1,1,function(i){sum(i>=1)>0}) 
summary(expr_genes_T1) #how many genes occur on average more than once per cell in at least one sample
#33,739 genes found in this study

countdata_T1<-read.table("Symbiodinium_C_Md.fa.transdecoder.cds.despliced",
                         sep="\t",head=F,row.names=1,as.is = T)
countdata_T1<-countdata_T1[,c(1:3,14:16)]
colnames(countdata_T1)<-header_T1
dim(countdata_T1)
countdata_T1[1:6,]

#remove genes that have no counts over all samples
countdata_T1<-countdata_T1[(rowSums(countdata_T1)>0),]
dim(countdata_T1) #1 genes with no counts

count_dist_T1 <- rowSums(countdata_T1)
summary(count_dist_T1)
boxplot(count_dist_T1)

#Variance Stabilization
library(DESeq2)

inf_T1<-read.table("sample_info.txt",sep="\t",head =T,row.names=1) # read the sample info in
inf_T1<-inf_T1[inf_T1$sampling=="T1",]
inf_T1<-inf_T1[names(countdata_T1),] #make sure that the ids in the info table have the same order as the columns in the count table
inf_T1$position<-factor(inf_T1$position)
inf_T1
dds_T1<-DESeqDataSetFromMatrix(countData=as.matrix(countdata_T1),colData=inf_T1,
                            design=~colony+tank)

dds_T1<-estimateSizeFactors(dds_T1) #library-size correction

sizeFactors(dds_T1)
dds_T1<-estimateDispersions(dds_T1)
save(dds_T1,file="symbiont_cds_dds_T1.RData ")
plotDispEsts(dds_T1)

vsd_T1<-varianceStabilizingTransformation(dds_T1) #calculate the variance stabilizing normalization
vsdMat_T1<-assay(vsd_T1) #extract a matrix of normalized counts from the S4 object vsd
save(vsd_T1,file="symbiont_T2_plot_info.RData ")

#comparing normalizations
library(vsn)

meanSdPlot(as.matrix(FPKM_T1)) #plot FPKM normalization
meanSdPlot(vsdMat_T1) #plot variance stabilization by DESeq2

#PCA with different variables explaining data
#ntop plots the specified number of most variable genes
plotPCA(vsd_T1,intgroup ="colony",ntop=500)
plotPCA(vsd_T1,intgroup ="tank",ntop=500)

####################################
#Differential Expressed Transcripts#
####################################

info_T1<-as.data.frame(colData(dds_T1))
head(info_T1)
countMatrix_T1<-counts(dds_T1,normalized =F)
design(dds_T1)

dds_T1<-DESeq(dds_T1,test="Wald", #Wald test: calculate difference in deviances (for glm what residuals are for the lm) and their variance and does a kind of t-test
           modelMatrixType="expanded" #to get an indicator for each factor level plus intercept, so that base level of factor is not absorbed in intercept
)
as.data.frame(mcols(mcols(dds_T1),use.names=T)) #check calculations done

treat_vs_control_T1<-results(dds_T1,contrast=c("tank","treat","control"), #specifiy contrast as a named c(covariate,level num,level)
                             independentFiltering=TRUE, #kick out genes, for which the test has no chance to become significant, e.g. due to few reads
                             addMLE=TRUE, #also report the unshrunken fold-change added
                             alpha=0.05 #p-value threshold for the Benjamini Hochberg
)
head(treat_vs_control_T1)
summary(treat_vs_control_T1) # 0 DEG

C2_vs_C1_T1<-results(dds_T1,contrast=c("colony","C2","C1"), #specifiy contrast as a named c(covariate,level num,level)
                             independentFiltering=TRUE, #kick out genes, for which the test has no chance to become significant, e.g. due to few reads
                             addMLE=TRUE, #also report the unshrunken fold-change added
                             alpha=0.05 #p-value threshold for the Benjamini Hochberg
)
head(C2_vs_C1_T1)
summary(C2_vs_C1_T1) # 0 DEG

######## T0 ########

M1A1p<-read.table("../201602_RSEM/Symbiont_CDS/M1A1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M1C1p<-read.table("../201602_RSEM/Symbiont_CDS/M1C1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M1E1p<-read.table("../201602_RSEM/Symbiont_CDS/M1E1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M1F1p<-read.table("../201602_RSEM/Symbiont_CDS/M1F1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M1H1p<-read.table("../201602_RSEM/Symbiont_CDS/M1H1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M1J2<-read.table("../201602_RSEM/Symbiont_CDS/M1J2_symbiont_cds_non-stranded.genes.results",
                 sep="\t",head=T)
M2B1p<-read.table("../201602_RSEM/Symbiont_CDS/M2B1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M2D1p<-read.table("../201602_RSEM/Symbiont_CDS/M2D1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M2G1p<-read.table("../201602_RSEM/Symbiont_CDS/M2G1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)
M2I1p<-read.table("../201602_RSEM/Symbiont_CDS/M2I1p_symbiont_cds_non-stranded.genes.results",
                  sep="\t",head=T)

header_T0<-c("M1A1p","M1C1p","M1E1p","M1F1p","M1H1p","M1J2","M2B1p","M2D1p","M2G1p","M2I1p")

FPKM_T0<-cbind(M1A1p$FPKM,
            M1C1p$FPKM,
            M1E1p$FPKM,
            M1F1p$FPKM,
            M1H1p$FPKM,
            M1J2$FPKM,
            M2B1p$FPKM,
            M2D1p$FPKM,
            M2G1p$FPKM,
            M2I1p$FPKM, deparse.level=2)
FPKM_T0[1:5,]

colnames(FPKM_T0)<-header_T0
rownames(FPKM_T0)<-M1B1$gene_id

FPKM_T0<-FPKM_T0[(rowSums(FPKM_T0)>0),]
dim(FPKM_T0)

transFPKM_T0<-t(FPKM_T0)
expMin_T0<-apply(transFPKM_T0,1,function(x){min(x[x >0])}) #determine the minimal value for each sample
log_fpkm_T0<-log10(transFPKM_T0+expMin_T0)
log_fpkm_T0[,1:5]
log_fpkm_T0<-data.frame(t(log_fpkm_T0))
log_fpkm_T0[1:5,]

fpkm_long_T0<-stack(log_fpkm_T0) # transform the data.frame to have all FPKMs in one column
head(fpkm_long_T0)

names(fpkm_long_T0)=c("FPKM","Sample") # change colnames
library(ggplot2)
p_T0<-ggplot(fpkm_long_T0,aes(x=Sample,y=FPKM))+geom_boxplot()
p_T0<-p_T0+theme(axis.text.x=element_text(angle=90)) #turn x-axis labels to 90 degrees
p_T0+ylab(expression(log[10](FPKM_T0))) #label the y-axis properly

#Transcript abundances with TPM
TPM_T0<-cbind(M1A1p$TPM,
           M1C1p$TPM,
           M1E1p$TPM,
           M1F1p$TPM,
           M1H1p$TPM,
           M1J2$TPM,
           M2B1p$TPM,
           M2D1p$TPM,
           M2G1p$TPM,
           M2I1p$TPM, deparse.level=2)
TPM_T0[1:5,]
colnames(TPM_T0)<-header
rownames(TPM_T0)<-M1B1$gene_id
expr_genes_T0<-apply(TPM_T0,1,function(i){sum(i>=1)>0}) 
summary(expr_genes_T0) #how many genes occur on average more than once per cell in at least one sample
# 33,833 genes found in this study

countdata_T0<-read.table("Symbiodinium_C_Md.fa.transdecoder.cds.despliced",
                      sep="\t",head=F,row.names=1,as.is = T)
countdata_T0<-countdata_T0[,4:13]
colnames(countdata_T0)<-header_T0
dim(countdata_T0)
countdata[1:6,]

#remove genes that have no counts over all samples
countdata_T0<-countdata_T0[(rowSums(countdata_T0)>0),]
dim(countdata_T0) #1 genes with no counts

count_dist_T0 <- rowSums(countdata_T0)
summary(count_dist_T0)
boxplot(count_dist_T0)

#Variance Stabilization
inf_T0<-read.table("sample_info.txt",sep="\t",head =T,row.names=1) # read the sample info in
inf_T0<-inf_T0[inf_T0$sampling=="T0",]
inf_T0<-inf_T0[names(countdata_T0),] #make sure that the ids in the info table have the same order as the columns in the count table
inf_T0$position<-factor(inf_T0$position)
inf_T0
dds_T0<-DESeqDataSetFromMatrix(countData=as.matrix(countdata_T0),colData=inf_T0,
                            design=~colony+portion+tank)

dds_T0<-estimateSizeFactors(dds_T0) #library-size correction

sizeFactors(dds_T0)
dds_T0<-estimateDispersions(dds_T0)
save(dds_T0,file="symbiont_cds_dds_T0.RData")
plotDispEsts(dds_T0)

vsd_T0<-varianceStabilizingTransformation(dds_T0) #calculate the variance stabilizing normalization
vsdMat_T0<-assay(vsd_T0) #extract a matrix of normalized counts from the S4 object vsd
save(vsd_T0,file="symbiont_T1_plot_info.RData")

#comparing normalizations
meanSdPlot(as.matrix(FPKM_T0)) #plot FPKM normalization
meanSdPlot(vsdMat_T0) #plot variance stabilization by DESeq2

#PCA with different variables explaining data
#ntop plots the specified number of most variable genes
plotPCA(vsd_T0,intgroup ="colony",ntop=500)
plotPCA(vsd_T0,intgroup ="portion",ntop=500)
plotPCA(vsd_T0,intgroup ="tank",ntop=500)

####################################
#Differential Expressed Transcripts#
####################################

info_T0<-as.data.frame(colData(dds_T0))
head(info_T0)
countMatrix_T0<-counts(dds_T0,normalized =F)
design(dds_T0)

dds_T0<-DESeq(dds_T0,test="Wald", #Wald test: calculate difference in deviances (for glm what residuals are for the lm) and their variance and does a kind of t-test
           modelMatrixType="expanded" #to get an indicator for each factor level plus intercept, so that base level of factor is not absorbed in intercept
)
as.data.frame(mcols(mcols(dds_T0),use.names=T)) #check calculations done

treat_vs_control_T0<-results(dds_T0,contrast=c("tank","treat","control"), #specifiy contrast as a named c(covariate,level num,level)
                          independentFiltering=TRUE, #kick out genes, for which the test has no chance to become significant, e.g. due to few reads
                          addMLE=TRUE, #also report the unshrunken fold-change added
                          alpha=0.05 #p-value threshold for the Benjamini Hochberg
)
head(treat_vs_control_T0)
summary(treat_vs_control_T0) # 0 DEG

To_C2_vs_C1<-results(dds,contrast=c("colony","C2","C1"), #specifiy contrast as a named c(covariate,level num,level)
                             independentFiltering=TRUE, #kick out genes, for which the test has no chance to become significant, e.g. due to few reads
                             alpha=0.05 #p-value threshold for the Benjamini Hochberg
)
head(To_C2_vs_C1)
summary(To_C2_vs_C1) # 0 DEG

###############################
######## INTERACTION T0 #######
###############################

dds_T0_int<-DESeqDataSetFromMatrix(countData=as.matrix(countdata_T0),colData=inf_T0,
                            design=~colony+tank+colony:tank)

dds_T0_int<-estimateSizeFactors(dds_T0_int) #library-size correction
dds_T0_int<-estimateDispersions(dds_T0_int)
plotDispEsts(dds_T0_int)
vsd_T0_int<-varianceStabilizingTransformation(dds_T0_int) #calculate the variance stabilizing normalization
vsdMat_T0_int<-assay(vsd_T0_int)
meanSdPlot(vsdMat_T0_int) #plot variance stabilization by DESeq2

#PCA with different variables explaining data
#ntop plots the specified number of most variable genes
plotPCA(vsd_T0_int,intgroup ="colony",ntop=500)
plotPCA(vsd_T0_int,intgroup ="tank",ntop=500)

info_T0_int<-as.data.frame(colData(dds_T0_int))
countMatrix_T0_int<-counts(dds_T0_int,normalized =F)
design(dds_T0_int)

dds_T0_int<-DESeq(dds_T0_int,test="LRT", reduced =~ colony + tank)
res <- results(dds_T0_int,independentFiltering = T,alpha = 0.05)
summary(res)
resultsNames(dds_T0_int)

### Condition effect for Colony 1 ###
T0_Treat_C1<-results(dds_T0_int,contrast=c("tank","treat","control"),independentFiltering = T,alpha = 0.05)
summary(T0_Treat_C1) # 0 DEG

### Condition effect for Colony 2 ###
T0_Treat_C2<-results(dds_T0_int,list(c("tank_treat_vs_control","colonyC2.tanktreat")),independentFiltering = T,alpha = 0.05)
summary(T0_Treat_C2) # 0 DEG

### Is the condition effect different across colonies? ###
T0_Treat_C2_vs_C1<-results(dds_T0_int,name="colonyC2.tanktreat",independentFiltering = T,alpha = 0.05)
summary(T0_Treat_C2_vs_C1) # 0 DEG

T0_T_vs_C_int <- results(dds_T0_int, name="tank_treat_vs_control",independentFiltering = T,alpha = 0.05)
summary(T0_T_vs_C_int) # 0 DEG
