
header<-c("M1A1p","M1C1p","M1E1p","M1F1p","M1H1p","M1J2","M2B1p","M2D1p","M2G1p","M2I1p")

countdata<-read.table("./coral_cds_non-stranded.counts",
                         sep="\t",head=F,row.names=1,as.is=T)
gene_ids<-row.names(countdata)
countdata<-countdata[,1:10]
row.names(countdata)<-gene_ids
colnames(countdata)<-header
dim(countdata)
countdata[1:5,]

#remove genes that have no counts over all samples
countdata<-countdata[(rowSums(countdata)>0),]
dim(countdata) #0 genes with no counts

#Variance Stabilization
library(DESeq2)

inf <- read.table("samples.txt",sep="\t",head =T,row.names=1) # read the sample info in
inf<-inf[names(countdata),] #make sure that the ids in the info table have the same order as the columns in the count table
inf$position<-factor(inf$position)
inf
dds<-DESeqDataSetFromMatrix(countData=as.matrix(countdata),colData=inf,
                                  design=~colony+tank)
dds<-estimateSizeFactors(dds) #library-size correction
sizeFactors(dds)
dds<-estimateDispersions(dds)
plotDispEsts(dds)
vsd<-varianceStabilizingTransformation(dds) #calculate the variance stabilizing normalization
vsdMat<-assay(vsd) #extract a matrix of normalized counts from the S4 object vsd

library(vsn)
meanSdPlot(vsdMat) #plot variance stabilization by DESeq2

#PCA with different variables explaining data
#ntop plots the specified number of most variable genes
plotPCA(vsd,intgroup ="colony",ntop=500)
plotPCA(vsd,intgroup ="tank",ntop=500)

####################################
#Differential Expressed Transcripts#
####################################

info<-as.data.frame(colData(dds))
head(info)
countMatrix<-counts(dds,normalized =F)
design(dds)

dds<-DESeq(dds,test="Wald", #Wald test: calculate difference in deviances (for glm what residuals are for the lm) and their variance and does a kind of t-test
                 betaPrior = T,
           modelMatrixType="expanded" #to get an indicator for each factor level plus intercept, so that base level of factor is not absorbed in intercept
)

C1_vs_C2<-results(dds,contrast=c("colony","C1","C2"), #specifiy contrast as a named c(covariate,level num,level)
                             independentFiltering=TRUE, #kick out genes, for which the test has no chance to become significant, e.g. due to few reads
                             addMLE=TRUE, #also report the unshrunken fold-change added
                             alpha=0.05 #p-value threshold for the Benjamini Hochberg
)
head(C1_vs_C2)
summary(C1_vs_C2) # 69 DEG

#ordering the genes by fold-change
col_fc <- abs(C1_vs_C2$log2FoldChange)
col_ordered <- C1_vs_C2[order(col_fc, decreasing = T), ] #make a reordered DataFrame
head(col_ordered)

#extracting the list of differentially expressed genes
col_genes <- as.data.frame(col_ordered)
col_genes <- na.omit(col_genes)
col_deg <- col_genes[col_genes$padj <= 0.05, ]
col_deg[1:10,]
write.table(col_deg,file="./coral_col_T1_deg.txt",sep="\t")

T_vs_C<-results(dds_noint,contrast=c("tank","treat","control"), #specifiy contrast as a named c(covariate,level num,level)
                  independentFiltering=TRUE, #kick out genes, for which the test has no chance to become significant, e.g. due to few reads
                  addMLE=TRUE, #also report the unshrunken fold-change added
                  alpha=0.05 #p-value threshold for the Benjamini Hochberg
)
head(T_vs_C)
summary(T_vs_C) # 18 DEG

#ordering the genes by fold-change
cond_fc<-abs(T_vs_C$log2FoldChange)
cond_ordered<-T_vs_C[order(cond_fc,decreasing=T),] #make a reordered DataFrame
head(cond_ordered)

#extracting the list of differentially expressed genes
cond_genes<-as.data.frame(cond_ordered)
cond_genes<-na.omit(cond_genes)
cond_deg<-cond_genes[cond_genes$padj<=0.05,]
cond_deg[1:10,]
write.table(cond_deg,file="./coral_cond_T1_deg.txt",sep="\t")
write.table(cond_genes,file="./coral_cond_T1_genes.txt",sep="\t")

###################################################
# Gene Ontology Enrichment for Biological Process #
###################################################

# TopGo now wants a numeric or factor Vector that signifies the grouping
cond_deg <- (T_vs_C$padj < 0.05)
cond_list <- rep(1, dim(cond_genes)[1])
names(cond_list) <- row.names(cond_genes)
cond_list[cond_deg] <- 0.01

cond_GOs <- read.table("./GO_process.tsv",
                     sep = "\t", header = F)
cond_GOs <- cond_GOs[cond_GOs$V3 != "", 2:3]
cond_GOs <- unique(cond_GOs)

write.table(cond_GOs,
            file = "./Condition_T1_BP_deg1.tab",
            quote = F, sep = "\t", row.names = F, col.names = F)
library(topGO)
cond_GOs <- readMappings(file = "./Condition_T1_BP_deg1.tab")
str(head(cond_GOs))
topDiffGenes <- function(x){return(x < 0.05)}

cond_TG <- new("topGOdata",
             description = "deg",
             ontology = "BP",
             allGenes = cond_list, # named vector with adj. p-values
             geneSel = topDiffGenes, #our function
             annot = annFUN.gene2GO,
             gene2GO = cond_GOs)
numGenes(cond_TG)
numSigGenes(cond_TG)

#Run a Fishers-Exact Test, using the eliminination Method to correct for the dependence structure
cond_GO_res <- runTest(cond_TG, algorithm = "elim", statistic = "fisher")
sum(score(cond_GO_res) <= 0.05)
cond_table <- GenTable(cond_TG, Fisher = cond_GO_res,
                     orderBy = "Fisher", ranksOf = "Fisher", topNodes = 137)
write.table(cond_table,
            "./Condition_T1_GO_results.tab",
            sep = "\t")

#how the significant nodes relate to one another within the GO-network
showSigOfNodes(cond_TG, # topGO data object
               score(cond_GO_res), # significance scores
               firstSigNodes = 5, # display the 5 most significant nodes
               useInfo = 'all' # what type of info should be displayed
)
