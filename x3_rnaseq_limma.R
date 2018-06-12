### A PARTIR DEL PROTOCOL DEL MASTER UPF: http://functionalgenomics.upf.edu/courses/IEO/DEanalysis2/index.html#21
### i HTSCOUNTS

# setwd("/home/mtormo/ssh_gpfs_projects_lab_genomica/Externs/180413_ACCC6EANXX_ccorral")

out <- "test_DEG/"




###################################################################
### 1. PREPARE THE DATA
### load data:
# library(SummarizedExperiment)
# se <- readRDS("countsSOunionFrg.STAR.GRCm38.84.gtf.ignore-true.rds")

library(DESeq2)
### get counts from htseq-count
htseq.path <- "htscount/"
### phenotype data
pfile <- read.csv(file = "pheno.csv",sep = "\t",header = T)

sampleFiles <- grep("csv",list.files(htseq.path),value=TRUE)

sampleTable <- data.frame(sampleName = pfile$iden,
                          fileName = sampleFiles,
                          condition = pfile$sample)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = htseq.path,
                                       design= ~ condition)
ddsHTSeq


### Add gene metadata: get gene information from GTF
require(rtracklayer)
gtf_tags <- readGFF("gtf/Schizosaccharomyces_pombe.ASM294v2.35.chr.gtf", tags = c("gene_id", "gene_name", "gene_biotype"))
genemdata <- data.frame(gtf_tags$gene_id,gtf_tags$gene_name,gtf_tags$gene_biotype)
### get uniques
genemdata <- unique(genemdata)
rownames(genemdata) <- genemdata[,1]
### assign rownames (gene_id)
genemdata[,1] <- NULL
genemdata <- genemdata[ order(rownames(genemdata)) , ]
head(genemdata)

nrow(ddsHTSeq)
nrow(genemdata)

mcols(rowRanges(ddsHTSeq)) <- genemdata

### create a dgelist object, including metadata
library(edgeR)
dge <- DGEList(counts = assays(ddsHTSeq)$counts, group = ddsHTSeq$condition, genes = as.data.frame(mcols(ddsHTSeq)))
dim(dge)
names(dge)


### 2. QUALITY ASSESSMENT

###################################################################
### Sequencing across lanes (ONLY WITH NEXTSEQ)

# png(sprintf("%s1_lanes_distribution.png",out))
# plotMDS(dge, labels=colnames(ddsHTSeq), col=as.integer(ddsHTSeq$condition)+1, cex=0.8, las=1, main ="Distribution of samples across lanes")
# 
# # plotMDS(dge, labels=as.integer(as.factor(colnames(ddsHTSeq))), col=as.integer(ddsHTSeq$condition)+1, cex=0.8, las=1, main ="Distribution of samples across lanes")
# # plotMDS(dge, labels=as.integer(se$PatientId), col=as.integer(se$Lane)+1, cex=0.8, las=1)
# # legend("topleft", levels(se$Lane), fill=2:(nlevels(se$Lane)+1))
# # legend("topleft", colnames(ddsHTSeq), fill=2:(nlevels(ddsHTSeq$condition)+1))
# dev.off()
# 
# ### they cluster together, collapse replicates then
# dds <- collapseReplicates(ddsHTSeq,groupby = ddsHTSeq$Lane)


#### with summarized experiment
# ### they cluster together by ID, let's sum counts per ID
# ### create unique names (without lane id)
# f <- factor(sub("\\_L00[0-9]", "", colnames(se)))
# ### create a matrix with indexes
# idx <- matrix(rep(f, times=nrow(se)), ncol=ncol(se), byrow=TRUE)
# ### get index of first's occurence sample and create colnames
# mt <- which(!duplicated(idx[1, ]))
# names(mt) <- idx[1, !duplicated(idx[1, ])]
# ### count by sample: create a matrix with all values of each sample, sum counts, and put the name
# cntby <- split(assays(se)$counts, idx)
# cntby <- lapply(cntby, matrix, nrow=nrow(se))
# cntby <- sapply(cntby, rowSums)
# rownames(cntby) <- rownames(se)
# cntby <- cntby[, names(mt)]
# ### get the colData without Lane
# cdata <- colData(se)[mt, -match("Lane", colnames(colData(se)))]
# rownames(cdata) <- names(mt)
# ### get the new summarizedExperiment
# se2 <- SummarizedExperiment(assays=list(counts=cntby),
#                             rowRanges=rowRanges(se),
#                             colData=cdata)
# se2
# se <- se2
# rm(se2,cdata,cntby,idx,pdata,pfile)
# 
# ### create a dgelist object
# dge <- DGEList(counts = assays(se)$counts)
# dim(dge)
# ### Calculate log2-CPM values and store them into the SummarizedExperiment object.
# assays(se)$logCPM <- cpm(dge, log = TRUE, prior.count=0.5)

# ### create a dgelist object
# dge <- DGEList(counts = assays(dds)$counts, genes = as.data.frame(mcols(dds)))
# dim(dge)

### without lane concatenation, impute dds object
dds <- ddsHTSeq

### Calculate log2-CPM values and store them into the DDS object.
assays(dds)$logCPM <- cpm(dge, log = TRUE, prior.count=0.5)

###################################################################
### Library sizes (total number of sequence read counts per sample)

count_per_sample <- as.data.frame(colSums(assays(dds)$counts))
colnames(count_per_sample) <- "counts"
png(sprintf("%s2_reads_per_sample.png",out))
barplot(sort(count_per_sample$counts/1000000), xlab = "Samples", ylab = "M Reads",main = "Reads per sample", names.arg = rownames(count_per_sample),cex.names = 0.6,las=2)
dev.off()

###################################################################
### Distribution of expression levels among samples

library(lattice)
library(reshape2)

cpm_log <- melt(assays(dds)$logCPM)
png(sprintf("%s3_expression_samples.png",out))
densityplot(~ value, groups = Var2, data = cpm_log, xlab="log2-CPM",main="Expression among samples",plot.points=FALSE,auto.key=TRUE)
dev.off()


###################################################################
### Filtering of lowly-expressed genes

avgexp <- rowMeans(assays(dds)$logCPM)
png(sprintf("%s4_expression_genes.png",out))
hist(avgexp, xlab="log2 CPM", main="Distribution of average expression level per gene", las=1)
abline(v=1, col="red", lwd=2)
dev.off()

### Filter out lowly-expressed genes.
mask <- avgexp > 1
# dim(dds)
# [1] 47729    20
dds.filt <- dds[mask, ]
# dim(dds.filt)
# [1] 13342    20

# dim(dge)
# [1] 47729    20
dge.filt <- dge[mask, ]
# dim(dge.filt)
# [1] 13342    20


###################################################################
### 3. NORMALIZATION

### calculate normalization factors
### The Trimmed Mean of M-values (TMM) method addresses the issue of the different RNA composition of the samples by estimating a scaling factor for each library
dge.filt <- calcNormFactors(dge.filt)
### and replace cpm by normalized cpm
assays(dds.filt)$logCPM <- cpm(dge.filt, normalized.lib.sizes=TRUE,log=TRUE, prior.count=0.5)


###################################################################
### MA-plots: one condition VS another

png(sprintf("%s5_MA-plots.png",out))
par(mfrow = c(1, 2), oma=c(0,0,1,0))

### first: between samples
dge.filt$samples$group <- dds.filt$condition
table(dge.filt$samples$group)
plotSmear(dge.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
abline(h = 0, col = "blue", lwd = 2)

# ### second: between treatments
# dge.filt$samples$group <- dds.filt$treatment
# table(dge.filt$samples$group)
# plotSmear(dge.filt, lowess = TRUE, las = 1, cex.lab = 1.5, cex.axis = 1.2)
# abline(h = 0, col = "blue", lwd = 2)

title("MA plots", outer = TRUE)
dev.off()


###################################################################
### Batch identification

### EDITAR!!!

png(sprintf("%s6_batch_identification1.png",out))
par(mfrow = c(1, 2), oma=c(0,0,1,0))
plotMDS(dge.filt, labels=colnames(dds.filt), col=as.integer(dds.filt$condition)+1, cex=0.8, las=1)
legend("topleft", levels(dds.filt$treatment), fill=2:(nlevels(dds.filt$treatment)+1))
# plotMDS(dge.filt, labels=dds.filt$condition, col=as.integer(dds.filt$sample)+1, cex=0.8, las=1)
# legend("topleft", levels(se.filt$sample), fill=2:(nlevels(se.filt$sample)+1))
title("Batch identification", outer = TRUE)
dev.off()
png(sprintf("%s6_batch_identification2.png",out))
par(mfrow = c(1, 1), oma=c(0,0,1,0))
plotMDS(dge.filt, labels=colnames(dds.filt), col=as.integer(dds.filt$condition)+1, cex=0.8, las=1, main="Batch identification by sample & treatment")
legend("topleft", levels(dds.filt$samp_treat), fill=2:(nlevels(dds.filt$samp_treat)+1),cex=0.6)
dev.off()

### cluster dendrogram
tcounts <- t(assays(dds.filt)$logCPM)
clusters <- hclust(dist(tcounts))
plot(clusters)


###################################################################
### 4. DIFFERENTIAL EXPRESSION


###################################################################
### FACTORIAL DESIGNS
### Factorial designs, also known as factorial experiments, are experimental designs where the outcome of interest consists of two or more factors.
### The limma Users's Guide, Section 9.5 and elsewhere, provides guidance and specific examples on how to analyze a factorial design. 
### There are several ways to approach this kind of analysis. One that facilitates the interpretation is building a single factor out of the combination
### of the factors under study and set the intercept term to zero.
### Fit a linear regression model using this new factor variable without intercept term.
mod <- model.matrix(~0 + dds.filt$condition, colData(dds.filt))
colnames(mod) <- levels(dds.filt$condition)
mod0 <- model.matrix(~1, colData(dds.filt))
### We can adjust for unkown covariates using surrogate variable analysis (SVA). 
### First, estimate surrogate variables (SVs) from the log-CPM values calculated by voom:
library(sva)
sv <- sva(assays(dds.filt)$logCPM, mod = mod, mod0 = mod0)
### Second, add these (SVs) to the design matrix:
mod <- cbind(mod, sv$sv)
colnames(mod) <- c(colnames(mod)[1:3], paste0("SV", 1:sv$n))
### Third, calculate again mean-variance weights:
v <- voom(dge.filt, mod, plot = FALSE)

### adjust for correlations with technical repilcats: NO!
# corfit <- duplicateCorrelation(v, mod, block = se.filt$Iden)
# v <- voom(dge, mod, block = se$celline, correlation = corfit$consensus, plot = FALSE)
# fit6 <- lmFit(v, mod, block = se$cellline, correlation = corfit$consensus)

### Fourth, fit again the linear models for each gene with the updated design matrix:
fit6 <- lmFit(v, mod)


### EDITAR!!
### Build the contrast matrix that specifies the contrasts of interest between a set of coefficients estimated from a linear regression model. 
#cont.matrix <- makeContrasts(wt2=wt_rim-wt_v, ko2=fx_rim-fx_v, treat=fx_rim-wt_rim, untreated=fx_v-wt_v, levels=mod)
cont.matrix <- makeContrasts(can=CAN, hs=HS, mm=MM, levels=mod)
### Estimate coefficients for a given set of contrasts.
fit6 <- contrasts.fit(fit6, cont.matrix)
### Calculate moderated t-statistics.
fit6 <- eBayes(fit6)
head(fit6$t)
### Fetch table of results for each coefficient.
# wt2 = topTable(fit6, coef = "wt2", n=Inf)
# treated = topTable(fit6, coef = "treat", n=Inf)
# untreated = topTable(fit6, coef = "untreated", n=Inf)
cantt = topTable(fit6, coef = "can", n=Inf)

# nrow(wt2[wt2$adj.P.Val<0.1,])
# nrow(ko2[ko2$adj.P.Val<0.1,])
# nrow(treated[treated$adj.P.Val<0.1,])
# nrow(untreated[untreated$adj.P.Val<0.1,])
nrow(cantt[cantt$adj.P.Val<0.05,])

### Explore graphically the overlap of DE genes between contrasts of interest.
res <- decideTests(fit6, p.value = 0.1)

png(sprintf("%s7_pvalues_distribution.png",out))
par(mfrow = c(2, 2), oma=c(0,0,1,0))
# hist(wt2$P.Value, xlab = "Raw P-values", main = "WT", las = 1)
# hist(ko2$P.Value, xlab = "Raw P-values", main = "KO", las = 1)
# hist(treated$P.Value, xlab = "Raw P-values", main = "TREATED", las = 1)
# hist(untreated$P.Value, xlab = "Raw P-values", main = "UNTREATED", las = 1)
hist(cantt$P.Value, xlab = "Raw P-values", main = "CAN", las = 1)
title("Raw P-values distribution",outer=T)
dev.off()

summary(res)
head(fit6$t["ENSMUSG00000000838",])

png(sprintf("%s8_venn_diagram.png",out))
vennDiagram(res)
dev.off()




### volcanoplot
# volcanoplot(fit6, coef = 2, highlight = 7, names = rownames(fit6$t), main = "DE genes",las = 1)
# ### 
# top <- order(fit6$lods[, 2], decreasing = TRUE)[1:7]
# limma::plotMA(fit6, coef = 2, status = rownames(fit6$lods), legend = FALSE,main = "Model 5", hl.pch = 46, hl.cex = 4, bg.pch = 46, bg.cex = 3, las = 1)
# text(fit6$Amean[top], fit6$coef[top, 2],  cex = 0.5, pos = 4)

save.image()
load(".RData")

