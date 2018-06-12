# setwd("/home/mtormo/ssh_gpfs_projects_lab_genomica/Externs/180413_ACCC6EANXX_ccorral")
library(DESeq2)

### get counts from htseq-count
htseq.path <- "htscount/"

pfile <- read.csv(file = "pheno.csv",sep = "\t",header = T)

sampleFiles <- grep("csv",list.files(htseq.path),value=TRUE)

sampleTable <- data.frame(sampleName = pfile$iden,
                          fileName = sampleFiles,
                          condition = pfile$sample)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = htseq.path,
                                       design= ~ condition)
ddsHTSeq

### collapse replicates
ddsHTSeq2 <- collapseReplicates(ddsHTSeq,groupby = ddsHTSeq$condition)


library(GenomicAlignments)
library(GenomicFeatures)

gtf <- "gtf/Schizosaccharomyces_pombe.ASM294v2.35.chr.gtf"
txdb <- makeTxDbFromGFF(gtf,format="gtf")

### get counts per gene
knownGenes <- exonsBy(txdb, by="gene")

### get fpkm from counts
# then for each gene, reduce all the exons to a set of non overlapping exons, calculate their lengths (widths) and sum then
exonic.gene.sizes <- lapply(knownGenes,function(x){sum(width(reduce(x)))})
sizes <- data.frame(matrix(unlist(exonic.gene.sizes)),row.names = names(exonic.gene.sizes))
colnames(sizes) <- "length"

countToFpkm <- function(counts, effLen){
    N <- sum(counts)
    print(c(counts,N))
    exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}


### join individual and collapsed counts
#mycounts <- merge(assays(ov)$counts,assays(ov2)$counts,by=0)
mycounts <- merge(assays(ddsHTSeq)$counts,assays(ddsHTSeq2)$counts,by=0)
rownames(mycounts)=mycounts$Row.names; mycounts=mycounts[2:length(mycounts)]
colnames(mycounts) <- paste(colnames(mycounts), "counts", sep = "_")
mycounts <- merge(mycounts,sizes,by="row.names",all.x=T)

head(mycounts)

for (i in 1:(ncol(mycounts))){
    # print(i)
    # print(colnames(mycounts)[i])
    if(colnames(mycounts)[i] != "length" & colnames(mycounts)[i] != "Row.names" ){
        # print(head(mycounts[,i]))
        fpkm_col <- sub("_counts", "_fpkm", colnames(mycounts)[i])
        # print(fpkm_col)
        mycounts[,fpkm_col] <- with(mycounts, countToFpkm(mycounts[,i],length))
        # print(head(mycounts[,fpkm_col]))
    }
}

colnames(mycounts)[1] <- "Gene"

write.csv(mycounts,file="ccorral_counts_and_fpkm.csv", quote=F,row.names = F)


################### PLOT MDS to see if they are together
library(edgeR)

dge <- DGEList(counts = assays(ddsHTSeq)$counts, group = colnames(ddsHTSeq))

png("samples_distribution.png")
plotMDS(dge, labels=colnames(ddsHTSeq), col=as.integer(ddsHTSeq$condition)+1, cex=0.8, las=1, main ="Distribution of samples")
dev.off()

