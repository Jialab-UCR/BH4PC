

############################################################################

source('GS_functions.R')
source('HAT_functions.R')
source('Commercial_Panels.R')

library(ggplot2)

############## DATA PROCESSING

### load phenotype data
phe <- read.csv('Data/TCGA-PRAD.phe.clinical.final.csv', header=T, row.names=1)
head(phe)

phe5 <- phe[! is.na(phe$pfr_10yr), ]
dim(phe5)


### load RNAseq, miRNAseq, and Methylatiion data
load('Data/methy450.noNA.TCGA-PRAD.rda')
dim(methy450)

load('Data/rnaExpr495F.rda')
dim(rnaExpr495F)

load('Data/mirExpr.rda')
dim(mirExpr)


samples <- Reduce(intersect, list(rownames(phe5), colnames(mirExpr), colnames(methy450), colnames(rnaExpr495F)))
samples


###

phe5 <- phe5[samples,]
gene <- rnaExpr495F[,samples]
#gene <- mirExpr[,samples]
gene <- as.matrix(t(gene))
gene[1:5,1:5]
dim(gene)
dim(phe5)


gene <- scale(gene)


### pfr_10yr
pheno <- as.matrix(phe5$pfr_10yr, drop=FALSE)
y <- as.numeric(pheno)


corr <- abs(apply(gene, 2, function(v) cor.test(v,y)$estimate))
corrVal <- apply(gene, 2, function(v) cor.test(v,y)$estimate)
o <- order(corr, decreasing=T)
corrDa <- data.frame(corrVAL=corrVal[o],corrABS=corr[o], rank=1:length(o))

write.table(corrDa, file='Results/gene_phe_correlation.pfr10yr.txt', sep='\t', quote=F)
write.table(corrDa[oncotype,], file='Results/gene_phe_correlation.pfr10yr.oncotype.txt', sep='\t', quote=F)
write.table(corrDa[decipher,], file='Results/gene_phe_correlation.pfr10yr.decipher.txt', sep='\t', quote=F)
write.table(corrDa[prolaris,], file='Results/gene_phe_correlation.pfr10yr.prolaris.txt', sep='\t', quote=F)




####################### MODEL COMPARISON

panels <- list(oncotype, decipher, prolaris)

h <- c()
for (i in 1:3) {
    selected <- panels[[i]]
    
    geno <- gene[,selected]
    kk <- kinship(gen=geno)
    
    kk <- kk[[1]]
    kk <- kk[,-c(1,2)]
    kk <- as.matrix(kk)
    
    result1 <- gsm.1d.simple(mydata=y, mykin=kk)
    hat <- result1$predic.HAT
    h <- c(h,hat)
    
}



nGene <- length(o)

for (i in seq(2,nGene,1)) {
    print ('====================================')
    print (i)
    selected <- o[1:i]
    
    geno<-gene[,selected]
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    result1 <- gsm.1d.simple(mydata=y, mykin=kk)
    hat <- result1$predic.HAT
    h <- c(h,hat)
    
}


#p <- data.frame(x=c('oncotype','decipher','prolaris', seq(2,nGene,1)), y=h)
#write.table(p, file='Results/hat.pfr10yr.step1.txt', sep='\t', quote=F, row.names = FALSE)




############################ PLOT

p <- data.frame(x=seq(2,nGene,1), y=h[-c(1:3)])
points <- p[match(c(12,18,31),p$x),]
panel <- factor(c('Oncotype', 'Decipher', 'Prolaris'),
                levels=c('Oncotype', 'Decipher', 'Prolaris'))


pdf('Results/hat.pfr_10yr.step1.pdf', paper='us')
ggplot() + geom_line(data=p, aes(x=x, y=y), size=1) + 
    geom_hline(aes(yintercept = h[1:3], color=panel),linetype=2, size=1) +
    xlab('Number of genes') + ylab('Predictability') + ylim(0,0.5) +
    geom_point(data=points, aes(x=x, y=y, color=panel), size=5, shape=20) +
    theme_bw()+theme(axis.line = element_line(colour = "black"),
                     #panel.grid.minor = element_blank(),
                     axis.text = element_text(size=14),
                     axis.title = element_text(size=16),
                     legend.text = element_text(size=12),
                     legend.title = element_text(size=0),
                     panel.background = element_blank(),
                     panel.border = element_rect(color='black')) 

dev.off()

p[p$y==max(p$y),]
