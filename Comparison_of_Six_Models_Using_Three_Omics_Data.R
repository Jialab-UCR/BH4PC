
######################################################################

######################################################################

source('GS_functions.R')
source('HAT_functions.R')
source('Commercial_Panels.R')

library(ggplot2)
library(foreach)
library(glmnet)
library(kernlab)
library(BGLR)
library(randomForest)
library(pls)


trait <- 'pfr10yr'
phe <- read.csv('Data/TCGA-PRAD.phe.clinical.final.csv', header=T, row.names=1)
head(phe)

phe5 <- phe[! is.na(phe$pfr_10yr), ]
dim(phe5)

load('Data/methy450.noNA.TCGA-PRAD.rda')
dim(methy450)

load('Data/rnaExpr495F.rda')
dim(rnaExpr495F)

load('Data/mirExpr.rda')
dim(mirExpr)

samples <- Reduce(intersect, list(rownames(phe5), colnames(mirExpr), colnames(methy450), colnames(rnaExpr495F)))


#samples <- intersect(rownames(phe), colnames(rnaExpr490F))
#samples
phe5<-phe5[samples,]

### pfr_10yr
pheno <- as.matrix(phe5$pfr_10yr, drop=FALSE)
y <- as.numeric(pheno)

rna <- rnaExpr495F[,samples]
mir <- mirExpr[,samples]
methy <- methy450[,samples]

dim(rna)

load('Data/foldidID378.10.rda')

nfold <- 10


############################## RNAseq #############################

geno <- as.matrix(t(rna))
geno <- scale(geno)
dim(geno)

res<-NULL
ypp1<-NULL
yoo1<-NULL
ypp2<-NULL
yoo2<-NULL
ypp3<-NULL
yoo3<-NULL
ypp4<-NULL
yoo4<-NULL
ypp5<-NULL
yoo5<-NULL
ypp6<-NULL
yoo6<-NULL


for (z in 1:10) {
    cat ('============================================\n')
    cat (rep(c(z,'-'),10))
    foldid <- foldidID[[z]]
    
    
    ############################### BLUP ################################
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    n<-length(pheno)
    x<-matrix(1,n,1)

    blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
    br2<- as.numeric(blup[[1]])
    BLUP <- sqrt(br2)
    
    ypp1<-cbind(ypp1,blup[[2]]$yhat)
    yoo1<-cbind(yoo1,blup[[2]]$yobs)
    
    
    ############################### LASSO ################################
    
    
    fit<-cv.glmnet(x=geno,y=pheno,foldid=foldid)
    theta<-coef(fit,s="lambda.min")
    beta<-theta[1]
    gamma<-theta[-1]
    yhat<-predict(fit,newx=geno,s="lambda.min")
    lambda<-fit$lambda.min
    index<-match(lambda,fit$lambda)
    mse.cv<-fit$cvm[index]
    nonzero<-fit$nzero[index]
    sigma2<-sum((pheno-yhat)^2)/nrow(pheno)
    vp<-sum((pheno-mean(pheno))^2)/nrow(pheno)
    #goodness<-drop(cor(y,yhat)^2)
    pred_mse<-1-mse.cv/vp
    pred_r2<-pearson(geno,pheno,foldid)
    R2<- pred_r2[[1]]
    LASSO <- sqrt(R2)
    
    ypp2<-cbind(ypp2,pred_r2[[2]]$yyhat)
    yoo2<-cbind(yoo2,pred_r2[[2]]$yyobs)
    
    
    ############################### PLS ################################

    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        pls.fit <- plsr(y1~x1,ncomp=5,validation="CV")
        #plot(RMSEP(pls.fit), legendpos = "topright")
        nn <-as.numeric( which.min(tt <-RMSEP(pls.fit)$val[1,,][-1]))
        ## predicting
        yhat<-predict(pls.fit, newdata=x2,ncomp=nn)
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    PLS<-cor(yo,yp)
    
    ypp3<-cbind(ypp3,yp)
    yoo3<-cbind(yoo3,yo)
    
    
    ############################### SSVS ################################
    
    x<-rep(1,length(pheno))
    x<-as.matrix(x)
    
    ETA<-list( list(X=x,model='FIXED'),
               list(X=geno, model='BayesB'))
    fm<-BGLR(y=pheno,ETA=ETA, nIter=1500, burnIn=500)
    uHat<- fm$ETA[[2]]$b
    yHat<-fm$yHat
    r2.fit<-cor(y,yHat)^2

    
    yobs<-NULL
    yhat<-NULL
    id<-NULL
    #nfold<-6
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-pheno[i1,,drop=F]
        z1<-geno[i1,,drop=F]
        ETA<-list( list(X=x1,model='FIXED'),
                   list(X=z1, model='BayesB'))
        fm<-BGLR(y=y1,ETA=ETA, nIter=1500, burnIn=500)
        b<-fm$ETA[[1]]$b
        u <- fm$ETA[[2]]$b
        x2<-x[i2,,drop=F]
        y2<-pheno[i2,,drop=F]
        z2<-geno[i2,,drop=F]
        y3<-x2*b+z2%*%u+fm$mu
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        #fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    
    BayesB<-cor(yobs,yhat)
    
    ypp4<-cbind(ypp4,yhat)
    yoo4<-cbind(yoo4,yobs)
    
    
    ############################### SVM-RBF ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "rbfdot",kpar = "automatic",cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMRBF<-cor(yo,yp)
    
    ypp5<-cbind(ypp5,yp)
    yoo5<-cbind(yoo5,yo)
    
    
    ############################### SVM-POLY ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "poly",kpar = "automatic", cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMPOLY<-cor(yo,yp)
    
    ypp6<-cbind(ypp6,yp)
    yoo6<-cbind(yoo6,yo)
    
    result<- data.frame(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
    res<-rbind(res,result)
}

mean<- apply(res,2,mean)
sd<- apply(res,2,sd)
out<-rbind(mean,sd)
res2=rbind(res,out)
res2
write.table(res2,file=paste('Results/six_methods.RNAseq.scale.', trait, '.txt', sep=''),
            sep='\t', quote=F, row.names=F)



############################## miRNA #############################

geno <- as.matrix(t(mir))
geno <- scale(geno)
geno[1:5,1:5]

res<-NULL
ypp1<-NULL
yoo1<-NULL
ypp2<-NULL
yoo2<-NULL
ypp3<-NULL
yoo3<-NULL
ypp4<-NULL
yoo4<-NULL
ypp5<-NULL
yoo5<-NULL
ypp6<-NULL
yoo6<-NULL


for (z in 1:10) {
    cat ('============================================\n')
    cat (rep(c(z,'-'),10))
    foldid <- foldidID[[z]]
    
    
    ############################### BLUP ################################
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    n<-length(pheno)
    x<-matrix(1,n,1)
    
    blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
    br2<- as.numeric(blup[[1]])
    BLUP <- sqrt(br2)
    
    ypp1<-cbind(ypp1,blup[[2]]$yhat)
    yoo1<-cbind(yoo1,blup[[2]]$yobs)
    
    
    ############################### LASSO ################################
    
    
    fit<-cv.glmnet(x=geno,y=pheno,foldid=foldid)
    theta<-coef(fit,s="lambda.min")
    beta<-theta[1]
    gamma<-theta[-1]
    yhat<-predict(fit,newx=geno,s="lambda.min")
    lambda<-fit$lambda.min
    index<-match(lambda,fit$lambda)
    mse.cv<-fit$cvm[index]
    nonzero<-fit$nzero[index]
    sigma2<-sum((pheno-yhat)^2)/nrow(pheno)
    vp<-sum((pheno-mean(pheno))^2)/nrow(pheno)
    #goodness<-drop(cor(y,yhat)^2)
    pred_mse<-1-mse.cv/vp
    pred_r2<-pearson(geno,pheno,foldid)
    R2<- pred_r2[[1]]
    LASSO <- sqrt(R2)
    
    ypp2<-cbind(ypp2,pred_r2[[2]]$yyhat)
    yoo2<-cbind(yoo2,pred_r2[[2]]$yyobs)
    
    
    ############################### PLS ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        pls.fit <- plsr(y1~x1,ncomp=5,validation="CV")
        #plot(RMSEP(pls.fit), legendpos = "topright")
        nn <-as.numeric( which.min(tt <-RMSEP(pls.fit)$val[1,,][-1]))
        ## predicting
        yhat<-predict(pls.fit, newdata=x2,ncomp=nn)
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    PLS<-cor(yo,yp)
    
    ypp3<-cbind(ypp3,yp)
    yoo3<-cbind(yoo3,yo)
    
    
    ############################### SSVS ################################
    
    x<-rep(1,length(pheno))
    x<-as.matrix(x)
    
    ETA<-list( list(X=x,model='FIXED'),
               list(X=geno, model='BayesB'))
    fm<-BGLR(y=pheno,ETA=ETA, nIter=1500, burnIn=500)
    uHat<- fm$ETA[[2]]$b
    yHat<-fm$yHat
    r2.fit<-cor(y,yHat)^2
    
    
    yobs<-NULL
    yhat<-NULL
    id<-NULL
    #nfold<-6
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-pheno[i1,,drop=F]
        z1<-geno[i1,,drop=F]
        ETA<-list( list(X=x1,model='FIXED'),
                   list(X=z1, model='BayesB'))
        fm<-BGLR(y=y1,ETA=ETA, nIter=1500, burnIn=500)
        b<-fm$ETA[[1]]$b
        u <- fm$ETA[[2]]$b
        x2<-x[i2,,drop=F]
        y2<-pheno[i2,,drop=F]
        z2<-geno[i2,,drop=F]
        y3<-x2*b+z2%*%u+fm$mu
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        #fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    
    BayesB<-cor(yobs,yhat)
    
    ypp4<-cbind(ypp4,yhat)
    yoo4<-cbind(yoo4,yobs)
    
    
    ############################### SVM-RBF ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "rbfdot",kpar = "automatic",cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMRBF<-cor(yo,yp)
    
    ypp5<-cbind(ypp5,yp)
    yoo5<-cbind(yoo5,yo)
    
    
    ############################### SVM-POLY ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "poly",kpar = "automatic", cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMPOLY<-cor(yo,yp)
    
    ypp6<-cbind(ypp6,yp)
    yoo6<-cbind(yoo6,yo)
    
    result<- data.frame(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
    res<-rbind(res,result)
}

mean<- apply(res,2,mean)
sd<- apply(res,2,sd)
out<-rbind(mean,sd)
res2=rbind(res,out)
res2
write.table(res2,file=paste('Results/six_methods.miRNAs.scale.', trait, '.txt', sep=''),
            sep='\t', quote=F, row.names=F)




############################## Methylation #############################

geno <- as.matrix(t(methy))
geno <- scale(geno)
dim(geno)

res<-NULL
ypp1<-NULL
yoo1<-NULL
ypp2<-NULL
yoo2<-NULL
ypp3<-NULL
yoo3<-NULL
ypp4<-NULL
yoo4<-NULL
ypp5<-NULL
yoo5<-NULL
ypp6<-NULL
yoo6<-NULL


for (z in 1:10) {
    cat ('============================================\n')
    cat (rep(c(z,'-'),10))
    foldid <- foldidID[[z]]
    
    
    ############################### BLUP ################################
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    n<-length(pheno)
    x<-matrix(1,n,1)
    
    blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
    br2<- as.numeric(blup[[1]])
    BLUP <- sqrt(br2)
    
    ypp1<-cbind(ypp1,blup[[2]]$yhat)
    yoo1<-cbind(yoo1,blup[[2]]$yobs)
    
    
    ############################### LASSO ################################
    
    
    fit<-cv.glmnet(x=geno,y=pheno,foldid=foldid)
    theta<-coef(fit,s="lambda.min")
    beta<-theta[1]
    gamma<-theta[-1]
    yhat<-predict(fit,newx=geno,s="lambda.min")
    lambda<-fit$lambda.min
    index<-match(lambda,fit$lambda)
    mse.cv<-fit$cvm[index]
    nonzero<-fit$nzero[index]
    sigma2<-sum((pheno-yhat)^2)/nrow(pheno)
    vp<-sum((pheno-mean(pheno))^2)/nrow(pheno)
    #goodness<-drop(cor(y,yhat)^2)
    pred_mse<-1-mse.cv/vp
    pred_r2<-pearson(geno,pheno,foldid)
    R2<- pred_r2[[1]]
    LASSO <- sqrt(R2)
    
    ypp2<-cbind(ypp2,pred_r2[[2]]$yyhat)
    yoo2<-cbind(yoo2,pred_r2[[2]]$yyobs)
    
    
    ############################### PLS ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        pls.fit <- plsr(y1~x1,ncomp=5,validation="CV")
        #plot(RMSEP(pls.fit), legendpos = "topright")
        nn <-as.numeric( which.min(tt <-RMSEP(pls.fit)$val[1,,][-1]))
        ## predicting
        yhat<-predict(pls.fit, newdata=x2,ncomp=nn)
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    PLS<-cor(yo,yp)
    
    ypp3<-cbind(ypp3,yp)
    yoo3<-cbind(yoo3,yo)
    
    
    ############################### SSVS ################################
    
    x<-rep(1,length(pheno))
    x<-as.matrix(x)
    
    ETA<-list( list(X=x,model='FIXED'),
               list(X=geno, model='BayesB'))
    fm<-BGLR(y=pheno,ETA=ETA, nIter=1500, burnIn=500)
    uHat<- fm$ETA[[2]]$b
    yHat<-fm$yHat
    r2.fit<-cor(y,yHat)^2
    
    
    yobs<-NULL
    yhat<-NULL
    id<-NULL
    #nfold<-6
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-pheno[i1,,drop=F]
        z1<-geno[i1,,drop=F]
        ETA<-list( list(X=x1,model='FIXED'),
                   list(X=z1, model='BayesB'))
        fm<-BGLR(y=y1,ETA=ETA, nIter=1500, burnIn=500)
        b<-fm$ETA[[1]]$b
        u <- fm$ETA[[2]]$b
        x2<-x[i2,,drop=F]
        y2<-pheno[i2,,drop=F]
        z2<-geno[i2,,drop=F]
        y3<-x2*b+z2%*%u+fm$mu
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        #fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    
    BayesB<-cor(yobs,yhat)
    
    ypp4<-cbind(ypp4,yhat)
    yoo4<-cbind(yoo4,yobs)
    
    
    ############################### SVM-RBF ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "rbfdot",kpar = "automatic",cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMRBF<-cor(yo,yp)
    
    ypp5<-cbind(ypp5,yp)
    yoo5<-cbind(yoo5,yo)
    
    
    ############################### SVM-POLY ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "poly",kpar = "automatic", cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMPOLY<-cor(yo,yp)
    
    ypp6<-cbind(ypp6,yp)
    yoo6<-cbind(yoo6,yo)
    
    result<- data.frame(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
    res<-rbind(res,result)
}

mean<- apply(res,2,mean)
sd<- apply(res,2,sd)
out<-rbind(mean,sd)
res2=rbind(res,out)
res2
write.table(res2,file=paste('Results/six_methods.methylation.scale.', trait, '.txt', sep=''),
            sep='\t', quote=F, row.names=F)





############################## RNAseq + miRNAs #############################
geno <- as.matrix(t(rbind(rna, mir)))
geno <- scale(geno)
dim(geno)

res<-NULL
ypp1<-NULL
yoo1<-NULL
ypp2<-NULL
yoo2<-NULL
ypp3<-NULL
yoo3<-NULL
ypp4<-NULL
yoo4<-NULL
ypp5<-NULL
yoo5<-NULL
ypp6<-NULL
yoo6<-NULL


for (z in 1:10) {
    cat ('============================================\n')
    cat (rep(c(z,'-'),10))
    foldid <- foldidID[[z]]
    
    
    ############################### BLUP ################################
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    n<-length(pheno)
    x<-matrix(1,n,1)
    
    blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
    br2<- as.numeric(blup[[1]])
    BLUP <- sqrt(br2)
    
    ypp1<-cbind(ypp1,blup[[2]]$yhat)
    yoo1<-cbind(yoo1,blup[[2]]$yobs)
    
    
    ############################### LASSO ################################
    
    
    fit<-cv.glmnet(x=geno,y=pheno,foldid=foldid)
    theta<-coef(fit,s="lambda.min")
    beta<-theta[1]
    gamma<-theta[-1]
    yhat<-predict(fit,newx=geno,s="lambda.min")
    lambda<-fit$lambda.min
    index<-match(lambda,fit$lambda)
    mse.cv<-fit$cvm[index]
    nonzero<-fit$nzero[index]
    sigma2<-sum((pheno-yhat)^2)/nrow(pheno)
    vp<-sum((pheno-mean(pheno))^2)/nrow(pheno)
    #goodness<-drop(cor(y,yhat)^2)
    pred_mse<-1-mse.cv/vp
    pred_r2<-pearson(geno,pheno,foldid)
    R2<- pred_r2[[1]]
    LASSO <- sqrt(R2)
    
    ypp2<-cbind(ypp2,pred_r2[[2]]$yyhat)
    yoo2<-cbind(yoo2,pred_r2[[2]]$yyobs)
    
    
    ############################### PLS ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        pls.fit <- plsr(y1~x1,ncomp=5,validation="CV")
        #plot(RMSEP(pls.fit), legendpos = "topright")
        nn <-as.numeric( which.min(tt <-RMSEP(pls.fit)$val[1,,][-1]))
        ## predicting
        yhat<-predict(pls.fit, newdata=x2,ncomp=nn)
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    PLS<-cor(yo,yp)
    
    ypp3<-cbind(ypp3,yp)
    yoo3<-cbind(yoo3,yo)
    
    
    ############################### SSVS ################################
    
    x<-rep(1,length(pheno))
    x<-as.matrix(x)
    
    ETA<-list( list(X=x,model='FIXED'),
               list(X=geno, model='BayesB'))
    fm<-BGLR(y=pheno,ETA=ETA, nIter=1500, burnIn=500)
    uHat<- fm$ETA[[2]]$b
    yHat<-fm$yHat
    r2.fit<-cor(y,yHat)^2
    
    
    yobs<-NULL
    yhat<-NULL
    id<-NULL
    #nfold<-6
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-pheno[i1,,drop=F]
        z1<-geno[i1,,drop=F]
        ETA<-list( list(X=x1,model='FIXED'),
                   list(X=z1, model='BayesB'))
        fm<-BGLR(y=y1,ETA=ETA, nIter=1500, burnIn=500)
        b<-fm$ETA[[1]]$b
        u <- fm$ETA[[2]]$b
        x2<-x[i2,,drop=F]
        y2<-pheno[i2,,drop=F]
        z2<-geno[i2,,drop=F]
        y3<-x2*b+z2%*%u+fm$mu
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        #fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    
    BayesB<-cor(yobs,yhat)
    
    ypp4<-cbind(ypp4,yhat)
    yoo4<-cbind(yoo4,yobs)
    
    
    ############################### SVM-RBF ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "rbfdot",kpar = "automatic",cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMRBF<-cor(yo,yp)
    
    ypp5<-cbind(ypp5,yp)
    yoo5<-cbind(yoo5,yo)
    
    
    ############################### SVM-POLY ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "poly",kpar = "automatic", cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMPOLY<-cor(yo,yp)
    
    ypp6<-cbind(ypp6,yp)
    yoo6<-cbind(yoo6,yo)
    
    result<- data.frame(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
    res<-rbind(res,result)
}

mean<- apply(res,2,mean)
sd<- apply(res,2,sd)
out<-rbind(mean,sd)
res2=rbind(res,out)
res2
write.table(res2,file=paste('Results/six_methods.RNAseq_and_miRNAs.', trait, '.txt', sep=''),
            sep='\t', quote=F, row.names=F)






############################## RNAseq + methylation #############################


geno <- as.matrix(t(rbind(rna, methy)))
geno <- scale(geno)
dim(geno)

res<-NULL
ypp1<-NULL
yoo1<-NULL
ypp2<-NULL
yoo2<-NULL
ypp3<-NULL
yoo3<-NULL
ypp4<-NULL
yoo4<-NULL
ypp5<-NULL
yoo5<-NULL
ypp6<-NULL
yoo6<-NULL


for (z in 1:10) {
    cat ('============================================\n')
    cat (rep(c(z,'-'),10))
    foldid <- foldidID[[z]]
    
    
    ############################### BLUP ################################
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    n<-length(pheno)
    x<-matrix(1,n,1)
    
    blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
    br2<- as.numeric(blup[[1]])
    BLUP <- sqrt(br2)
    
    ypp1<-cbind(ypp1,blup[[2]]$yhat)
    yoo1<-cbind(yoo1,blup[[2]]$yobs)
    
    
    ############################### LASSO ################################
    
    
    fit<-cv.glmnet(x=geno,y=pheno,foldid=foldid)
    theta<-coef(fit,s="lambda.min")
    beta<-theta[1]
    gamma<-theta[-1]
    yhat<-predict(fit,newx=geno,s="lambda.min")
    lambda<-fit$lambda.min
    index<-match(lambda,fit$lambda)
    mse.cv<-fit$cvm[index]
    nonzero<-fit$nzero[index]
    sigma2<-sum((pheno-yhat)^2)/nrow(pheno)
    vp<-sum((pheno-mean(pheno))^2)/nrow(pheno)
    #goodness<-drop(cor(y,yhat)^2)
    pred_mse<-1-mse.cv/vp
    pred_r2<-pearson(geno,pheno,foldid)
    R2<- pred_r2[[1]]
    LASSO <- sqrt(R2)
    
    ypp2<-cbind(ypp2,pred_r2[[2]]$yyhat)
    yoo2<-cbind(yoo2,pred_r2[[2]]$yyobs)
    
    
    ############################### PLS ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        pls.fit <- plsr(y1~x1,ncomp=5,validation="CV")
        #plot(RMSEP(pls.fit), legendpos = "topright")
        nn <-as.numeric( which.min(tt <-RMSEP(pls.fit)$val[1,,][-1]))
        ## predicting
        yhat<-predict(pls.fit, newdata=x2,ncomp=nn)
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    PLS<-cor(yo,yp)
    
    ypp3<-cbind(ypp3,yp)
    yoo3<-cbind(yoo3,yo)
    
    
    ############################### SSVS ################################
    
    x<-rep(1,length(pheno))
    x<-as.matrix(x)
    
    ETA<-list( list(X=x,model='FIXED'),
               list(X=geno, model='BayesB'))
    fm<-BGLR(y=pheno,ETA=ETA, nIter=1500, burnIn=500)
    uHat<- fm$ETA[[2]]$b
    yHat<-fm$yHat
    r2.fit<-cor(y,yHat)^2
    
    
    yobs<-NULL
    yhat<-NULL
    id<-NULL
    #nfold<-6
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-pheno[i1,,drop=F]
        z1<-geno[i1,,drop=F]
        ETA<-list( list(X=x1,model='FIXED'),
                   list(X=z1, model='BayesB'))
        fm<-BGLR(y=y1,ETA=ETA, nIter=1500, burnIn=500)
        b<-fm$ETA[[1]]$b
        u <- fm$ETA[[2]]$b
        x2<-x[i2,,drop=F]
        y2<-pheno[i2,,drop=F]
        z2<-geno[i2,,drop=F]
        y3<-x2*b+z2%*%u+fm$mu
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        #fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    
    BayesB<-cor(yobs,yhat)
    
    ypp4<-cbind(ypp4,yhat)
    yoo4<-cbind(yoo4,yobs)
    
    
    ############################### SVM-RBF ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "rbfdot",kpar = "automatic",cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMRBF<-cor(yo,yp)
    
    ypp5<-cbind(ypp5,yp)
    yoo5<-cbind(yoo5,yo)
    
    
    ############################### SVM-POLY ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "poly",kpar = "automatic", cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMPOLY<-cor(yo,yp)
    
    ypp6<-cbind(ypp6,yp)
    yoo6<-cbind(yoo6,yo)
    
    result<- data.frame(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
    res<-rbind(res,result)
}

mean<- apply(res,2,mean)
sd<- apply(res,2,sd)
out<-rbind(mean,sd)
res2=rbind(res,out)
res2
write.table(res2,file=paste('Results/six_methods.RNAseq_and_methylation.', trait, '.txt', sep=''),
            sep='\t', quote=F, row.names=F)







############################## miRNAs + methylation #############################


geno <- as.matrix(t(rbind(mir, methy)))
geno <- scale(geno)
dim(geno)

res<-NULL
ypp1<-NULL
yoo1<-NULL
ypp2<-NULL
yoo2<-NULL
ypp3<-NULL
yoo3<-NULL
ypp4<-NULL
yoo4<-NULL
ypp5<-NULL
yoo5<-NULL
ypp6<-NULL
yoo6<-NULL


for (z in 1:10) {
    cat ('============================================\n')
    cat (rep(c(z,'-'),10))
    foldid <- foldidID[[z]]
    
    
    ############################### BLUP ################################
    kk<-kinship(gen=geno)
    
    kk <- kk[[1]]
    kk<-kk[,-c(1,2)]
    kk<-as.matrix(kk)
    
    n<-length(pheno)
    x<-matrix(1,n,1)
    
    blup<-cv.mixed(x=x,y=pheno,kk=kk,nfold=nfold,foldid=foldid)
    br2<- as.numeric(blup[[1]])
    BLUP <- sqrt(br2)
    
    ypp1<-cbind(ypp1,blup[[2]]$yhat)
    yoo1<-cbind(yoo1,blup[[2]]$yobs)
    
    
    ############################### LASSO ################################
    
    
    fit<-cv.glmnet(x=geno,y=pheno,foldid=foldid)
    theta<-coef(fit,s="lambda.min")
    beta<-theta[1]
    gamma<-theta[-1]
    yhat<-predict(fit,newx=geno,s="lambda.min")
    lambda<-fit$lambda.min
    index<-match(lambda,fit$lambda)
    mse.cv<-fit$cvm[index]
    nonzero<-fit$nzero[index]
    sigma2<-sum((pheno-yhat)^2)/nrow(pheno)
    vp<-sum((pheno-mean(pheno))^2)/nrow(pheno)
    #goodness<-drop(cor(y,yhat)^2)
    pred_mse<-1-mse.cv/vp
    pred_r2<-pearson(geno,pheno,foldid)
    R2<- pred_r2[[1]]
    LASSO <- sqrt(R2)
    
    ypp2<-cbind(ypp2,pred_r2[[2]]$yyhat)
    yoo2<-cbind(yoo2,pred_r2[[2]]$yyobs)
    
    
    ############################### PLS ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        pls.fit <- plsr(y1~x1,ncomp=5,validation="CV")
        #plot(RMSEP(pls.fit), legendpos = "topright")
        nn <-as.numeric( which.min(tt <-RMSEP(pls.fit)$val[1,,][-1]))
        ## predicting
        yhat<-predict(pls.fit, newdata=x2,ncomp=nn)
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    PLS<-cor(yo,yp)
    
    ypp3<-cbind(ypp3,yp)
    yoo3<-cbind(yoo3,yo)
    
    
    ############################### SSVS ################################
    
    x<-rep(1,length(pheno))
    x<-as.matrix(x)
    
    ETA<-list( list(X=x,model='FIXED'),
               list(X=geno, model='BayesB'))
    fm<-BGLR(y=pheno,ETA=ETA, nIter=1500, burnIn=500)
    uHat<- fm$ETA[[2]]$b
    yHat<-fm$yHat
    r2.fit<-cor(y,yHat)^2
    
    
    yobs<-NULL
    yhat<-NULL
    id<-NULL
    #nfold<-6
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-pheno[i1,,drop=F]
        z1<-geno[i1,,drop=F]
        ETA<-list( list(X=x1,model='FIXED'),
                   list(X=z1, model='BayesB'))
        fm<-BGLR(y=y1,ETA=ETA, nIter=1500, burnIn=500)
        b<-fm$ETA[[1]]$b
        u <- fm$ETA[[2]]$b
        x2<-x[i2,,drop=F]
        y2<-pheno[i2,,drop=F]
        z2<-geno[i2,,drop=F]
        y3<-x2*b+z2%*%u+fm$mu
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        #fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    
    BayesB<-cor(yobs,yhat)
    
    ypp4<-cbind(ypp4,yhat)
    yoo4<-cbind(yoo4,yobs)
    
    
    ############################### SVM-RBF ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "rbfdot",kpar = "automatic",cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMRBF<-cor(yo,yp)
    
    ypp5<-cbind(ypp5,yp)
    yoo5<-cbind(yoo5,yo)
    
    
    ############################### SVM-POLY ################################
    
    yp<-NULL
    yo<-NULL
    for(k in 1:nfold){
        id1<-which(foldid!=k)
        id2<-which(foldid==k)
        x1<-geno[id1,]
        x2<-geno[id2,]
        y1<-pheno[id1]
        y2<-pheno[id2]
        kern<- ksvm(x=x1,y=y1,type='eps-svr',kernel = "poly",kpar = "automatic", cross=0)
        ## predicting
        yhat<-predict(kern,x2,type='decision')
        yp<-c(yp,yhat)
        yo<-c(yo,y2)
    }
    SVMPOLY<-cor(yo,yp)
    
    ypp6<-cbind(ypp6,yp)
    yoo6<-cbind(yoo6,yo)
    
    result<- data.frame(BLUP, LASSO, PLS, BayesB, SVMRBF, SVMPOLY)
    res<-rbind(res,result)
}

mean<- apply(res,2,mean)
sd<- apply(res,2,sd)
out<-rbind(mean,sd)
res2=rbind(res,out)
res2
write.table(res2,file=paste('Results/six_methods.miRNAs_and_methylation.', trait, '.txt', sep=''),
            sep='\t', quote=F, row.names=F)

