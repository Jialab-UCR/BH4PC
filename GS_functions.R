

####################### for BLUP ###########################


mixed <- function(x,y,kk,method="REML",eigen=FALSE){
    
    loglike<-function(theta){
        lambda<-exp(theta)
        logdt<-sum(log(lambda*delta+1))
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,s,1)
        xx<-matrix(0,s,s)
        for(i in 1:s){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:s){
                xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
        }
        xx
        if(method=="REML"){
            loglike<- -0.5*logdt-0.5*(n-s)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
        } else {
            loglike<- -0.5*logdt-0.5*n*log(yy-t(yx)%*%solve(xx)%*%yx)
        }
        return(-loglike)
    }
    
    
    fixed<-function(lambda){
        h<-1/(lambda*delta+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,s,1)
        xx<-matrix(0,s,s)
        for(i in 1:s){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:s){
                xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
        } 
        beta<-solve(xx,yx)
        if(method=="REML"){
            sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-s)
        } else {
            sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/n
        }
        var<-diag(solve(xx)*drop(sigma2))
        stderr<-sqrt(var)
        return(c(beta,stderr,sigma2))
    }
    
    n<-length(y)
    qq<-eigen(kk,symmetric=TRUE)
    delta<-qq[[1]]
    uu<-qq[[2]]
    s<-ncol(x)
    yu<-t(uu)%*%y
    xu<-t(uu)%*%x
    theta<-0
    #parm<-optim(par=theta,fn=loglike,NULL,hessian = TRUE, method="L-BFGS-B",lower=-10,upper=10)
    parm<-optim(par=theta,fn=loglike,NULL,hessian = TRUE, method="Brent",lower=-10,upper=10)
    lambda<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    hess<-parm$hessian
    parmfix<-fixed(lambda)
    beta<-parmfix[1:s]
    stderr<-parmfix[(s+1):(2*s)]
    ve<-parmfix[2*s+1]
    lod<-lrt/4.61
    p_value<-1-pchisq(lrt,1)
    va<-lambda*ve
    h2<-va/(va+ve)
    par<-data.frame(method,beta,stderr,va,ve,lambda,h2,conv,fn1,fn0,lrt,lod,p_value)
    if(eigen){
        return(list(par,qq))
    } else {
        return(list(par))
    }
}


kinship <- function(gen){
    m<-ncol(gen)
    n<-nrow(gen)

    kk <- crossprod(t(gen))
    
    #kk<-matrix(0,n,n)
    #for(k in 1:m){
    #    z<-gen[,k,drop=F]
    #    kk<-kk+z%*%t(z)
    #}
    cc<-mean(diag(kk))
    kk<-kk/cc
    parm<-matrix(1,n,1)
    row<-matrix(1:n,n,1)
    kk<-data.frame(parm,row,kk)
    aa<-c("parm","row")
    for(k in 1:n){
        a<-paste("col",k,sep="")
        aa<-c(aa,a)
    }
    names(kk)<-aa
    cc<-data.frame(cc)
    return(list(kk,cc))
}




cv.mixed<-function(x,y,kk,nfold=nfold,foldid=NULL){
    n<-length(y)
    #foldid<-read.csv (file="yan\\foldid.csv",row.names= 1)
    #if(length(foldid)>0){
    #nfold<-max(foldid)
    #}  else {
    # foldid<-sample(rep(1:nfold,ceiling(n/nfold))[1:n])
    #}
    yobs<-NULL
    yhat<-NULL
    fold<-NULL
    id<-NULL
    for(k in 1:nfold){
        i1<-which(foldid!=k)
        i2<-which(foldid==k)
        x1<-x[i1,,drop=F]
        y1<-y[i1,,drop=F]
        k11<-kk[i1,i1]
        parm<-mixed(x=x1,y=y1,kk=k11,method="REML",eigen=TRUE)
        qq<-parm[[2]]
        delta<-qq[[1]]
        u<-qq[[2]]
        beta<-as.matrix(parm[[1]]$beta,ncol(x),1)
        va<-parm[[1]]$va[1]
        ve<-parm[[1]]$ve[1]
        x2<-x[i2,,drop=F]
        y2<-y[i2,,drop=F]
        k21<-kk[i2,i1]
        h<-1/(delta*va+ve)
        y3<-x2%*%beta+va*k21%*%u%*%diag(h)%*%t(u)%*%(y1-x1%*%beta)
        #y3<-x2%*%beta+va*k21%*%solve(k11*va+diag(length(y1))*ve)%*%(y1-x1%*%beta)
        fold<-c(fold,rep(k,length(y2)))
        yobs<-c(yobs,y2)
        yhat<-c(yhat,y3)
        id<-c(id,i2)
    }
    pred<-data.frame(fold,id,yobs,yhat)
    cv.r2<-cor(pred$yobs,pred$yhat)^2 
    return(list(data.frame(r2=cv.r2),pred))
}


####################################


