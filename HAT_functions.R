gsm.1d.simple <- function(mydata, mykin){
    
    ### Negative nloglikelihood Function
    nloglik.REML.1d <- function(theta){
        
        lambda<-exp(theta)
        logdt<-sum(log(lambda*vv+1))
        h<-1/(lambda*vv+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,s,1)
        xx<-matrix(0,s,s)
        for(i in 1:s){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:s){
                xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
        }
        loglike<- -0.5*logdt-0.5*(n-s)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
        return(-loglike)
    }
    
    
    ### Henderson Equation and Estimation
    henderson.FUN <- function(t, y, kin){
        
        n <- length(y)
        x.design <- matrix(1, n, 1)
        
        phi <- t[1]
        sigma <- t[2]
        lambda <- phi/sigma
        
        TL <- t(x.design)%*%x.design
        BL <- x.design
        TR <- t(x.design)
        
        vv<-eigen(kin)[[1]]
        thres <- min(abs(vv))
        
        if (thres > 1e-8) {
            BR <- diag(n) + solve(kin)/lambda
        } else {
            BR <- diag(n) + solve(kin+diag(1e-8,nrow(kin)))/lambda
        }
        
        #BR <- diag(n) + ginv(kin)/lambda ## hat value reported is not reasonable
        
        v.est <- solve(cbind(rbind(TL, BL), rbind(TR, BR)))
        est <- v.est%*%matrix(c(t(x.design)%*%y, y), n+1, 1)
        
        beta_hat <- est[1, 1]
        eta_hat <- est[-1, 1]
        
        return(list(beta=beta_hat, eta=eta_hat, v=v.est))
    }
    
    
    x<-matrix(1,length(mydata),1)
    s<-ncol(x)
    n<-length(mydata)
    uu<-eigen(mykin)[[2]]
    vv<-eigen(mykin)[[1]]
    xu<-t(uu)%*%x
    yu<-t(uu)%*%mydata
    
    
    
    ### 1st search
    coef.space <- seq(-1, 1, 0.1)
    
    coef.comb <- NULL
    for(cg in coef.space){
        coef.comb <- rbind(coef.comb, cg)
    }
    #print(coef.comb)
    nloglik.REML.1d.batch <- function(v){
        theta <- v
        return(nloglik.REML.1d(theta))
    }
    
    nloglik.batch <- apply(coef.comb, 1, nloglik.REML.1d.batch)
    nloglik.optimal <- min(nloglik.batch)
    sel.id <- which(nloglik.batch==nloglik.optimal)[1]
    theta <- coef.comb[sel.id]
    junk <- optim(par=theta,fn=nloglik.REML.1d,NULL,hessian = TRUE, method="L-BFGS-B",lower=-10,upper=10)
    lambda <- exp(junk$par) 
    #print(lambda)
    h<-1/(lambda*vv+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,s,1)
    xx<-matrix(0,s,s)
    for(i in 1:s){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:s){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
    }
    
    
    sigma2 <- (yy-t(yx)%*%solve(xx)%*%yx)/(n-s)
    par <- c(lambda*sigma2, sigma2)
    var.para <- par
    #print(sigma2)
    
    junk2 <- henderson.FUN(par, mydata, mykin)
    beta.est <- junk2$beta
    eta.est <- junk2$eta
    v.est <- junk2$v
    
    
    #print("Starting HAT algorithm ......")
    
    v.phi <- var.para[1]*mykin
    v.sigma <- var.para[2]*diag(n)
    v <- v.phi+v.sigma
    
    hatMatrix <- v.phi%*%solve(v)
    residual <- mydata - rep(beta.est, n) - eta.est
    
    press <- 0
    for(i in 1:n){
        residual.tmp <- residual[i]
        hatMatrix.tmp <- hatMatrix[i, i]
        residual_pred.tmp <- residual.tmp/(1-hatMatrix.tmp)
        press <- press + residual_pred.tmp**2
    }
    
    ss <- sum((mydata-mean(mydata))**2)
    
    predic.HAT <- 1-press/ss
    
    res <- list(beta.est=beta.est, var.para=var.para, eta.est=eta.est, v.est=v.est, predic.HAT=predic.HAT)
    
    return(res)
    
}




###### simplified 1d function ######
gsm.1d.simple.m <- function(data, kin, foldid){
    p0 <- which(foldid==0)
    p1 <- which(foldid==1)
    mydata <- data[p0]
    mykin <- kin[p0,p0]
    
    
    ### Negative nloglikelihood Function
    nloglik.REML.1d <- function(theta){
        
        lambda<-exp(theta)
        logdt<-sum(log(lambda*vv+1))
        h<-1/(lambda*vv+1)
        yy<-sum(yu*h*yu)
        yx<-matrix(0,s,1)
        xx<-matrix(0,s,s)
        for(i in 1:s){
            yx[i]<-sum(yu*h*xu[,i])
            for(j in 1:s){
                xx[i,j]<-sum(xu[,i]*h*xu[,j])
            }
        }
        loglike<- -0.5*logdt-0.5*(n-s)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
        return(-loglike)
    }
    
    
    ### Henderson Equation and Estimation
    henderson.FUN <- function(t, y, kin){
        
        n <- length(y)
        x.design <- matrix(1, n, 1)
        
        phi <- t[1]
        sigma <- t[2]
        lambda <- phi/sigma
        
        TL <- t(x.design)%*%x.design
        BL <- x.design
        TR <- t(x.design)
        
        vv<-eigen(kin)[[1]]
        thres <- min(abs(vv))
        
        if (thres > 1e-8) {
            BR <- diag(n) + solve(kin)/lambda
        } else {
            BR <- diag(n) + solve(kin+diag(1e-8,nrow(kin)))/lambda
        }
        
        #BR <- diag(n) + ginv(kin)/lambda ## hat value reported is not reasonable
        
        v.est <- solve(cbind(rbind(TL, BL), rbind(TR, BR)))
        est <- v.est%*%matrix(c(t(x.design)%*%y, y), n+1, 1)
        
        beta_hat <- est[1, 1]
        eta_hat <- est[-1, 1]
        
        return(list(beta=beta_hat, eta=eta_hat, v=v.est))
    }
    
    

    x<-matrix(1,length(mydata),1)
    s<-ncol(x)
    n<-length(mydata)
    uu<-eigen(mykin)[[2]]
    vv<-eigen(mykin)[[1]]
    xu<-t(uu)%*%x
    yu<-t(uu)%*%mydata
    
    
    
    ### 1st search
    coef.space <- seq(-1, 1, 0.1)
    
    coef.comb <- NULL
    for(cg in coef.space){
        coef.comb <- rbind(coef.comb, cg)
    }
    #print(coef.comb)
    nloglik.REML.1d.batch <- function(v){
        theta <- v
        return(nloglik.REML.1d(theta))
    }
    
    nloglik.batch <- apply(coef.comb, 1, nloglik.REML.1d.batch)
    nloglik.optimal <- min(nloglik.batch)
    sel.id <- which(nloglik.batch==nloglik.optimal)[1]
    theta <- coef.comb[sel.id]
    junk <- optim(par=theta,fn=nloglik.REML.1d,NULL,hessian = TRUE, method="L-BFGS-B",lower=-10,upper=10)
    lambda <- exp(junk$par) 
    #print(lambda)
    h<-1/(lambda*vv+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,s,1)
    xx<-matrix(0,s,s)
    for(i in 1:s){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:s){
            xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
    }
    
    
    sigma2 <- (yy-t(yx)%*%solve(xx)%*%yx)/(n-s)
    par <- c(lambda*sigma2, sigma2)
    var.para <- par
    #print(sigma2)
    
    junk2 <- henderson.FUN(par, mydata, mykin)
    beta.est <- junk2$beta
    eta.est <- junk2$eta
    v.est <- junk2$v
    
    
    #print("Starting HAT algorithm ......")
    
    v.phi <- var.para[1]*mykin
    v.sigma <- var.para[2]*diag(n)
    v <- v.phi+v.sigma
    
    lambda <- (var.para[1])/(var.para[2])
    pred <- (lambda*kin[p1,p0]%*%uu%*%diag(h)%*%t(uu)%*%(mydata-rep(beta.est, length(p0)))) + rep(beta.est, length(p1))
    
    
    
    
    hatMatrix <- v.phi%*%solve(v)
    residual <- mydata - rep(beta.est, n) - eta.est
    
    press <- 0
    for(i in 1:n){
        residual.tmp <- residual[i]
        hatMatrix.tmp <- hatMatrix[i, i]
        residual_pred.tmp <- residual.tmp/(1-hatMatrix.tmp)
        press <- press + residual_pred.tmp**2
    }
    
    ss <- sum((mydata-mean(mydata))**2)
    
    predic.HAT <- 1-press/ss
    
    r2 <- cor(pred, data[p1])^2
    #print(r2)
    res <- list(beta.est=beta.est, var.para=var.para, eta.est=eta.est, v.est=v.est, predic.HAT=predic.HAT, pred=pred, r2=r2)
    
    return(res)
    
}




