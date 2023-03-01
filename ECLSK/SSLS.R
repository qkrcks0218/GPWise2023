############################################
# Last update : Dec 19 2021
# Code for Calculating Estimators
############################################

PLM <- function(M,Y,A,PS,gX,g1,g0,adj=FALSE,target="ATE",PSrange=c(0.0,1.0),CID=NULL){
    
    G <- dim(M)[2]
    
    pos <- which(PSrange[1]<=PS&PS<=PSrange[2])
    if(G==1){
        M <- matrix(M[pos,],length(pos),1)
    } else {
        M <- M[pos,]
    }
    Y <- Y[pos]
    A <- A[pos]
    PS <- PS[pos]
    gX <- gX[pos]
    g1 <- g1[pos]
    g0 <- g0[pos]
    CID <- CID[pos]
    
    Est <- list()
    Sigma <- list()
    
    
    y1 <- Y - gX
    if(target=="ATE"){
        phi <- A*(Y-g1)/PS - (1-A)*(Y-g0)/(1-PS) + (g1-g0)
    } else if (target=="ATT"){
        phi <- (A - (1-A)*PS/(1-PS))*(Y-g0)/mean(A)
    }
    
    a1 <- (A - PS)*M
    a2 <- M
    LM1 <- lm(y1~0+a1)
    LM2 <- lm(phi~0+a2)
    XX1 <- t(a1)%*%(a1)
    XX2 <- t(a2)%*%(a2)
    
    if(is.null(CID)){
        MEAT1 <- t(a1)%*%diag(LM1$residuals^2)%*%(a1)
        MEAT2 <- t(a2)%*%diag(LM2$residuals^2)%*%(a2)
        MEAT12 <- t(a1)%*%diag(LM1$residuals*LM2$residuals)%*%(a2)
        V1 <- solve(XX1)%*%MEAT1%*%solve(XX1)
        V2 <- solve(XX2)%*%MEAT2%*%solve(XX2)
        V12 <- solve(XX1)%*%MEAT12%*%solve(XX2)
        
    } else {
        MEAT1 <- MEAT2 <- MEAT12 <- matrix(0,G,G)
        for(cid in unique(CID)){
            A1 <- matrix(a1[cid==CID,],sum(cid==CID),G)
            A2 <- matrix(a2[cid==CID,],sum(cid==CID),G)
            res1 <- matrix(LM1$residuals[cid==CID],sum(cid==CID),1)
            res2 <- matrix(LM2$residuals[cid==CID],sum(cid==CID),1)
            
            MEAT1  <- MEAT1 +t(A1)%*%t(t(res1))%*%t(res1)%*%(A1)
            MEAT2  <- MEAT2 +t(A2)%*%t(t(res2))%*%t(res2)%*%(A2)
            MEAT12 <- MEAT12+t(A1)%*%t(t(res1))%*%t(res2)%*%(A2)
        }
        V1 <- solve(XX1)%*%MEAT1%*%solve(XX1)
        V2 <- solve(XX2)%*%MEAT2%*%solve(XX2)
        V12 <- solve(XX1)%*%MEAT12%*%solve(XX2)
    }
    
    
    
    RESULT <- list()
    
    RESULT$Est <- matrix(c(LM1$coefficients,
                           LM2$coefficients),2*G,1)
    
    RESULT$Sigma <- matrix(0,2*G,2*G)
    RESULT$Sigma[1:G,1:G] <- V1
    RESULT$Sigma[1:G,G+1:G] <- V12
    RESULT$Sigma[G+1:G,1:G] <- t(V12)
    RESULT$Sigma[G+1:G,G+1:G] <- V2
    
    rownames(RESULT$Est) <- rownames(RESULT$Sigma) <- colnames(RESULT$Sigma) <- 
        c(paste("SSLS_",1:G,sep=""),paste("EIF_",1:G,sep=""))
    
    RESULT$N <- length(Y)
    
    return(RESULT)
}


LinComb <- function(M,PLM.EFF,w.range=c(0,1)){
    
    G <- dim(M)[2]
    N <- PLM.EFF$N
    EST <- PLM.EFF$Est
    Sigma <- PLM.EFF$Sigma
    
    
    S11 <- Sigma[1:G,1:G]
    S12 <- Sigma[1:G,G+(1:G)]
    S21 <- Sigma[G+(1:G),1:G]
    S22 <- Sigma[G+(1:G),G+(1:G)]
    
    output <- list()
    output$N <- apply(M,2,sum)
    output$EST <- matrix(0,3,G)
    output$Sigma <- matrix(0,3,G)
    
    w <- w2 <- rep(0,G)
    for(gg in 1:G){
        
        OPT <- optimize(function(w){ matrix(c(w,1-w),1,2)%*%Sigma[c(gg,gg+G),c(gg,gg+G)]%*%matrix(c(w,1-w),2,1) },
                        lower=w.range[1],upper=w.range[2])
        w[gg] <- (Sigma[gg+G,gg+G]-Sigma[gg,gg+G])/(Sigma[gg+G,gg+G]-2*Sigma[gg,gg+G]+Sigma[gg,gg])
        w[gg] <- as.numeric(w[gg]<w.range[1])*w.range[1] + 
            as.numeric( w.range[1] <= w[gg] & w[gg] <= w.range[2])*w[gg] + 
            as.numeric(w[gg]>w.range[2])*w.range[2]
    }
    
    W <- cbind(diag(w),diag(1-w))
    output$EST <- matrix(rbind(EST,W%*%EST),3,G,byrow=T)
    output$Sigma <- matrix(sqrt(c( diag(Sigma), diag(W%*%Sigma%*%t(W)))),3,G,byrow=T)
    output$FullSigma <- W%*%Sigma%*%t(W)
    output$W <- W
    
    return(output)
}













