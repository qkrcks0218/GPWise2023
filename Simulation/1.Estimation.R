############################################
# Last update : Dec 19 2021
# Code for Section 4
############################################

############################################
# Load Source Files and Packages
############################################

source("../MySL.R")
source("../SSLS.R")

library(grf)
library(dplyr)

############################################
# Parameter Setups
############################################

N <- 2000
TotalIter <- 5
SS.Total <- 5

ParaGrid <- expand.grid(c(1:100),c(0,1,2,3),c(0,0.75),c("V","C")) 

SL.hpara <- list()                          # Super Learner parameters
SL.hpara$SLL <- 1# c(1,2,3,4,5,6,7,9,10)
# Superlearner basic learning algorithms: 
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 9: gbm
# 10: 1-layer MLP
SL.hpara$MTRY <- c(2,4)                # random forest parameters
SL.hpara$MLPL <- c(2,4)                # number of nodes in MLP
SL.hpara$NMN <- 25                     # gbm parameter
SL.hpara$MLPdecay <- 10^c(-3,-4,-5)    # MLP decay parameter

pos.AX <- 2:6    # position of (A,X)
pos.X <- 3:6     # position of X

############################################
# Recommend to use parallel computing across BATCH in 1:1600
############################################

for(BATCH in 1:1600){   
    
    iteration <- ParaGrid[BATCH,1]
    beta <- ParaGrid[BATCH,2]
    rho <- ParaGrid[BATCH,3]
    PST <- ParaGrid[BATCH,4]
    sA <- 5*as.numeric(rho>0)*as.numeric(PST=="C") + 
        2.5*as.numeric(rho>0)*as.numeric(PST=="V")
    ClSize <- 10
    OT <- gaussian()
    
    
    Thr <- c(-Inf,seq(-1,1,1),Inf)
    G <- length(Thr)-1
    
    Mf <- function(x,jj){ as.numeric( Thr[jj] <= x & x < Thr[jj+1] )} 
    
    expit <- function(x){ exp(x)/(1+exp(x)) }
    
    tau <- function(x){
        (beta/8)*abs(x^2)*as.numeric(x> -1) + (beta*abs(x)/4-beta/8)*as.numeric(x<= -1) + 1
    }
    
    PS.Ing <- function(x){
        if( is.null(dim(x)) ){
            lps <- 0.25*x[1] + 0.3*x[2]^2 - 0.3*abs(x[3])*x[4]
        } else {
            lps <- 0.25*x[,1] + 0.3*x[,2]^2 - 0.3*abs(x[,3])*x[,4]
        }
        return( lps )
    }
    
    
    
    ps <- function(x,type="V",sA){
        if(type=="V"){
            if(is.null(dim(x))){
                lps <- PS.Ing(x[1:4]) + sA*x[5]
            } else {
                lps <- PS.Ing(x[,1:4]) + sA*x[,5]
            }
            exp(lps)/(1+exp(lps))
        } else {
            if(is.null(dim(x))){
                lps <- 0.5+sA*x[5]
            } else {
                lps <- 0.5+sA*x[,5]
            }
            exp(lps)/(1+exp(lps))
        }
    }
    
    ps.int <- function(x,type="V",sA){
        if(type=="V"){
            if(sA>0){
                ps.temp <- function(u){
                    if(is.null(dim(x))){
                        lps <- PS.Ing(x[1:4]) + u
                    } else {
                        lps <- PS.Ing(x[,1:4]) + u
                    }
                    exp(lps)/(1+exp(lps))*dnorm(u,mean=0,sd=sA)
                }
                return( integrate(ps.temp,-10,10)$value )
            } else {
                if(is.null(dim(x))){
                    lps <- PS.Ing(x[1:4])
                } else {
                    lps <- PS.Ing(x[,1:4])
                }
                return( exp(lps)/(1+exp(lps)) )
            }
        } else {
            if(sA>0){
                ps.temp <- function(u){
                    if(is.null(dim(x))){
                        lps <- 0.5 + u
                    } else {
                        lps <- rep(0.5,length(x[,1])) + u
                    }
                    exp(lps)/(1+exp(lps))*dnorm(u,mean=0,sd=sA)
                }
                return( integrate(ps.temp,-10,10)$value )
            } else {
                if(is.null(dim(x))){
                    lps <- 0.5
                } else {
                    lps <- rep(0.5,length(x[,1]))
                }
                return( exp(lps)/(1+exp(lps)) )
            }
            return(lps)
        }
    }
    
    
    G0 <- function(x){
        X1 <- x[,1] ; X2 <- x[,2] ; X3 <- x[,3] ; X4 <- x[,4] ; U <- x[,5]
        (X1-0.5*X2^2+X3*X4)
    }
    G1 <- function(x){
        X1 <- x[,1] ; X2 <- x[,2] ; X3 <- x[,3] ; X4 <- x[,4] ; U <- x[,5]
        (X1-0.5*X2^2+X3*X4) + tau(X1)
    }
    tauH <- function(x){ 
        G1(x)-G0(x)
    }
    
    
    ############################################
    # True Groupwise Effects 
    ############################################
    
    N0 <- 10^7
    CS <- 1*as.numeric(rho==0) + ClSize*as.numeric(rho>0)
    C0 <- N0/CS
    X1 <- rnorm(N0)
    X2 <- rnorm(N0)
    X3 <- rep(rnorm(C0),each=CS)
    X4 <- rbinom(N0,1,0.5)
    U <-  rep(rnorm(C0),each=CS)
    V <-  rep(rnorm(C0),each=CS)
    PS <- ps(cbind(X1,X2,X3,X4,V),type=PST,sA)
    
    tg <- tov <- rep(0,length(Thr)-1)
    g1 <- G1(cbind(X1,X2,X3,X4,U))
    g0 <- G0(cbind(X1,X2,X3,X4,U))
    TH <- tauH(cbind(X1,X2,X3,X4,U))
    for(jj in 1:(length(Thr)-1)){
        tg[jj] <- sum(Mf(X1,jj)*TH)/sum(Mf(X1,jj))
        tov[jj] <- mean(Mf(X1,jj)*PS*(1-PS)*TH)/mean(Mf(X1,jj)*PS*(1-PS))
        print(jj)
    }
    
    cbind(tg,tov)
    
    rm(list=c("X1","X2","X3","X4","PS"))
    
    ###########################
    # DGP
    ###########################
    
    EST  <- BIAS  <- SE  <- COVER <-  COVER.S <-  matrix(0,TotalIter,4*G)
    EST2 <- BIAS2 <- SE2 <- COVER2 <- COVER2.S <- matrix(0,TotalIter,4*G)
    EST.CL  <- BIAS.CL  <- SE.CL  <- COVER.CL <-  COVER.S.CL <-  matrix(0,TotalIter,4*G)
    EST2.CL <- BIAS2.CL <- SE2.CL <- COVER2.CL <- COVER2.S.CL <- matrix(0,TotalIter,4*G)
    
    for(kkk in 1:TotalIter){
        
        set.seed(beta*10000+iteration*100+kkk)
        
        X1 <- rnorm(N)
        X2 <- rnorm(N)
        X3 <- rep(rnorm(N/CS),each=CS)
        X4 <- rbinom(N,1,0.5)
        U  <- rep(rnorm(N/CS),each=CS)
        V  <- rep(rnorm(N/CS),each=CS)
        
        PS <- ps(cbind(X1,X2,X3,X4,V),type=PST,sA)
        PS.noV <- apply(cbind(X1,X2,X3,X4),1,function(t){ps.int(t,type=PST,sA)})
        
        
        A <- rbinom(N,1,PS)
        CID <- rep(1:(N/CS),each=CS)
        
        tauX <- tauH( cbind(X1,X2,X3,X4,U) )
        g0 <- G0( cbind(X1,X2,X3,X4,U) )
        g1 <- g0 + tauX
        gX <- g0 + tauX*PS.noV
        
        pY <- g0+tauX*A
        
        Y <- g0 + tauX*A + rho*U + rnorm(N)
        
        M <- cbind(Mf(X1,1),Mf(X1,2))
        for(jj in 3:(length(Thr)-1)){
            M <- cbind(M,Mf(X1,jj))
        }
        
        PS.Original <- PS.noV
        gX.Original <- gX
        g1.Original <- g1
        g0.Original <- g0
        
        
        ############# EST
        
        Data <- data.frame( cbind(Y,A,X1,X2,X3,X4) )
        colnames(Data) <- c("Y","A","X1","X2","X3","X4")
        
        SS.Est <- list()
        SS.Var <- list()
        
        for(gg in 1:G){
            SS.Est[[gg]] <- matrix(0,SS.Total,2)
            SS.Var[[gg]] <- list()
        }
        
        EST.temp <- matrix(0,SS.Total,2*G)
        VAR.temp <- matrix(0,SS.Total,4*G^2)
        
        EST.temp.CL <- matrix(0,SS.Total,2*G)
        VAR.temp.CL <- matrix(0,SS.Total,4*G^2)
        
        for(ss.iter in 1:SS.Total){
            
            MSc <- sort(sample(unique(CID),length(unique(CID))/2))
            ASc <- unique(CID)[-MSc]
            
            DDD <- data.frame(id=1:N,cid=CID)
            
            DDD %>% filter(cid %in% MSc) -> DDD.temp
            MS.Cl <- DDD.temp$id      # split sample 1
            AS.Cl <- (1:N)[-MS.Cl]    # split sample 2
            
            OR.Fit <- OR.Fit.X <- PS.Fit <- list()
            
            OR.Fit$MS <-   # Estimation of E(Y|A,X) using split sample 1
                MySL(Data[MS.Cl,], 1, pos.AX, Ydist=OT,
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
            OR.Fit$AS <-   # Estimation of E(Y|A,X) using split sample 2
                MySL(Data[AS.Cl,], 1, pos.AX, Ydist=OT,
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
            
            OR.Fit.X$MS <- # Estimation of E(Y|X) using split sample 1
                MySL(Data[MS.Cl,], 1, pos.X, Ydist=OT,
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
            OR.Fit.X$AS <- # Estimation of E(Y|X) using split sample 2
                MySL(Data[AS.Cl,], 1, pos.X, Ydist=OT,
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
            
            OR.EST <- matrix(0,N,4)
            
            TData <- Data[MS.Cl,]
            OR.EST[MS.Cl,3] <-  # E(Y|A,X) estimates in split sample 2
                predict(OR.Fit$AS,newdata=TData[,pos.AX],onlySL=TRUE)$pred
            OR.EST[MS.Cl,4] <-  # E(Y|X) estimates in split sample 2
                predict(OR.Fit.X$AS,newdata=TData[,pos.X],onlySL=TRUE)$pred
            
            TData <- Data[MS.Cl,] ; TData$A <- 1
            OR.EST[MS.Cl,1] <-  # E(Y|A=1,X) estimates in split sample 2
                predict(OR.Fit$AS,newdata=TData[,pos.AX],onlySL=TRUE)$pred
            
            TData <- Data[MS.Cl,] ; TData$A <- 0
            OR.EST[MS.Cl,2] <-   # E(Y|A=0,X) estimates in split sample 2
                predict(OR.Fit$AS,newdata=TData[,pos.AX],onlySL=TRUE)$pred
            
            TData <- Data[AS.Cl,]
            OR.EST[AS.Cl,3] <-  # E(Y|A,X) estimates in split sample 1
                predict(OR.Fit$MS,newdata=TData[,pos.AX],onlySL=TRUE)$pred
            OR.EST[AS.Cl,4] <-  # E(Y|X) estimates in split sample 1
                predict(OR.Fit.X$MS,newdata=TData[,pos.X],onlySL=TRUE)$pred
            
            TData <- Data[AS.Cl,] ; TData$A <- 1
            OR.EST[AS.Cl,1] <-  # E(Y|A=1,X) estimates in split sample 1
                predict(OR.Fit$MS,newdata=TData[,pos.AX],onlySL=TRUE)$pred
            
            TData <- Data[AS.Cl,] ; TData$A <- 0
            OR.EST[AS.Cl,2] <-  # E(Y|A=0,X) estimates in split sample 1
                predict(OR.Fit$MS,newdata=TData[,pos.AX],onlySL=TRUE)$pred
            
            gX <- OR.EST[,4] ; g1 <- OR.EST[,1] ; g0 <- OR.EST[,2]
            
            if(PST=="V"){
                PS.Fit$MS <-  # P(A=1|X) estimation using split sample 1
                    MySL(Data[MS.Cl,], 2, pos.X, Ydist=binomial(),
                         SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
                PS.Fit$AS <-  # P(A=1|X) estimation using split sample 2
                    MySL(Data[AS.Cl,], 2, pos.X, Ydist=binomial(),
                         SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, NMN=SL.hpara$NMN,MLPdecay=SL.hpara$MLPdecay)
                PS.EST <- matrix(0,N,1)
                TData <- Data[MS.Cl,]
                PS.EST[MS.Cl,] <-  # P(A=1|X) estimates in split sample 2
                    predict(PS.Fit$AS,newdata=TData[,pos.X],onlySL=TRUE)$pred
                TData <- Data[AS.Cl,]
                PS.EST[AS.Cl,] <-  # P(A=1|X) estimates in split sample 1
                    predict(PS.Fit$MS,newdata=TData[,pos.X],onlySL=TRUE)$pred
                PS <- PS.EST
            } else {
                PS <- PS.EST <- PS.Original
            }
            
            
            
            
            
            PLM.Result.Original <- PLM(M,Y,A,PS,gX,g1,g0,adj=FALSE,target="ATE")
            
            EST.temp[ss.iter,] <- as.numeric( PLM.Result.Original$Est )
            VAR.temp[ss.iter,] <- as.numeric( PLM.Result.Original$Sigma )
            
            
            
            PLM.Result.Original <- PLM(M,Y,A,PS,gX,g1,g0,adj=FALSE,target="ATE",CID=CID)
            
            EST.temp.CL[ss.iter,] <- as.numeric( PLM.Result.Original$Est )
            VAR.temp.CL[ss.iter,] <- as.numeric( PLM.Result.Original$Sigma )
            
            
        }
        
        PLM.MEDIAN <- PLM.Result.Original
        PLM.MEDIAN$Est[,1] <- apply(EST.temp,2,median)
        
        
        dev.mat <- (EST.temp[,1:(2*G)]-matrix( PLM.MEDIAN$Est[,1],SS.Total,2*G, byrow=T ))
        NORM <- rep(0,SS.Total)
        for(ss.iter in 1:SS.Total){
            NORM[ss.iter] <- norm(t(t(dev.mat[ss.iter,]))%*%t(dev.mat[ss.iter,]) + matrix(VAR.temp[ss.iter,],2*G,2*G))
        }
        ss.iter <- which.min(abs(NORM-median(NORM)))
        PLM.MEDIAN$Sigma <- t(t(dev.mat[ss.iter,]))%*%t(dev.mat[ss.iter,]) + matrix(VAR.temp[ss.iter,],2*G,2*G)
        
        LinResult01 <- LinComb(M,PLM.MEDIAN,w.range=c(0,1))
        
        GRF <- causal_forest(cbind(X1,X2,X3,X4),Y,A,tune.parameters="all")
        ATE.GRF <- cbind( average_treatment_effect(GRF,subset=M[,1]==1),
                          average_treatment_effect(GRF,subset=M[,2]==1),
                          average_treatment_effect(GRF,subset=M[,3]==1),
                          average_treatment_effect(GRF,subset=M[,4]==1) )
        
        
        PLM.MEDIAN.CL <- PLM.Result.Original
        PLM.MEDIAN.CL$Est[,1] <- apply(EST.temp.CL,2,median)
        
        
        dev.mat.CL <- (EST.temp.CL[,1:(2*G)]-matrix( PLM.MEDIAN.CL$Est[,1],SS.Total,2*G, byrow=T ))
        NORM <- rep(0,SS.Total)
        for(ss.iter in 1:SS.Total){
            NORM[ss.iter] <- norm(t(t(dev.mat.CL[ss.iter,]))%*%t(dev.mat.CL[ss.iter,]) + matrix(VAR.temp.CL[ss.iter,],2*G,2*G))
        }
        ss.iter <- which.min(abs(NORM-median(NORM)))
        PLM.MEDIAN.CL$Sigma <- t(t(dev.mat.CL[ss.iter,]))%*%t(dev.mat.CL[ss.iter,]) + matrix(VAR.temp.CL[ss.iter,],2*G,2*G)
        
        LinResult201 <- LinComb(M,PLM.MEDIAN.CL,w.range=c(0,1))
        
        GRF.CL <- causal_forest(cbind(X1,X2,X3,X4),Y,A,clusters=CID,tune.parameters="all")
        ATE.GRF.CL <- cbind( average_treatment_effect(GRF.CL,subset=M[,1]==1),
                             average_treatment_effect(GRF.CL,subset=M[,2]==1),
                             average_treatment_effect(GRF.CL,subset=M[,3]==1),
                             average_treatment_effect(GRF.CL,subset=M[,4]==1) )
        
        
        
        EST2[kkk,1:(3*G)] <- as.numeric( t(LinResult01$EST) )
        EST2[kkk,1:G+(3*G)] <- ATE.GRF[1,]
        BIAS2[kkk,] <- EST2[kkk,] - rep(tg,4)
        
        SE2[kkk,1:(3*G)] <- as.numeric( t(LinResult01$Sigma) )
        SE2[kkk,1:G+(3*G)] <- ATE.GRF[2,]
        
        for(ggg in 1:G){
            for(estt in 0:3){
                COVER2[kkk,ggg+G*estt] <- 
                    as.numeric( EST2[kkk,ggg+(estt*G)]-qnorm(0.975)*SE2[kkk,ggg+(estt*G)] <= tg[ggg] &
                                    EST2[kkk,ggg+(estt*G)]+qnorm(0.975)*SE2[kkk,ggg+(estt*G)] >= tg[ggg] )
                COVER2.S[kkk,ggg+G*estt] <- 
                    as.numeric( EST2[kkk,ggg+(estt*G)]-qnorm(1-0.05/2/G)*SE2[kkk,ggg+(estt*G)] <= tg[ggg] &
                                    EST2[kkk,ggg+(estt*G)]+qnorm(1-0.05/2/G)*SE2[kkk,ggg+(estt*G)] >= tg[ggg] )
            }
            
        }
        
        
        EST2.CL[kkk,1:(3*G)] <- as.numeric( t(LinResult201$EST) )
        EST2.CL[kkk,1:G+(3*G)] <- ATE.GRF.CL[1,]
        BIAS2.CL[kkk,] <- EST2.CL[kkk,] - rep(tg,4)
        
        SE2.CL[kkk,1:(3*G)] <- as.numeric( t(LinResult201$Sigma) )
        SE2.CL[kkk,1:G+(3*G)] <- ATE.GRF.CL[2,]
        
        for(ggg in 1:G){
            for(estt in 0:3){
                COVER2.CL[kkk,ggg+G*estt] <- 
                    as.numeric( EST2.CL[kkk,ggg+(estt*G)]-qnorm(0.975)*SE2.CL[kkk,ggg+(estt*G)] <= tg[ggg] &
                                    EST2.CL[kkk,ggg+(estt*G)]+qnorm(0.975)*SE2.CL[kkk,ggg+(estt*G)] >= tg[ggg] )
                COVER2.S.CL[kkk,ggg+G*estt] <- 
                    as.numeric( EST2.CL[kkk,ggg+(estt*G)]-qnorm(1-0.05/2/G)*SE2.CL[kkk,ggg+(estt*G)] <= tg[ggg] &
                                    EST2.CL[kkk,ggg+(estt*G)]+qnorm(1-0.05/2/G)*SE2.CL[kkk,ggg+(estt*G)] >= tg[ggg] )
            }
            
        }
        
        
        print(kkk)
        
    }
    
    
    
    
    
    RRR2 <- cbind(BIAS2,SE2,COVER2,COVER2.S)
    RRR2.CL <- cbind(BIAS2.CL,SE2.CL,COVER2.CL,COVER2.S.CL)
    colnames(RRR2) <- colnames(RRR2.CL) <- c(paste("Bias_GP_",1:G,sep=""),
                                             paste("Bias_DR_",1:G,sep=""),
                                             paste("Bias_Lin01_",1:G,sep=""),
                                             paste("Bias_GRF_",1:G,sep=""),
                                             paste("SE_GP_",1:G,sep=""),
                                             paste("SE_DR_",1:G,sep=""),
                                             paste("SE_Lin01_",1:G,sep=""),
                                             paste("SE_GRF_",1:G,sep=""),
                                             paste("Cover_GP_",1:G,sep=""),
                                             paste("Cover_DR_",1:G,sep=""),
                                             paste("Cover_Lin01_",1:G,sep=""),
                                             paste("Cover_GRF_",1:G,sep=""),
                                             paste("Cover_Simul_GP_",1:G,sep=""),
                                             paste("Cover_Simul_DR_",1:G,sep=""),
                                             paste("Cover_Simul_Lin01_",1:G,sep=""),
                                             paste("Cover_Simul_GRF_",1:G,sep=""))
    
    write.csv(RRR2,sprintf("Result_%sPS_%sOR_Est_rho%0.3d_B%0.4d_SubB%0.4d.csv",
                           PST,OST,round(rho*100),beta,iteration),row.names=FALSE)
    write.csv(RRR2.CL,sprintf("Result_%sPS_%sOR_CL_rho%0.3d_B%0.4d_SubB%0.4d.csv",
                              PST,OST,round(rho*100),beta,iteration),row.names=FALSE)
    
}
