############################################
# Load Source Files and Packages
############################################

source("SSLS.R")

TrueData <- read.csv("ECLSK.csv")
SSI.Index <- read.csv("SSInd_Index.csv")
SSC.Index <- read.csv("SSCL_Index.csv")
TrueData$S_age <- scale(TrueData$S_age)
TrueData$S_sestatus <- scale(TrueData$S_sestatus)
TrueData0 <- TrueData
TrueData0[,"A"] <- 0
TrueData1 <- TrueData
TrueData1[,"A"] <- 1
Y <- TrueData$Y
A <- TrueData$A
GP <- TrueData$GP
pos.X  <- 3:16
pos.AX <- 2:16
N <- dim(TrueData)[1]
X <- TrueData[,pos.X]

############################################
# Read NF
############################################

TNF <- 100
NData <- list()
for(ii in 1:TNF){
  NData[[ii]] <- read.csv(sprintf("NData/Result_BATCH%0.4d.csv",ii))
}

NF.Mat <- array(unlist(NData), c(dim(NData[[1]])[1], dim(NData[[1]])[2], TNF))
NF.median <- apply(NF.Mat,c(1:2),median) # 20-42
NF.median <- data.frame(NF.median)
colnames(NF.median) <- colnames(NData[[1]])

############################################
# Define Subgroups
############################################

CC <- C1 <- 1*TrueData$C_location_city +  2*TrueData$C_location_rural + (1-TrueData$C_location_city-TrueData$C_location_rural)*3
C1.Name <- c("Center City","Rural","Urban")
Gname <- apply(expand.grid(C1.Name),1,function(v){paste(v,collapse="\n")})
 

G <- length(unique(CC))
M <- matrix(0,N,G)

for(bb in 1:N){
  M[bb,CC[bb]] <- 1
}

############################################
# Semiparametric/Nonparametric Estimators
############################################

EFF.CL <- matrix(0,TNF,2*G)
VAR.CL <- matrix(0,TNF,4*G^2)

for(ii in 1:TNF){

  PLM.EFF2 <- PLM(M,
                  Y,
                  A,
                  PS=NData[[ii]]$PS,
                  gX=NData[[ii]]$ORX,
                  g1=NData[[ii]]$OR1,
                  g0=NData[[ii]]$OR0,
                  CID=GP)
  
  EFF.CL[ii,] <- t(PLM.EFF2$Est)
  VAR.CL[ii,] <- as.numeric(PLM.EFF2$Sigma)
  
  print(ii)
}

# write.csv(EFF.CL,"EFF_CL_3.csv",row.names=F)
# write.csv(VAR.CL,"VAR_CL_3.csv",row.names=F)

############################################
# Semiparametric/Nonparametric Estimators
############################################

# EFF.CL <- as.matrix(read.csv("EFF_CL_3.csv")[1:TNF,])
# VAR.CL <- as.matrix(read.csv("VAR_CL_3.csv")[1:TNF,])
for(bb in 1:dim(EFF.CL)[1]){
  EFF.CL[bb,] <- as.numeric(EFF.CL[bb,])
  VAR.CL[bb,] <- as.numeric(VAR.CL[bb,])
}

RES.CL <- (EFF.CL - matrix(apply(EFF.CL,2,median),TNF,2*G,byrow=T))
VAR.NEW.CL <- VAR.CL
MM <- matrix(0,2*G,2*G)
diag(MM) <- 1
for(gg in 1:G){
  MM[gg,gg+G] <- MM[G+gg,gg] <- 1
}

############################################
# Median adjustment across 100 cross-fitting estimates
############################################

for(jj in 1:TNF){
  VAR.NEW.CL[jj,] <- as.numeric( VAR.CL[jj,])  + 
    matrix(as.numeric(RES.CL[jj,]),2*G,1)%*%matrix(as.numeric(RES.CL[jj,]),1,2*G)
  print(jj)
}

NORM.CL <- apply(VAR.NEW.CL,1,function(x){norm(matrix(x,2*G,2*G),type="2")})

PLM.EFF.Median.CL <- PLM(M,
                         Y,
                         A,
                         PS=NData[[1]]$PS,
                         gX=NData[[1]]$ORX,
                         g1=NData[[1]]$OR1,
                         g0=NData[[1]]$OR0,
                         CID=GP)
PLM.EFF.Median.CL$Est[,1] <- as.numeric( apply(EFF.CL,2,median) )
PLM.EFF.Median.CL$Sigma[1:(2*G),1:(2*G)] <- matrix(VAR.NEW.CL[which.min(abs(median(NORM.CL)-NORM.CL)),],2*G,2*G)

############################################
# Obtain Weighted Estimate
############################################

LC <- LinComb(M,PLM.EFF.Median.CL)

############################################
# Assumption Check: Semiparametric Model
############################################
 
## Falsification Test
## overlap gp effect = gp effect, marginal

CONT <- cbind(diag(rep(1,G)),diag(rep(-1,G)))

(CONT%*%PLM.EFF.Median.CL$Est)^2 / diag(CONT%*%PLM.EFF.Median.CL$Sigma%*%t(CONT))
qchisq(0.95,1)

############################################
# Assumption Check: Balance: Use GLMM
############################################

No.Adj.T <- No.Adj.T2 <- No.Adj.T3 <- matrix(0,dim(X)[2],TNF)
Adj.T <- Adj.T2 <- Adj.T3 <- matrix(0,dim(X)[2],TNF)
fff <- function(x){sample(x,1)}

randmat <- data.frame(r=1:N,CID=TrueData$GP)
MM <- aggregate(randmat$r~randmat$CID,FUN="length")[,2]

for(jj in 1:TNF){
  No.Adj.X <- A*X/mean(A) - (1-A)*X/mean(1-A)
  Adj.X <- A*X/NData[[jj]]$PS - (1-A)*X/(1-NData[[jj]]$PS)
  
  for(kk in 1:dim(X)[2]){
    No.Adj.T3[kk,jj] <- summary( lme4::lmer(No.Adj.X[,kk]~(1|TrueData$GP)) )$coefficients[3]
    Adj.T3[kk,jj] <- summary( lme4::lmer(Adj.X[,kk]~(1|TrueData$GP)) )$coefficients[3]
  }
  print(jj)
}
# write.csv(No.Adj.T3,"Balance_No_Adj_X.csv",row.names=F)
# write.csv(   Adj.T3,"Balance_Adj_X.csv",row.names=F)

T.stat <- cbind(colnames(X)[1:7]," & ",
                round(apply(No.Adj.T3,1,median),2)[1:7]," & ",
                round(apply(Adj.T3,1,median),2)[1:7]," & " ,
                colnames(X)[7+1:7]," & ",
                round(apply(No.Adj.T3,1,median),2)[7+1:7]," & ",
                round(apply(Adj.T3,1,median),2)[7+1:7]," \\\\ \\hline " )

print(data.frame(T.stat),row.names=FALSE)

############################################
# Assumption Check: Overlap
############################################

par(mar=c(3,3,1,1))
H1 <- hist( NF.median$PS[A==1],xlim=c(0.4,1),xlab="",ylim=c(0,10),ylab="",prob=TRUE,
            breaks=seq(0,1,0.01),border=FALSE,col=rgb(1,0,0,0.2),main="")
title(xlab="Estimated PS",ylab="Density",line=2)
par(new=T)
H2 <- hist( NF.median$PS[A==0],xlim=c(0.4,1),xlab="",ylim=c(0,10),ylab="",prob=TRUE,xaxt='n',yaxt='n',
            breaks=seq(0,1,0.01),border=FALSE,col=rgb(0,0,1,0.2),main="")
legend("topright",legend=c("A=1","A=0"),pch=c(15,15),col=c(rgb(1,0,0,0.2),rgb(0,0,1,0.2)))

text(0.45,15-6,sprintf("Range of the PS under treatment (A=1)\n= [%0.4f,%0.4f]",min(NF.median$PS[A==1]),max(NF.median$PS[A==1])),pos=4)
text(0.45,13-6,sprintf("Range of the PS under control (A=0)\n= [%0.4f,%0.4f]",  min(NF.median$PS[A==0]),max(NF.median$PS[A==0])),pos=4)

############################################
# Comparison to GRF
############################################

GRF <- grf::causal_forest(X,Y,A,num.trees = 1000,clusters = GP) # change 1000 to a larger number

C1 <- 1*TrueData$C_location_city +  2*TrueData$C_location_rural + (1-TrueData$C_location_city-TrueData$C_location_rural)*3
C1.Name <- c("Center City","Rural","Urban")
Gname <- apply(expand.grid(C1.Name),1,function(v){paste(v,collapse="\n")})
CC <- (C1)

G <- length(unique(CC))
M <- matrix(0,N,G)

for(bb in 1:N){
  M[bb,CC[bb]] <- 1
}

RGRF <- sapply(1:3,function(tt){grf::average_treatment_effect(GRF,subset=M[,tt]==1)})
colnames(RGRF) <- Gname

# write.csv(RGRF,
#           "GRF_3.csv",row.names = F)

Comb.Est <- rbind(as.matrix(LC$EST),matrix(RGRF[1,],1,3))
Comb.SE  <- rbind(as.matrix(LC$Sigma),matrix(RGRF[2,],1,3))
Comb.Est <- t(cbind(Comb.Est[,1],Comb.Est[,3],Comb.Est[,2]))
Comb.SE  <- t(cbind(Comb.SE[,1],Comb.SE[,3],Comb.SE[,2]))
Comb.Est.Print <- Comb.Est
Comb.SE.Print  <- Comb.SE
for(ii in 1:3){
  for(jj in 1:4){
    Comb.Est.Print[ii,jj] <- sprintf("%0.3f",Comb.Est[ii,jj])
    Comb.SE.Print[ii,jj]  <- sprintf("%0.3f",Comb.SE[ii,jj])
  }
}

CI <- function(est,se,alpha=0.05,G=1){
  sprintf("(%0.3f,%0.3f)",
          est-qnorm(1-(1-(1-0.05)^(1/G))/2)*se,
          est+qnorm(1-(1-(1-0.05)^(1/G))/2)*se)
}

Block <- function(jj){
  rbind(Comb.Est.Print[jj,],
        Comb.SE.Print[jj,],
        sapply(1:4,function(tt){CI(Comb.Est[jj,tt],Comb.SE[jj,tt],alpha=0.05,G=1)}),
        sapply(1:4,function(tt){CI(Comb.Est[jj,tt],Comb.SE[jj,tt],alpha=0.05,G=3)}) )
}

F.Result <- data.frame( cbind(c("\\multirow{4}{*}{Central City}","","","",
                                "\\multirow{4}{*}{Urban}","","","",
                                "\\multirow{4}{*}{Rural}","","",""),
                              rep(c("Estimate","SE","95% CI","95% SCI"),3),
                              rbind( Block(1),
                                     Block(2),
                                     Block(3))) )

FF.Result <- matrix(0,12,12)
for(bb in 1:5){
  FF.Result[,2*bb-1] <- F.Result[,bb]
  FF.Result[,2*bb] <- " & "
}
bb <- 6
FF.Result[,2*bb-1] <- F.Result[,bb]
FF.Result[,2*bb] <- rep(c("\\\\ \\cline{2-6}","\\\\ \\cline{2-6}","\\\\ \\cline{2-6}","\\\\ \\hline"),3)


print(data.frame(FF.Result),row.names=F)
