############################################
# Last update : Dec 19 2021
# Code for Section 4
############################################

PARA <- expand.grid(c(0,1),c(0,75),c(0,1,2,3))

FOLDER <- "Result"
FN <- "SIM_GRF"



F1 <- F2 <- list()

for(iter in 1:dim(PARA)[1]){
    
    G <- 4
    
    PST <- if(PARA[iter,1]==0){"CPS"} else {"VPS"}
    
    
    CN <- colnames(read.csv(sprintf("%s/Result_%s_COR_Est_rho%0.3d_B%0.4d_SubB%0.4d.csv",
                                    FOLDER,PST,PARA[iter,2],PARA[iter,3],1)))
    F1[[iter]] <- F2[[iter]] <- data.frame(matrix(0,500,80))
    colnames(F1[[iter]]) <- colnames(F2[[iter]]) <- CN
    
    for(jj in 1:100){
        file1 <- sprintf("%s/Result_%s_COR_Est_rho%0.3d_B%0.4d_SubB%0.4d.csv",
                         FOLDER,PST,PARA[iter,2],PARA[iter,3],jj)
        file2 <- sprintf("%s/Result_%s_COR_CL_rho%0.3d_B%0.4d_SubB%0.4d.csv",
                         FOLDER,PST,PARA[iter,2],PARA[iter,3],jj)
        F1[[iter]][(jj-1)*5+1:5,] <- read.csv(file1)
        F2[[iter]][(jj-1)*5+1:5,] <- read.csv(file2)
        
        
    }
    
    write.csv(F1[[iter]],sprintf("%s/Result_%s_COR_Est_rho%0.3d_B%0.4d.csv",
                                 "Summary_Result",PST,PARA[iter,2],PARA[iter,3]),
              row.names=FALSE)
    
    write.csv(F2[[iter]],sprintf("%s/Result_%s_COR_CL_rho%0.3d_B%0.4d.csv",
                                 "Summary_Result",PST,PARA[iter,2],PARA[iter,3]),
              row.names=FALSE)
    
}


for(iter in 1:dim(PARA)[1]){
    
    G <- 4
    
    PST <- if(PARA[iter,1]==0){"CPS"} else {"VPS"}
    
    
    CN <- colnames(read.csv(sprintf("%s/Result_%s_COR_Est_rho%0.3d_B%0.4d_SubB%0.4d.csv",
                                    FOLDER,PST,PARA[iter,2],PARA[iter,3],1)))
    
    F1[[iter]] <- read.csv(sprintf("%s/Result_%s_COR_Est_rho%0.3d_B%0.4d.csv",
                                   "Summary_Result",PST,PARA[iter,2],PARA[iter,3]))
    
    F2[[iter]] <- read.csv(sprintf("%s/Result_%s_COR_CL_rho%0.3d_B%0.4d.csv",
                                   "Summary_Result",PST,PARA[iter,2],PARA[iter,3]))
    
}









PLOT5 <- function(GP=1){
    
    COL <- c(rgb(0,0,0,0.5),4,2,"orange") ; LTY <- c(1,1,1,1) ; RNC <- TRUE ; PCH <- c(15,17,3,4) ; CEX <- c(1.25,1.25,1.25,1.25)
    
    LL <- rbind( cbind(c(30,31,31,31,32,32,32,33), 
                       c(34,35,36,37,38,39,40,41),
                       rbind(26:29, rbind( matrix(1:24,6,4), 25))),
                 c(42,42,43:46) )
    
    
    LLf <- rbind( LL[1:4,], 47 ,LL[5:9,] ) 
    
    layout( LLf,
            widths=c(0.5,0.5,3,3,3,3),
            heights=c(1, 3,3,3, 0.5, 3,3,3, 0.5,1))
    
    MAR <- c(2,2,0.25,0.25)
    
    par(mar=MAR)
    
    Ps2 <- function(FFFF,pos,type=NULL){
        
        
        if(type=="B"){
            CCpos <- c(which(colnames(FFFF[[1]])==sprintf("Bias_GP_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("Bias_DR_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("Bias_Lin01_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("Bias_GRF_%0.1d",GP)))
            
            CCpos2 <- c(which(colnames(FFFF[[1]])==sprintf("SE_GP_%0.1d",GP)),
                        which(colnames(FFFF[[1]])==sprintf("SE_DR_%0.1d",GP)),
                        which(colnames(FFFF[[1]])==sprintf("SE_Lin01_%0.1d",GP)),
                        which(colnames(FFFF[[1]])==sprintf("SE_GRF_%0.1d",GP)))
            
            CC <- list()
            
            for(gg in 1:4){
                CC[[gg]] <- cbind( c((FFFF[[pos[gg]]][,CCpos[1]]/FFFF[[pos[gg]]][,CCpos2[1]]),
                                     (FFFF[[pos[gg]]][,CCpos[2]]/FFFF[[pos[gg]]][,CCpos2[2]]),
                                     (FFFF[[pos[gg]]][,CCpos[3]]/FFFF[[pos[gg]]][,CCpos2[3]]),
                                     (FFFF[[pos[gg]]][,CCpos[4]]/FFFF[[pos[gg]]][,CCpos2[4]])),
                                   rep(1:4+(gg-1)*4,each=dim(FFFF[[pos[[gg]]]])[1]) )
            }
            
            CC.bind <- rbind(CC[[1]],CC[[2]],CC[[3]],CC[[4]])
            
        } else if (type=="V"){
            CCpos <- c(which(colnames(FFFF[[1]])==sprintf("SE_GP_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("SE_DR_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("SE_Lin01_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("SE_GRF_%0.1d",GP)))
            
            CC <- list()
            
            for(gg in 1:4){
                CC[[gg]] <- cbind( c((FFFF[[pos[gg]]][,CCpos[1]]/FFFF[[pos[gg]]][,CCpos[3]]),
                                     (FFFF[[pos[gg]]][,CCpos[2]]/FFFF[[pos[gg]]][,CCpos[3]]),
                                     (FFFF[[pos[gg]]][,CCpos[3]]/FFFF[[pos[gg]]][,CCpos[3]]),
                                     (FFFF[[pos[gg]]][,CCpos[4]]/FFFF[[pos[gg]]][,CCpos[3]])),
                                   rep(1:4+(gg-1)*4,each=dim(FFFF[[pos[[gg]]]])[1]) )
            }
            
            CC.bind <- rbind(CC[[1]],CC[[2]],CC[[3]],CC[[4]])
            
        } else if (type=="C"){
            CCpos <- c(which(colnames(FFFF[[1]])==sprintf("Cover_GP_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("Cover_DR_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("Cover_Lin01_%0.1d",GP)),
                       which(colnames(FFFF[[1]])==sprintf("Cover_GRF_%0.1d",GP)))
            
            CC <- list()
            
            for(gg in 1:4){
                CC[[gg]] <- cbind( c(mean(FFFF[[pos[gg]]][,CCpos[1]]),
                                     mean(FFFF[[pos[gg]]][,CCpos[2]]),
                                     mean(FFFF[[pos[gg]]][,CCpos[3]]),
                                     mean(FFFF[[pos[gg]]][,CCpos[4]])),
                                   c(1:4+(gg-1)*4) )
            }
            
            CC.bind <- rbind(CC[[1]],CC[[2]],CC[[3]],CC[[4]])
            
        }
        
        XX <- c(c(1,2,3,4),
                c(1,2,3,4)+8,
                c(1,2,3,4)+8*2,
                c(1,2,3,4)+8*3)
        Xt <- apply(matrix(XX,4,4),2,mean)
        Xr <- max(XX)+1
        XAXIS.size <- 0.5
        YAXIS.size <- 0.5
        Xline <- Yline <- 0.4
        TICK.L <- -0.04
        
        if(type=="B"){ 
            CC.agg <- aggregate(CC.bind[,1]~CC.bind[,2],FUN="mean")[,c(2,1)]
            
            IND1 <- ifelse(CC.agg[,1]>  1.005,TRUE,FALSE)
            IND2 <- ifelse(CC.agg[,1]< -1.005,TRUE,FALSE)
            
            Ygap <- ifelse( CC.agg[,1] > 1.005, 1.005, CC.agg[,1] )
            Ygap <- ifelse( Ygap < -1.005, -1.005, Ygap )
            
            plot( XX[c(0,4,8,12)+1], Ygap[c(0,4,8,12)+1], xlab="",ylab="",
                  xlim=c(0,Xr),xaxt='n',yaxt='n',col=COL[1], pch=PCH[1],
                  ylim=c(-1.00,1.00), cex=CEX[1])
            
            axis(1,at=Xt,label=NA,tck=TICK.L)
            mtext(c(0,1,2,3),side=1,line=Xline,at=Xt,cex=XAXIS.size)
            
            
            Yt <- c(-1,0,1)
            Ylab <- c(expression(""<="-1"),Yt[2],expression("">="1"))
            axis(2,at=Yt[1:3],label=NA,tck=TICK.L)
            axis(2,at=Yt,label=NA,tck=TICK.L)
            mtext(Ylab,side=2,line=Yline,at=Yt,cex=YAXIS.size)
            
            points( XX[c(0,4,8,12)+2], Ygap[c(0,4,8,12)+2], col=COL[2], pch=PCH[2], cex=CEX[2])
            points( XX[c(0,4,8,12)+3], Ygap[c(0,4,8,12)+3], col=COL[3], pch=PCH[3], cex=CEX[3])
            points( XX[c(0,4,8,12)+4], Ygap[c(0,4,8,12)+4], col=COL[4], pch=PCH[4], cex=CEX[4])
            
            
            if(sum(IND1)>0){
                text( XX[which(IND1)], Ygap[which(IND1)]-0.2, 
                      sprintf("%0.2f",CC.agg[which(IND1),1]), col=rep(COL,4)[which(IND1)], cex=0.8)
            }
            
            if(sum(IND2)>0){
                text( XX[which(IND2)], Ygap[which(IND2)]+0.2, 
                      sprintf("%0.2f",CC.agg[which(IND2),1]), col=rep(COL,4)[which(IND2)], cex=0.8)
            }
            
            
            abline(h=0.0,col=3,lty=3)
            
            
        } else if (type=="V"){
            
            CC.agg <- aggregate(CC.bind[,1]~CC.bind[,2],FUN="mean")[,c(2,1)]
            
            
            IND1 <- ifelse(CC.agg[,1]>1.21,TRUE,FALSE)
            IND2 <- ifelse(CC.agg[,1]<0.89,TRUE,FALSE)
            
            
            Ygap <- ifelse( CC.agg[,1] > 1.21, 1.21, CC.agg[,1] )
            Ygap <- ifelse( Ygap < 0.89, 0.89, Ygap )
            
            plot( XX[c(0,4,8,12)+1], Ygap[c(0,4,8,12)+1], xlab="",ylab="",
                  xlim=c(0,Xr),xaxt='n',yaxt='n',col=COL[1], pch=PCH[1],
                  ylim=range(0.9,1.22), cex=CEX[1])
            
            axis(1,at=Xt,label=NA,tck=TICK.L)
            mtext(c(0,1,2,3),side=1,line=Xline,at=Xt,cex=XAXIS.size)
            
            Yt <- c(0.9, 1, 1.1, 1.2)
            Ylab <- c(expression(""<="0.9"),Yt[2:3],expression("">="1.2"))
            axis(2,at=Yt[1:4],label=NA,tck=TICK.L)
            axis(2,at=Yt,label=NA,tck=TICK.L)
            mtext(Ylab,side=2,line=Yline,at=Yt,cex=YAXIS.size)
            
            points( XX[c(0,4,8,12)+2], Ygap[c(0,4,8,12)+2], col=COL[2], pch=PCH[2], cex=CEX[2])
            points( XX[c(0,4,8,12)+4], Ygap[c(0,4,8,12)+4], col=COL[4], pch=PCH[4], cex=CEX[4])
            
            
            if(sum(IND1)>0){
                text( XX[which(IND1)], Ygap[which(IND1)]-0.04, 
                      sprintf("%0.2f",CC.agg[which(IND1),1]), col=rep(COL,4)[which(IND1)], cex=0.8)
            }
            
            if(sum(IND2)>0){
                text( XX[which(IND2)], Ygap[which(IND2)]+0.04, 
                      sprintf("%0.2f",CC.agg[which(IND2),1]), col=rep(COL,4)[which(IND2)], cex=0.8)
            }
            
            
            
            abline(h=1,col=3,lty=3)
            
            
            
            # XXX <- seq(-2,Xr+2,length=101)
            # YYY <- sin(XXX)*0.01+1.20
            # lines(XXX,YYY-0.005,col=rgb(0,0,0,0.25))
            # lines(XXX,YYY+0.005,col=rgb(0,0,0,0.25))
            
            
        } else if (type=="C"){
            CC.agg <- aggregate(CC.bind[,1]~CC.bind[,2],FUN="mean")[,c(2,1)]
            
            IND1 <- ifelse(CC.agg[,1]<  0.85,TRUE,FALSE)
            
            Ygap <- ifelse( CC.agg[,1]< 0.85, 0.85, CC.agg[,1] )
            
            plot( XX[c(0,4,8,12)+1], Ygap[c(0,4,8,12)+1], xlab="",ylab="",
                  xlim=c(0,Xr),xaxt='n',yaxt='n',col=COL[1], pch=PCH[1],
                  ylim=c(0.85,1), cex=CEX[1])
            
            axis(1,at=Xt,label=NA,tck=TICK.L)
            mtext(c(0,1,2,3),side=1,line=Xline,at=Xt,cex=XAXIS.size)
            
            
            Yt <- c(0.85,0.9,0.95,1)
            Ylab <- c(expression(""<="0.85"),Yt[2:4])
            axis(2,at=Yt[1:4],label=NA,tck=TICK.L)
            axis(2,at=Yt,label=NA,tck=TICK.L)
            mtext(Ylab,side=2,line=Yline,at=Yt,cex=YAXIS.size)
            
            points( XX[c(0,4,8,12)+2], Ygap[c(0,4,8,12)+2], col=COL[2], pch=PCH[2], cex=CEX[2])
            points( XX[c(0,4,8,12)+3], Ygap[c(0,4,8,12)+3], col=COL[3], pch=PCH[3], cex=CEX[3])
            points( XX[c(0,4,8,12)+4], Ygap[c(0,4,8,12)+4], col=COL[4], pch=PCH[4], cex=CEX[4])
            
            
            if(sum(IND1)>0){
                text( XX[which(IND1)], Ygap[which(IND1)]+0.02, 
                      sprintf("%0.2f",CC.agg[which(IND1),1]), col=rep(COL,4)[which(IND1)], cex=0.8)
            }
            
            
            abline(h=0.95,col=3,lty=3)
            
        }
        
    }
    
    for(pst in 0:1){
        for(rhot in unique(PARA[,2])){
            
            pos <- which( PARA[,1]==pst & PARA[,2]==rhot )
            
            Ps2(F1,pos,type="B")
            Ps2(F1,pos,type="V")
            Ps2(F1,pos,type="C")
            
            Ps2(F2,pos,type="B")
            Ps2(F2,pos,type="V")
            Ps2(F2,pos,type="C")
            
            print(c(pst,rhot))
        }
    }
    
    par(mar=c(0,MAR[2],0,MAR[4]))
    plot.new() # 25 x-axis
    text(0.5,0.5,expression(beta),cex=1.2)
    
    plot.new() # 26 head 1
    text(0.5,0.5,expression("IID, CPS"),cex=1.15)
    
    plot.new() # 27 head 2
    text(0.5,0.5,expression("CL, CPS"),cex=1.15)
    
    plot.new() # 28 head 3
    text(0.5,0.5,expression("IID, VPS"),cex=1.15)
    
    plot.new() # 29 head 4
    text(0.5,0.5,expression("CL, VPS"),cex=1.15)
    
    par(mar=c(0,0,0,0))
    
    plot.new() # 30 empty
    
    par(mar=c(MAR[1],0,MAR[3],0))
    
    plot.new() # 31 left panel
    text(0.5,0.5,"Results Using Non-Cluster-Robust\nStandard Errors",cex=1.2,srt=90)
    
    plot.new() # 32 right panel
    text(0.5,0.5,"Results Using Cluster-Robust\nStandard Errors",cex=1.2,srt=90)
    
    par(mar=c(0,0,0,0))
    
    plot.new() # 33 empty
    plot.new() # 34 empty
    for(jj in 1:2){
        par(mar=c(MAR[1],0,MAR[3],0))
        plot.new() # 35,38 RBI
        text(0.3,0.5,"Bias",cex=0.95,srt=90)
        text(0.7,0.5,"(in z-score units)",cex=0.8,srt=90)
        plot.new() # 36,39 RBI
        text(0.5,0.5,expression("SEs"),cex=0.95,srt=90)
        plot.new() # 37,40 Coverage
        text(0.5,0.5,"Coverages",cex=0.95,srt=90)
        
    }
    
    par(mar=c(0,0,0,0))
    
    plot.new() # 41 empty
    
    plot.new() # 42 empty
    
    par(mar=c(0,MAR[2],0,MAR[4]))
    for(jj in 1:4){
        plot.new() # 43-46 legend
        # polygon(c(1,4,4,1,1)/10,c(8,8,6,6,8)/10,col=COL[jj])
        # segments(0.25,0.8,0.25,0.6,lwd=1.5)
        points(0.25,0.5,cex=2,pch=PCH[jj],col=COL[jj])
        if(jj==1){ text(0.5,0.5,expression(hat(tau)["SSLS"]),   pos=4,cex=1.2) } 
        if(jj==2){ text(0.5,0.5,expression(hat(tau)["EIF"]) ,   pos=4,cex=1.2) }
        if(jj==3){ text(0.5,0.5,expression(hat(tau)["W"]),      pos=4,cex=1.2) }
        if(jj==4){ text(0.5,0.5,expression(hat(tau)["GRF"]),    pos=4,cex=1.2) }
        
    }
}




png("plot/Simulation_GP1.png",height=7,width=7,unit="in",res=500)
PLOT5(GP=1)
dev.off()

png("plot/Simulation_GP2.png",height=7,width=7,unit="in",res=500)
PLOT5(GP=2)
dev.off()

png("plot/Simulation_GP3.png",height=7,width=7,unit="in",res=500)
PLOT5(GP=3)
dev.off()

png("plot/Simulation_GP4.png",height=7,width=7,unit="in",res=500)
PLOT5(GP=4)
dev.off()
