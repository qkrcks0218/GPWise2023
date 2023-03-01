####################################################
# Source
####################################################

source("MySL.R")
source("SSLS.R")
library(readstata13)

####################################################
# Cleaning 
# Download data from https://nces.ed.gov/ecls/dataproducts.asp
# Make ECLSK_Kto8_child_STATA.dta
####################################################

# RAWDATA <- read.dta13("ECLSK_Kto8_child_STATA.dta") # <- This takes too long time; recommend to run the following stata .do file
# set maxvar 30000
# use ECLSK_Kto8_child_STATA.dta 
# keep CHILDID PARENTID S1_ID S2_ID S4_ID T1_ID T2_ID T4_ID CREGION KURBAN_R GENDER R1_KAGE C1CMOTOR P1CENTER P1HFAMIL WKPARED WKRACETH WKSESL S2KSCTYP C1R4MTSC C1R4RTSC C2R4MTSC C2R4RTSC C4R4MTSC C4R4RTSC 
# save ${path}/short_ECLSK.dta, replace

RAWDATA <- read.dta13("short_ECLSK.dta") 
ShortData <- RAWDATA[,c("CHILDID","PARENTID","S1_ID","S2_ID","S4_ID","T1_ID","T2_ID","T4_ID",
                        "CREGION","KURBAN_R","GENDER","R1_KAGE","C1CMOTOR","P1CENTER","P1HFAMIL",
                        "WKPARED","WKRACETH","WKSESL","S2KSCTYP",
                        "C1R4RTSC","C1R4MTSC","C1RGTSCO","C1CMOTOR")]

Data <- ShortData[!is.na(ShortData$CREGION) & !is.na(ShortData$KURBAN_R) & !is.na(ShortData$S2KSCTYP) ,]
Data <- Data[!is.na(Data$GENDER) & 
               ( !is.na(Data$WKRACETH) & Data$WKRACETH!="NOT ASCERTAINED" & Data$WKRACETH!="NOT APPLICABLE" ) & 
               !is.na(Data$R1_KAGE) & (Data$R1_KAGE>0) & 
               !is.na(Data$P1HFAMIL) & 
               ( !is.na(Data$C1CMOTOR) & Data$C1CMOTOR>0 ) & !is.na(Data$WKSESL) & !is.na(Data$WKPARED) &
               ( !is.na(Data$P1CENTER) & Data$P1CENTER!="NOT ASCERTAINED"& Data$P1CENTER!="NOT APPLICABLE" ) ,]

Data <- Data[ ( !is.na(Data$C1R4RTSC) & Data$C1R4RTSC>0 ) , ]                ## READING score (1st semester)
Data <- Data[ ( !is.na(Data$C1R4MTSC) & Data$C1R4MTSC>0 ) , ]                ## MATH score (1st semester)
Data <- Data[ ( !is.na(Data$C1RGTSCO) & Data$C1RGTSCO>0 ) , ]                ## General score (1st semester)
Data <- Data[ ( !is.na(Data$C1CMOTOR) & Data$C1CMOTOR>0 ) , ]                ## Motor Skill (1st semester)
Data <- Data[ Data$S1_ID!="", ]

My.Data <- Data[,c("CHILDID","PARENTID","S1_ID","S2_ID","S4_ID","T1_ID","T2_ID","T4_ID",
                   "CREGION","KURBAN_R","S2KSCTYP","GENDER","R1_KAGE","WKRACETH",
                   "P1HFAMIL","WKPARED","WKSESL",
                   "P1CENTER",
                   "C1R4RTSC","C1R4MTSC","C1RGTSCO","C1CMOTOR")]
colnames(My.Data) <- c("Child.ID","Parent.ID","Kind.ID1","Kind.ID2","Kind.ID4","Teacher.ID1","Teacher.ID2","Teacher.ID4",
                       "C_region","C_location","C_public","S_sex","S_age","S_race",
                       "S_familytype","S_parentaledu","S_sestatus",
                       "A",
                       "Y_Read","Y_Math","Y_GenKnow","Y_Motor")


My.Data$C_region_NE <- as.numeric( My.Data$C_region=="NORTHEAST" )
My.Data$C_region_SOUTH <- as.numeric( My.Data$C_region=="SOUTH" )
My.Data$C_region_WEST <- as.numeric( My.Data$C_region=="WEST" )

My.Data$C_location_city <- as.numeric( My.Data$C_location=="CENTRAL CITY" )
My.Data$C_location_rural <- as.numeric( My.Data$C_location=="SMALL TOWN AND RURAL" )

My.Data$C_public_public <- as.numeric( My.Data$C_public=="PUBLIC" )

My.Data$C_M <- 0
for(jj in 1:dim(My.Data)[1]){
  My.Data$C_M[jj] <- sum(My.Data$Kind.ID1[jj]==My.Data$Kind.ID1)
}

My.Data <- My.Data[My.Data$C_M<=25,]

My.Data$S_sex_male <- as.numeric( My.Data$S_sex=="MALE" )

My.Data$S_race_H <- as.numeric( My.Data$S_race=="HISPANIC, RACE SPECIFIED" | My.Data$S_race=="HISPANIC, RACE NOT SPECIFIED" )
My.Data$S_race_B <- as.numeric( My.Data$S_race=="BLACK OR AFRICAN AMERICAN, NON-HISPANIC" )
My.Data$S_race_A <- as.numeric( My.Data$S_race=="ASIAN" )
My.Data$S_race_W <- as.numeric( My.Data$S_race=="WHITE, NON-HISPANIC" )

My.Data$S_familytype_bothparent <- as.numeric( My.Data$S_familytype=="2 PARENTS PLUS SIBLINGS" | My.Data$S_familytype=="2 PARENTS NO SIBLING" )

My.Data$S_parentaledu_college <- as.numeric( My.Data$S_parentaledu=="DOCTORATE OR PROFESSIONAL DEGREE" | 
                                               My.Data$S_parentaledu=="BACHELOR'S DEGREE" |
                                               My.Data$S_parentaledu=="MASTER'S DEGREE (MA, MS)" |
                                               My.Data$S_parentaledu=="SOME COLLEGE")
My.Data$S_age <- as.numeric(My.Data$S_age)
My.Data$S_sestatus <- as.numeric(My.Data$S_sestatus)
My.Data$A <- as.numeric(My.Data$A=="YES")
My.Data$GP <- My.Data$Kind.ID1

indmat <- cbind(1:length(unique( My.Data$GP )),sort(unique( My.Data$GP )))
for(ii in 1:dim(My.Data)[1]){
  My.Data$GP[ii] <- indmat[which(indmat[,2]==My.Data$GP[ii]),1]
}
My.Data$GP <- as.numeric(My.Data$GP)

Final.Data <- My.Data[,c("Y_Read","Y_Math","Y_GenKnow","Y_Motor",
                         "A", 
                         "C_M",
                         "C_region_NE","C_region_SOUTH","C_region_WEST",
                         "C_location_city","C_location_rural",
                         "C_public_public",
                         "S_sex_male","S_age","S_race_A","S_race_B","S_race_H","S_race_W",
                         "S_familytype_bothparent","S_parentaledu_college","S_sestatus",
                         "GP")]
write.csv(Final.Data,"ECLSK.csv",row.names=FALSE)

####################################################
# Read Data
####################################################

TrueData <- read.csv("ECLSK.csv")
TrueData$S_age <- scale(TrueData$S_age)
TrueData$S_sestatus <- scale(TrueData$S_sestatus)

TrueData0 <- TrueData
TrueData0[,"A"] <- 0
TrueData1 <- TrueData
TrueData1[,"A"] <- 1

N <- dim(TrueData)[1]

####################################################
# Split Sample for Cross Fitting
####################################################

SSI <- SSC <- matrix(0,N,100)
for(bb in 1:100){
  set.seed(bb)
  GPU <- unique(TrueData$GP)
  MSC <- sort(sample(GPU,length(GPU)/2))
  ASC <- (GPU)[-MSC]
  for(tt in MSC){
    pos <- which(TrueData$GP==tt)
    SSC[pos,bb] <- 1
    RS <- as.numeric(length(pos)>3)*3 + as.numeric(length(pos)<=3)*(length(pos))
    SSI[sample(pos,RS),bb] <- 1
  }
  for(tt in ASC){
    pos <- which(TrueData$GP==tt)
    SSC[pos,bb] <- 2
    RS <- as.numeric(length(pos)>3)*3 + as.numeric(length(pos)<=3)*(length(pos))
    SSI[sample(pos,RS),bb] <- 2
  }
}

write.csv(SSC,"SSCL_Index.csv",row.names=F)
write.csv(SSI,"SSInd_Index.csv",row.names=F)

####################################################
# Analysis: 
# Recommend to use cluster computing system by
# iterating BATCH in 1:100
####################################################

SL.hpara <- list()
SL.hpara$SLL <- c(1,2,3,4,5,6,7,9,10)
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
SL.hpara$MTRY <- c(5,10)                  # random forest parameters
SL.hpara$MLPL <- c(2,4)                   # number of nodes in MLP
SL.hpara$NMN <- 50                        # gbm parameter
SL.hpara$MLPdecay <- 10^c(-4,-5)          # MLP decay parameter

pos.X  <- 3:16
pos.AX <- 2:16

SSI <- list()
SSI$MS <- which(SSI.Index[,BATCH]==1)
SSI$AS <- which(SSI.Index[,BATCH]==2)

SSC <- list()
SSC$MS <- which(SSC.Index[,BATCH]==1)
SSC$AS <- which(SSC.Index[,BATCH]==2)

####################################################
# Machine learning estimators
####################################################

Y.Read.ML.Fit <- Y.Read.X.ML.Fit <- list()
A.ML.Fit <- list()

Y.Read.ML.Fit$MS <- MySL( Data=TrueData[SSI$MS,],
                          locY=1,
                          locX=pos.AX,
                          Ydist=gaussian(),
                          SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN)
Y.Read.ML.Fit$AS <- MySL( Data=TrueData[SSI$AS,],
                          locY=1,
                          locX=pos.AX,
                          Ydist=gaussian(),
                          SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN)

Y.Read.X.ML.Fit$MS <- MySL( Data=TrueData[SSI$MS,],
                          locY=1,
                          locX=pos.X,
                          Ydist=gaussian(),
                          SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN)
Y.Read.X.ML.Fit$AS <- MySL( Data=TrueData[SSI$AS,],
                          locY=1,
                          locX=pos.X,
                          Ydist=gaussian(),
                          SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN)

A.ML.Fit$MS <- MySL( Data=TrueData[SSI$MS,],
                     locY=2,
                     locX=pos.X,
                     Ydist=binomial(),
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN)
A.ML.Fit$AS <- MySL( Data=TrueData[SSI$AS,],
                     locY=2,
                     locX=pos.X,
                     Ydist=binomial(),
                     SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPdecay=SL.hpara$MLPdecay, NMN=SL.hpara$NMN)

####################################################
# Obtain Estimated Nuisance Functions
####################################################

Predict.ML.Mat <- matrix(0,N,6)
colnames(Predict.ML.Mat) <- c("SSIndex",
                              "OR1","OR0","OR","ORX","PS")

Y.Read.Predict <- Y.Read.X.Predict <- list()
A.Predict <- list()

Y.Read.Predict$AS <- predict(Y.Read.ML.Fit$MS,newdata=TrueData[SSC$AS,pos.AX],onlySL=TRUE)$pred
Y.Read.Predict$MS <- predict(Y.Read.ML.Fit$AS,newdata=TrueData[SSC$MS,pos.AX],onlySL=TRUE)$pred

Y.Read.Predict$AS1 <- predict(Y.Read.ML.Fit$MS,newdata=TrueData1[SSC$AS,pos.AX],onlySL=TRUE)$pred
Y.Read.Predict$MS1 <- predict(Y.Read.ML.Fit$AS,newdata=TrueData1[SSC$MS,pos.AX],onlySL=TRUE)$pred

Y.Read.Predict$AS0 <- predict(Y.Read.ML.Fit$MS,newdata=TrueData0[SSC$AS,pos.AX],onlySL=TRUE)$pred
Y.Read.Predict$MS0 <- predict(Y.Read.ML.Fit$AS,newdata=TrueData0[SSC$MS,pos.AX],onlySL=TRUE)$pred

Y.Read.X.Predict$AS0 <- predict(Y.Read.X.ML.Fit$MS,newdata=TrueData0[SSC$AS,pos.X],onlySL=TRUE)$pred
Y.Read.X.Predict$MS0 <- predict(Y.Read.X.ML.Fit$AS,newdata=TrueData0[SSC$MS,pos.X],onlySL=TRUE)$pred

A.Predict$AS <- predict(A.ML.Fit$MS,newdata=TrueData[SSC$AS,pos.X],onlySL=TRUE)$pred
A.Predict$MS <- predict(A.ML.Fit$AS,newdata=TrueData[SSC$MS,pos.X],onlySL=TRUE)$pred

Predict.ML.Mat[SSC$MS,] <- cbind(1,
                                 Y.Read.Predict$MS1,
                                 Y.Read.Predict$MS0,
                                 Y.Read.Predict$MS,
                                 Y.Read.X.Predict$MS,
                                 A.Predict$MS)
Predict.ML.Mat[SSC$AS,] <- cbind(2,
                                 Y.Read.Predict$AS1,
                                 Y.Read.Predict$AS0,
                                 Y.Read.Predict$AS,
                                 Y.Read.X.Predict$AS,
                                 A.Predict$AS)

Predict.ML.Mat <- data.frame(Predict.ML.Mat)

####################################################
# Save Results
####################################################

write.csv(Predict.ML.Mat,
          sprintf("NData/Result_BATCH%0.4d.csv",BATCH),
          row.names=F)

