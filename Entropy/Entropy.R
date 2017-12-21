rm(list=ls())
type <- "SEQ"
RNAname <- "16S"
FOLDER = "/Users/huali/Dropbox/New_Results_Review1/"
#FOLDER="//Users/huali/Desktop/Codes/Fluoride_riboswitch/"
Oritext <- read.table(paste(FOLDER,RNAname,".txt",sep=""))#change
NumStruct <- 1000
OriStructs <- Oritext[2:(2+NumStruct-1),]
RNALen = nchar(as.character(OriStructs[1]))
Structures = matrix(nrow=NumStruct,ncol=RNALen)
StructuresDot = vector(length=NumStruct)
FoldCount = vector(mode="numeric",length=NumStruct)



p <- scan(paste("/Users/huali/Desktop/Codes/Fluoride_riboswitch/",RNAname,type,"para.txt",sep=""), what = list(0));
p <- unlist(p);
r <- scan(paste("/Users/huali/Desktop/Codes/Fluoride_riboswitch/",RNAname,type,"Response.txt",sep=""), what = list(0));
r <- unlist(r);
x <- scan(paste("/Users/huali/Desktop/Codes/Fluoride_riboswitch/",RNAname,type,"DSMatrix.txt",sep=""), what = list(0));
x <- unlist(x);
m<-length(r)
n<-length(x)/m
idx<-NULL
FoldsNum <- p[2];
for(i in 3:length(p))
  idx<-c(idx,p[i]);

library("nnls")
X <- t(array(x,c(n,m)))
para <- coef(nnls(X, r))
rho <- para/sum(para)
id <- which(rho>0.01)
rhot <- rho[id]
rhot <- rhot/sum(rhot)
idxt <- idx[id]



StructPK1<-c(".....1234..56789a4321......bcd....dcb....a98765..............................") #PK1
StructPK2<-c(".....12345.6789ab4321......cde....edc5...ba9876..............................") #PK1+LR1
StructPK3<-c(".....12345.6789ab4321......cde....edc5.f.ba9876f.............................") #PK+2LRI
StructPK4<-c(".....1234..56789a4321......bcd....dcb....a98765.efghijk....kjihgfe...........") #PK1+T
StructPK5<-c(".....12345.6789ab4321......cde....edc5...ba9876.fghijkl....lkjihgf...........") #PK1+LR1+T
StructPK6<-c(".....12345.6789ab4321......cde....edc5.f.ba9876fghijklm....mlkjihg...........") #PK1+2LRI+T
StructPK1<-unlist(strsplit(StructPK1,""))
StructPK2<-unlist(strsplit(StructPK2,""))
StructPK3<-unlist(strsplit(StructPK3,""))
StructPK4<-unlist(strsplit(StructPK4,""))
StructPK5<-unlist(strsplit(StructPK5,""))
StructPK6<-unlist(strsplit(StructPK6,""))
StructPK<-rbind(StructPK1,StructPK2,StructPK3,StructPK4,StructPK5,StructPK6)
StructPK<-StructPK[,1:RNALen]
 

TPP80 = c(1,2,2,2,2,2,2,2,2,1,2,2,2,2,1,2,2,2,1,1,1,1,1,3,3,3,3,3,3,1,1,1,1,1,1,3,1,1,3,3,3,3,1,1,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,1,1,1,1,3,3,3,3,1,1,1,3,3,3,1,1,3,3,3,3,1)
TPPDot80 =".((((((((..(((.(((.....)))))).........)))).....(((...((((......))))...)))..))))."

ADD74 = c(1,1,1,2,2,2,2,2,2,2,2,2,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,3,3,3,3,3,3,1,1,1,1,1,1,1,1,2,2,2,2,2,2,1,1,1,1,1,1,1,3,3,3,3,3,3,1,1,3,3,3,3,3,3,3,3,3)
ADDDot74 = "...(((((((((...((((((.........))))))........((((((.......))))))..)))))))))"

Allele_A = c(1,2,2,2,2,2,1,1,1,1,1,1,1,3,3,3,3,3,1,1,1,3,3,3,1,3,3,3,3,3,3,3,1);
Allele_C = c(1,1,1,2,2,2,2,2,1,2,2,2,2,1,2,2,1,1,1,1,1,1,1,1,3,3,1,3,3,3,3,1,3);
Allele_ADot = ".((((((((((....(((((.......)))))...))).)))))))."

t <- 1
flag <- 0
for (s in  1:NumStruct)
{
  CurrentStructDot <- OriStructs[s]
  CurrentStructDot <- as.character(CurrentStructDot)
  CurrentStructDotSplit <- unlist(strsplit(CurrentStructDot,""))
  CurrentStruct <- vector(mode="numeric",length=RNALen)
  for (i in 1: RNALen)
  {
    if(CurrentStructDotSplit[i]==".")
      CurrentStruct[i] <- 1
    else if(CurrentStructDotSplit[i]=="(")
      CurrentStruct[i] <- 2
    else
      CurrentStruct[i] <- 3
  }
  if(s != 1)
  {
    for (j in 1:(t-1))
    {
      if(sum(Structures[j,]==CurrentStruct)==RNALen)#others
        #if(sum(Structures[j,15:47]==CurrentStruct[15:47])==33)#Allele_A
        #if(sum(Structures[j,8:40]==CurrentStruct[8:40])==33)#Allele_C
      {
        FoldCount[j] = FoldCount[j]+1  
        flag <- 1			
      }
    }		
  }
  if(flag == 0)
  {
    Structures[t,] <- CurrentStruct
    StructuresDot[t] <- CurrentStructDot
    FoldCount[t] <- FoldCount[t] + 1
    t <- t+1 
  }
  flag <- 0
  
}
FoldNum <- t-1
Structures <- Structures[1:FoldNum,]
StructuresDot <- StructuresDot[1:FoldNum]
FoldCount <- FoldCount[1:FoldNum]



if(length(which(idxt>FoldNum))>0)
{
  idxn<-idxt[-which(idxt>FoldNum)]
  rhon<-rhot[-which(idxt>FoldNum)]
  idxPK<-idxt[which(idxt>FoldNum)]-FoldNum
  rhoPK<-rhot[which(idxt>FoldNum)]
}
if(length(which(idxt>FoldNum))==0)
{
  idxn<-idxt
  rhon<-rhot
  idxPK <- c()
  rhoPK <- c()
}



  Xbp <- matrix(0,nrow=RNALen,ncol=RNALen)
  for (i in 1:length(idxn))
  {
    leftposi <- vector(mode="integer",length=RNALen)
    t <- 1
    for(j in 1:RNALen)
    {
      if(Structures[idxn[i],j]==2)
      {
        leftposi[t] = j
        t = t+1
      }
      if(Structures[idxn[i],j]==3)
      {
        t = t- 1
        Xbp[j,leftposi[t]] = Xbp[j,leftposi[t]]+rhon[i]
        leftposi[t] = 0;
      }
    }
  }
  
  Xn <- matrix(0,nrow=RNALen,ncol=RNALen)
  for(i in 1:length(idxPK))
  {
    for(j in 1:ncol(StructPK))
      if(StructPK[idxPK[i],j]!=".")
      {
        k<-which(StructPK[idxPK[i],]==StructPK[idxPK[i],j])
        StructPK[idxPK[i],k]="."
        Xbp[k[2],j] = Xbp[k[2],j]+rhoPK[i]
        Xn[k[2],j] = Xn[k[2],j]+rhoPK[i]
      }
  }
  library("gplots")
  x11()
  x_hm <- heatmap.2(Xbp, Rowv=FALSE, Colv=FALSE, col=rev(heat.colors(256)),tracecol="white")
  
  
  flag <- 0
  pp<- vector(mode="integer",length=RNALen)
  for(i in 1:RNALen)
  {
    pp[i] = sum(Xbp[i,])
    if(pp[i]<1)
      Xbp[i,i] = 1- pp[i];
    if(pp[i] >1.00000001)
      flag<-1 
    
  }
  
  
  
  H <- vector(mode="numeric",length=RNALen)
  for(i in 1:RNALen)
  {
    if(sum(Xbp[i])>0)
      Xbp[i,] = Xbp[i,]/sum(Xbp[i,])
    H[i] <- 0
    for(j in 1:RNALen)
    {
      if(Xbp[i,j] != 0)
        H[i] = H[i] - Xbp[i,j]*log2(Xbp[i,j])
    }
  }
  idxt
  rhot
  flag
  mean(H)
 