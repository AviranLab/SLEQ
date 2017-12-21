type <- "MAP"
RNAname <- "Allele"
#FOLDER = "/Users/huali/Desktop/Codes/Fluoride_riboswitch/"
FOLDER = "/Users/huali/Desktop/Codes/Human_riboSNitch_alleles/DMSMaPSeq/"
#FOLDER = "/Users/huali/Dropbox/New_Results_Review1/"
#FOLDER = "/Users/huali/Dropbox/"
#p <- scan(paste(FOLDER,RNAname,type,"para.txt",sep=""), what = list(0));
#p <- unlist(p);
r <- scan(paste(FOLDER,RNAname,type,"Response.txt",sep=""), what = list(0));
r <- unlist(r);
x <- scan(paste(FOLDER,RNAname,type,"DSMatrix.txt",sep=""), what = list(0));
x <- unlist(x); 
m<-length(r)
n<-length(x)/m
idx<-NULL
FoldsNum <- p[2];
for(i in 3:length(p))
  idx<-c(idx,p[i]);

library("nnls")
X <- t(array(x,c(n,m)))

dup_idx = which(duplicated(t(X))==1)
X = X[,-dup_idx]

dup_idx_inverse = 1:n
dup_idx_inverse = dup_idx_inverse[-dup_idx]

fit_results <- nnls(X, r)
para <- coef(fit_results)
SStot <- sum((r-mean(r))^2)
SSres <- sum((r-X%*%para)^2)
R2 <- 1- SSres/SStot
R2

rho <- para/sum(para)
id <- which(rho>0.01)
rhot <- rho[id]
rhot <- rhot/sum(rhot)
#idxt <- idx[id]
idxt <- dup_idx_inverse[id]


#select all duplicate structures of the most abundant selected structure A 
#sidx = c(3,6,8,10,13,15,18,24,25,29,30)
#StructuresA_select_group = StructuresA_tol[,15:47]
#StructuresA_select_group = StructuresA_select_group[,sidx]
#select_s = StructuresA[idxt[1],15:47]
#select_s = select_s[sidx]

#mark = c()
#for(i in 1:dim(StructuresA_select_group)[1])
#{
#  if(sum(StructuresA_select_group[i,]==select_s)==11)
#  {mark = c(mark,i)}
#}
#StructuresDotA_tol[mark]

#StructuresC_select_group = StructuresC_tol[,8:40]
#StructuresC_select_group = StructuresC_select_group[,sidx]
#select_s = StructuresC[idxt[5]-92,8:40]
#select_s = select_s[sidx]

#mark = c()
#for(i in 1:dim(StructuresC_select_group)[1])
#{
#  if(sum(StructuresC_select_group[i,]==select_s)==11)
#  {mark = c(mark,i)}
#}
#StructuresDotC_tol[mark]

idxt
rhot
#StructuresDot[idxt]
dim(X)
rhot[length(rhot)]+rhot[length(rhot)-1]

