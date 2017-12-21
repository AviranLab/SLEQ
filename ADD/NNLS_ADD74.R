type <- "SEQ"
RNAname <- "ADD74"
p <- scan(paste("/Users/huali/Desktop/Codes/ADD/",RNAname,type,"para.txt",sep=""), what = list(0));
p <- unlist(p);
r <- scan(paste("/Users/huali/Desktop/Codes/ADD/",RNAname,type,"Response.txt",sep=""), what = list(0));
r <- unlist(r);
x <- scan(paste("/Users/huali/Desktop/Codes/ADD/",RNAname,type,"DSMatrix.txt",sep=""), what = list(0));
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
idxt
rhot