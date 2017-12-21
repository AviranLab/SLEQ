r <- scan("/Users/huali/Desktop/Codes/Human_riboSNitch_alleles/DMSMaPSeq/AlleleMAPResponse_shortfolds.txt", what = list(0));
r <- unlist(r);
x <- scan("/Users/huali/Desktop/Codes/Human_riboSNitch_alleles/DMSMaPSeq/AlleleMAPDSMatrix_shortfolds.txt", what = list(0));
x <- unlist(x);
m<-length(r)
n<-length(x)/m
library("nnls")
X <- t(array(x,c(n,m)))
para <- coef(nnls(X, r))
rho <- para/sum(para)
id <- which(rho>0.01)
rhot <- rho[id]
rhot <- rhot/sum(rhot)
