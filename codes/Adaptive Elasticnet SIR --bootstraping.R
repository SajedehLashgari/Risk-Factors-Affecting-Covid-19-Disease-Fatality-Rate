### by Sajedeh Lashgari
# Adaptive Elastic-net Sliced Inverse Regression
###################################################

# Addressing
#setwd("")

# Setting showing 
options(scipen=20)

# Requirment libraries
## for ELR and RLR 
require(glmnet)
library(tictoc)

# Get data
df = read.csv("final dataset.csv")

# About data
df = df[,2:62]

x_ = as.matrix(df[,1:60])
y_ = df[,61]

########--------------------------  Alg  -------------------------########
##########################################################################
########### Create Penalized Sliced Inverse Regression Alg ###############
##########################################################################

PSIR = function (X, Y, H, n, p, data, r=1, choosing.d, alpha, no.dim=0, adaptive){
  if (no.dim != 0) 
    choosing.d="given"
  ORD = order(Y)
  X = X[ORD, ]
  Y = Y[ORD]
  ms = array(0, n)
  m = floor(n/H)
  c = n%%H
  M = matrix(0, nrow=H, ncol=n)
  if (c == 0) {
    M = diag(H) %x% matrix(1, nrow=1, ncol=m)/m
    ms = m+ms
  }
  else {
    for (i in 1:c) {
      M[i, ((m+1)*(i-1)+1):((m+1)*i)] = 1/(m+1)
      ms[((m+1)*(i-1)+1):((m+1)*i)] = m
    }
    for (i in (c+1):H) {
      M[i, ((m+1)*c+(i-c-1)*m+1):((m+1)*c+(i-c)*m)] = 1/m
      ms[((m+1)*c+(i-c-1)*m+1):((m+1)*c+(i-c)*m)] = m-1
    }
  }
  keep.ind = c(1:p)
  X = X[, keep.ind]
  X.H = matrix(0, nrow=H, ncol=dim(X)[2])
  grand.mean = matrix(apply(X, 2, mean), nrow=1, ncol=dim(X)[2])
  X.stand.ord = X 
  X.H = M %*% X.stand.ord
  svd.XH = svd(X.H, nv=p)
  res.eigen.value = array(0, p)
  res.eigen.value[1:dim(X.H)[1]] = (svd.XH$d)^2/H
  ###########
  # Kernel
  ########### 
  ## Kernel LR
  if (adaptive == "TRUE") {
    lr = lm(IFR~. ,data=data)
    betaK_lr = as.data.frame(lr$coefficients)[, 1]
    weightt = as.matrix(1/abs(betaK_lr[2:dim(data)[2]])^(r))
    weightt[weightt[,1] == Inf] = 999999999 ## Replacing values estimated as Infinite for 999999999
    weightt[is.na(weightt[,1])] = 0
  }
  else {
    weightt = c(1)
  }
  if (choosing.d == "manual") {
    plot(c(1:p), res.eigen.value, ylab="eigen values")
    no.dim = (function() {
      dim = readline("Choose the number of directions:   ")
      dim = as.numeric(unlist(strsplit(dim, ",")))
      return(dim)
    })()
  }
  if (choosing.d == "automatic") {
    beta.hat = array(0, c(p, min(p, H)))
    Y.tilde = array(0, c(n, min(p, H)))
    for (ii in 1:min(p, H)) {
      eii = matrix(0, nrow=dim(svd.XH$v)[2], ncol=1)
      eii[ii] = 1
      eigen.vec = solve(t(svd.XH$v), eii)
      Y.tilde[, ii] = t(M) %*% M %*% X.stand.ord %*% eigen.vec/(res.eigen.value[ii])*
        matrix(1/ms, nrow=n, ncol=1)
    }
    mus = array(0, min(p, H))
    for (ii in 1:min(p, H)) {
      lars.fit.cv = cv.glmnet(X.stand.ord, Y.tilde[, ii], nfolds=10, alpha=alpha)
      ## choose the one with the smallest cvm
      ind = max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind = 2
      lambda = lars.fit.cv$lambda[ind]
      mus[ii] = lambda
      lars.fit = glmnet(X.stand.ord, Y.tilde[, ii], lambda=lambda, penalty.factor=weightt, alpha=alpha)
      beta.hat[keep.ind, ii] = as.double(lars.fit$beta)
    }
    ## The statistic for determining d is ||beta_i||*lambda_i
    temp.2 = sqrt(apply(beta.hat^2, 2, sum))*res.eigen.value[1:H]
    temp = temp.2/temp.2[1]
    res.kmeans = kmeans(temp, centers=2)
    no.dim = min(sum(res.kmeans$cluster == 1), sum(res.kmeans$cluster == 2))
  }
  beta.hat = array(0, c(p, no.dim))
  Y.tilde = array(0, c(n, no.dim))
  for (ii in 1:no.dim) {
    eii = matrix(0, nrow=dim(t(svd.XH$v))[2], ncol=1)
    eii[ii] = 1
    eigen.vec = solve(t(svd.XH$v), eii)
    Y.tilde[, ii] = t(M) %*% M %*% X.stand.ord %*% eigen.vec/(res.eigen.value[ii])*
      matrix(1/ms, nrow=n, ncol=1)
  }
  mus = array(0, no.dim)
  if (no.dim == 1) {
    lars.fit.cv = cv.glmnet(X.stand.ord, Y.tilde, nfolds=10, alpha=alpha)
    
    ind = max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
    if (ind == 1) 
      ind=2
    lambda = lars.fit.cv$lambda[ind]
    lars.fit = glmnet(X.stand.ord, Y.tilde, lambda=lambda, penalty.factor=weightt, alpha=alpha)
    beta.hat[keep.ind] = as.double(lars.fit$beta)
  }
  else {
    for (ii in 1:no.dim) {
      lars.fit.cv = cv.glmnet(X.stand.ord, Y.tilde[, ii], nfolds=10, alpha=alpha)
      ind = max(which(lars.fit.cv$cvm == min(lars.fit.cv$cvm)))
      if (ind == 1) 
        ind = 2
      lambda = lars.fit.cv$lambda[ind]
      mus[ii] = lambda
      lars.fit = glmnet(X.stand.ord, Y.tilde[, ii], lambda=lambda, penalty.factor=weightt, alpha=alpha)
      beta.hat[keep.ind, ii] = as.double(lars.fit$beta)
    }
  }
  list(beta=beta.hat, eigen.value=res.eigen.value, eigen.vec=eigen.vec,
       no.dim=no.dim, Y.tilde=Y.tilde, lars.fit.cv=lars.fit.cv,
       lars.fit=lars.fit, X.stand.ord=X.stand.ord, weight=weightt, lambda=lambda)
}

##########################################################################
################# Create Penalized Linear Regression Alg #################
##########################################################################
pen_LR = function(X, Y, n, alpha, adaptive=FALSE, r, data){
  ###########
  # Kernel
  ###########
  if (adaptive == "TRUE") {
    lr = lm(IFR~. , data=data)
    betaK_lr = as.data.frame(lr$coefficients)[, 1]
    w = as.matrix(1/abs(betaK_lr[2:dim(data)[2]])^(r))
    w[w[,1] == Inf]=999999999 ## Replacing values estimated as Infinite for 999999999
    w[is.na(w[,1])]=0
  }
  else {
    w=c(1)
  }
  cv_P = cv.glmnet(X, Y, family ="gaussian", type.measure="mse", nfolds=10 , alpha=alpha)
  fit = glmnet(X, Y, family="gaussian", lambda=cv_P$lambda.min, alpha=alpha, penalty.factor=w)
  list(lambda=cv_P$lambda.min, fit=fit)
}

##########################################################################
########### find appropriate d
PSIR(x_, y_, H=12, n=dim(df)[1], p=dim(df)[2]-1,
     data=df, r=2, choosing.d='automatic', adaptive=TRUE, alpha=0.7)$no.dim

##########################################################################
############################## Bootstrapping##############################
#set.seed()

Bootstrapping = function(iter, alpha=1, A=TRUE, H=12, data, r=2, met='PSIR'){
  beta_boot = matrix(0,dim(data)[2]-1,iter)
  for (i in 1:iter){
    bootst = data[sample(nrow(data), nrow(data), TRUE),]
    boot_x = as.matrix(bootst[,1:dim(data)[2]-1])
    boot_y = as.matrix(bootst[,dim(data)[2]])
    
    if (met == 'PSIR') {
      beta_boot[,i] = PSIR(boot_x, boot_y, H=H, n=dim(bootst)[1], p=dim(bootst)[2]-1, 
                      data=bootst, r=r, choosing.d='given', adaptive=A, no.dim=1, alpha=alpha)$beta
    }
    if (met == 'PLR'){
      beta_boot[,i] = pen_LR(boot_x, boot_y, n=dim(bootst)[1], data=bootst, alpha=alpha, adaptive=A, r=2)$fit$beta[,1]
    }
    if (met == 'LR') {
      beta_boot = matrix(0,dim(data)[2],iter)
      for (i in 1:iter){
        bootst = data[sample(nrow(data), nrow(data), TRUE),]
        beta_boot[,i] = lm(IFR ~ ., bootst)$coefficients
      }
    } 
  }
  data.frame(beta_boot, sd=apply(beta_boot,1,sd))
}

##########################################################################
################################ Results #################################
########### Penalized Linear
{tic("LLR fitting")
LLR = Bootstrapping(iter=100, A=FALSE, data=df, met='PLR')
a = toc()
llr = a$toc - a$tic}

{tic("ALLR fitting")
ALLR = Bootstrapping(iter=100, data=df, met='PLR')
aA = toc()
allr = aA$toc - aA$tic}

{tic("ELR fitting")
ELR = Bootstrapping(iter=100, alpha=0.7, A=FALSE, data=df, met='PLR')
b = toc()
elr = b$toc - b$tic}

{tic("AELR fitting")
AELR = Bootstrapping(iter=100, alpha=0.7, data=df, met='PLR')
bB = toc()
aelr = bB$toc - bB$tic}

########### Penalized SIR
{tic("LSIR fitting")
LSIR = Bootstrapping(iter=100, A=FALSE, data=df)
c = toc()
lsir = c$toc - c$tic}

{tic("ALSIR fitting")
ALSIR = Bootstrapping(iter=100, data=df)
cC = toc()
alsir = cC$toc - cC$tic}

{tic("AESIR fitting")
AESIR = Bootstrapping(iter=100, alpha=0.7, data=df)
d = toc()
aesir = d$toc - d$tic}

tcpx_r2 = data.frame(AESIR=aesir, ALSIR=alsir, LSIR=lsir, AELR=aelr, ALLR=allr, ELR=elr, LLR=llr)
tcpx_r2

##########################################################################
m_cov = data.frame(AESIR=mean(AESIR$sd), ALSIR=mean(ALSIR$sd), LSIR=mean(LSIR$sd),
                AELR=mean(AELR$sd), ALLR=mean(ALLR$sd), ELR=mean(ELR$sd), LLR=mean(LLR$sd))
m_cov

sd_sd_cov = data.frame(AESIR=sd(AESIR$sd), ALSIR=sd(ALSIR$sd), LSIR=sd(LSIR$sd), 
                AELR=sd(AELR$sd), ALLR=sd(ALLR$sd), ELR=sd(ELR$sd), LLR=sd(LLR$sd))
round(cbind(sd_sd_cov[1:3]*(10^5), sd_sd_cov[4:7]),2)


sd_cov = data.frame(AESIR=AESIR$sd, ALSIR=ALSIR$sd, LSIR=LSIR$sd, AELR=AELR$sd, ALLR=ALLR$sd, 
                ELR=ELR$sd, LLR=LLR$sd)
sd_cov

AE_aml = PSIR(x_, y_, H=12, n=dim(df)[1], p=dim(df)[2]-1, data=df, r=2, 
              choosing.d='given', alpha=0.7, no.dim=1, adaptive=TRUE)
AL_aml = PSIR(x_, y_, H=12, n=dim(df)[1], p=dim(df)[2]-1, data=df, r=2, 
              choosing.d='given', alpha=1, no.dim=1, adaptive=TRUE)
L_aml = PSIR(x_, y_, H=12, n=dim(df)[1], p=dim(df)[2]-1, data=df, r=2, 
              choosing.d='given', alpha=1, no.dim=1, adaptive=FALSE)
l = pen_LR(x_, y_, alpha=1, adaptive=FALSE, r=2, data=df)
al = pen_LR(x_, y_, alpha=1, adaptive=TRUE, r=2, data=df)
e = pen_LR(x_, y_, alpha=0.7, adaptive=FALSE, r=2, data=df)
ae = pen_LR(x_, y_, alpha=0.7, adaptive=TRUE, r=2, data=df)

beta = cbind(data.frame('AESIR_beta'=AE_aml$beta), data.frame('ALSIR_beta'=AL_aml$beta), 
                data.frame('LSIR_beta'=AE_aml$beta), data.frame('AELR_beta'=ae$fit$beta[,1]),
                data.frame('ALLR_beta'=al$fit$beta[,1]), data.frame('ELR_beta'=e$fit$beta[,1]),
                data.frame('LLR_beta'=l$fit$beta[,1]))

##########################################################################
# save results
write.csv(sd_cov, 'sd_cov.csv',  row.names=FALSE)
write.csv(m_cov, 'mean_cov.csv',  row.names=FALSE)
write.csv(beta, 'beta_cov.csv',  row.names=FALSE)
