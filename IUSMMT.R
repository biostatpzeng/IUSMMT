IUSMMT <- function(M,G,X,OS){
  M <- as.matrix(data.frame(apply(t(na.omit(t(M))),2,scale)))
  G <- as.matrix(data.frame(apply(t(na.omit(t(G))),2,scale)))
  X <- as.matrix(data.frame(apply(t(na.omit(t(X))),2,scale)))
  OS <- data.frame(OS)
  colnames(OS) <- c('time','event')
  OS[,1] <- OS[,1]+0.001
  data1 = list(M = M,G=G,x=X,time = OS$time,event=OS$event)
  ## step 1   total effect and DE -----coxKM
  Gamma = coxph(Surv(OS$time,OS$event) ~ X)$coef
  p_gamma_total = coxKM(Z=M, U=OS$time, Delta=OS$event, X=X, gamma=Gamma, kernel="linear")$`p.value`
  Gamma = coxph(Surv(OS$time,OS$event) ~ cbind(X,G))$coef
  p_gamma_DE = coxKM(Z=M, U=OS$time, Delta=OS$event, X=cbind(X,G), gamma=Gamma, kernel="linear")$`p.value`
  ## step 2   alpha  -----SKAT
  K1 = M%*%t(M)/length(G)
  fit1 = emmreml(y = G,X = X,Z=diag(length(G)),K = K1)
  sigma_alpha = sum(fit1$uhat)
  obj = SKAT_Null_Model(G ~ X, out_type="C")
  p_alpha = SKAT(M,obj)$`p.value`
  ## step 3   beta  -----coxme
  if(dim(M)[1]<=dim(M)[2]){
    K=M%*%t(M)
    eig=eigen(K)
    evals<-eig$values
    index_evals = which(evals<0); if(length(index_evals)>0){evals[index_evals] <- 0}
    evecs<-eig$vectors
    K12=as.matrix(eig$vectors%*%diag(sqrt(evals))%*%t(eig$vectors))
    fit2 <- coxme(Surv(OS$time,OS$event) ~ G + X + (K12 | 1))
  } else{
    K12 = M
    fit2 <- coxme(Surv(OS$time,OS$event) ~ G + X + (K12 | 1))
  }
  beta = fixef(fit2)[1]
  p_beta = extract_coxme_table(fit2)$p[1]
  beta = fixef(fit2)[1]
  p_IUT = max(p_beta,p_alpha)

  Re <- c(P_gamma_total = p_gamma_total,
          P_gamma_DE = p_gamma_DE,
          sigma_alpha = sigma_alpha,
          P_alpha = p_alpha,
          beta = beta,
          P_beta = p_beta,
          P_IUT = p_IUT)
  return(Re)
}
