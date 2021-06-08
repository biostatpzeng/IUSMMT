Calculate_Pvalue <- function(M,G,X,OS){
  ####################################################################################################################################
  #########################################                   Function body                ###########################################
  ####################################################################################################################################
 ####### INPUT
  ####### M : High-dimensional exposures that can be either data.frame or matrix. Rows represent samples, columns represent variables
  ####### G : Independent variable that is a vector
  ####### X : covariables that can be either data.frame or matrix
  ####### OS : Survival data. Rows represents samples, the first column is the survival time, the second column is the survival status
  ####################################################################################################################################
  ####### Values 
  ####### P_gamma_TE : The P value of total effect.
  ####### P_gamma_DE : The P value of direct effect.
  ####### sigma_alpha : the coefficient can reflect the association of M –> G, note that the effect is adjusted by covariables when covariables
  #######	are not NULL.
  ####### P_alpha : The P value of alpha.
  ####### beta : the coefficient can reflect the association of G –> OS, note that the effect is adjusted by X. When covariables are not
  #######	NULL, the effect is adjusted by X and covariables. 
  ####### P_beta : The P value of beta.
  ####################################################################################################################################
  pkgs <- list("data.table","coxKM","SKAT", "coxme","EMMREML")
  checking<-unlist(lapply(pkgs, require, character.only = T))
  if(any(checking==F))
    stop("Please install the necessary packages first!")
  M <- as.matrix(data.frame(apply(t(na.omit(t(M))),2,scale)))
  G <- as.matrix(data.frame(apply(t(na.omit(t(G))),2,scale)))
  X <- as.matrix(data.frame(apply(t(na.omit(t(X))),2,scale)))
  OS <- data.frame(OS)
  colnames(OS) <- c('time','event')
  OS[,1] <- OS[,1]+0.001
  data1 = list(M = M,G=G,x=X,time = OS$time,event=OS$event)
  ## step 1   total effect and DE -----coxKM
  Gamma = coxph(Surv(OS$time,OS$event) ~ X)$coef
  P_gamma_TE = coxKM(Z=M, U=OS$time, Delta=OS$event, X=X, gamma=Gamma, kernel="linear")$`p.value`
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
  
  Re <- c(P_gamma_TE = p_gamma_TE,
          P_gamma_DE = p_gamma_DE,
          sigma_alpha = sigma_alpha,
          P_alpha = p_alpha,
          beta = beta,
          P_beta = p_beta)
  return(Re)
}

IUSMMT <- function(p_a, p_b){
  ####### INPUT
  ####### p_a : The P value of alpha.
  ####### p_b : The P value of beta.
  ####################################################################################################################################
  ####### Values 
  ####### p_IUSMMT : joint significant test for mediators.
  ####################################################################################################################################
  input_pvalues <- cbind(p_a, p_b)
  input_pvalues <- apply(input_pvalues,2,dup.fun)
  nullprop <- null_estimation(input_pvalues,lambda=0.5)
  p_IUSMMT <- fdr_est(nullprop$alpha00,nullprop$alpha01,nullprop$alpha10,nullprop$alpha1,
                    nullprop$alpha2,input_pvalues,exact=1)
  return(p_IUSMMT)
}

extract_coxme_table <- function (mod){
  beta <- fixef(mod)
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 10)
  p<- signif(1 - pchisq((beta/se)^2, 1), 10)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}

null_estimation <- function (input_pvalues, lambda = 0.5) 
{
  if (is.null(ncol(input_pvalues))) 
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2) 
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues)) 
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues), 
  ]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) < 
      1) 
    stop("input_pvalues doesn't have valid p-values")
  pcut <- seq(0.1, 0.8, 0.1)
  frac1 <- rep(0, 8)
  frac2 <- rep(0, 8)
  frac12 <- rep(0, 8)
  for (i in 1:8) {
    frac1[i] <- mean(input_pvalues[, 1] >= pcut[i])/(1 - 
                                                       pcut[i])
    frac2[i] <- mean(input_pvalues[, 2] >= pcut[i])/(1 - 
                                                       pcut[i])
    frac12[i] <- mean(input_pvalues[, 2] >= pcut[i] & input_pvalues[, 
                                                                    1] >= pcut[i])/(1 - pcut[i])^2
  }
  alpha00 <- min(frac12[pcut == lambda], 1)
  if (ks.test(input_pvalues[, 1], "punif", 0, 1, alternative = "greater")$p > 
      0.05) 
    alpha1 <- 1
  else alpha1 <- min(frac1[pcut == lambda], 1)
  if (ks.test(input_pvalues[, 2], "punif", 0, 1, alternative = "greater")$p > 
      0.05) 
    alpha2 <- 1
  else alpha2 <- min(frac2[pcut == lambda], 1)
  if (alpha00 == 1) {
    alpha01 <- 0
    alpha10 <- 0
    alpha11 <- 0
  }
  else {
    if (alpha1 == 1 & alpha2 == 1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
      alpha00 <- 1
    }
    if (alpha1 == 1 & alpha2 != 1) {
      alpha10 <- 0
      alpha11 <- 0
      alpha01 <- alpha1 - alpha00
      alpha01 <- max(0, alpha01)
      alpha00 <- 1 - alpha01
    }
    if (alpha1 != 1 & alpha2 == 1) {
      alpha01 <- 0
      alpha11 <- 0
      alpha10 <- alpha2 - alpha00
      alpha10 <- max(0, alpha10)
      alpha00 <- 1 - alpha10
    }
    if (alpha1 != 1 & alpha2 != 1) {
      alpha10 <- alpha2 - alpha00
      alpha10 <- max(0, alpha10)
      alpha01 <- alpha1 - alpha00
      alpha01 <- max(0, alpha01)
      if ((1 - alpha00 - alpha01 - alpha10) < 0) {
        alpha11 <- 0
        alpha10 <- 1 - alpha1
        alpha01 <- 1 - alpha2
        alpha00 <- 1 - alpha10 - alpha01
      }
      else {
        alpha11 <- 1 - alpha00 - alpha01 - alpha10
      }
    }
  }
  alpha.null <- list(alpha10 = alpha10, alpha01 = alpha01, 
                     alpha00 = alpha00, alpha1 = alpha1, alpha2 = alpha2)
  return(alpha.null)
}

fdr_est <- function (alpha00, alpha01, alpha10, alpha1, alpha2, input_pvalues, exact = 0) 
{
  if (is.null(ncol(input_pvalues))) 
    stop("input_pvalues should be a matrix or data frame")
  if (ncol(input_pvalues) != 2) 
    stop("inpute_pvalues should have 2 column")
  input_pvalues <- matrix(as.numeric(input_pvalues), nrow = nrow(input_pvalues))
  if (sum(complete.cases(input_pvalues)) < nrow(input_pvalues)) 
    warning("input_pvalues contains NAs to be removed from analysis")
  input_pvalues <- input_pvalues[complete.cases(input_pvalues), 
  ]
  if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues) < 
      1) 
    stop("input_pvalues doesn't have valid p-values")
  pmax <- apply(input_pvalues, 1, max)
  nmed <- length(pmax)
  efdr1 <- rep(0, nmed)
  if (exact == 0) {
    for (i in 1:nmed) {
      fdr11 <- (pmax[i] * alpha01)/mean(pmax <= pmax[i])
      fdr12 <- (pmax[i] * alpha10)/mean(pmax <= pmax[i])
      fdr2 <- (pmax[i] * pmax[i] * alpha00)/mean(pmax <= 
                                                   pmax[i])
      efdr1[i] <- fdr11 + fdr12 + fdr2
    }
  }
  if (exact == 1) {
    nmed <- nrow(input_pvalues)
    cdf12 <- input_pvalues
    xx1 <- c(0, input_pvalues[order(input_pvalues[, 1]), 
                              1])
    yy1 <- c(0, seq(1, nmed, by = 1)/nmed)
    gfit1 <- gcmlcm(xx1, yy1, type = "lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots) * gfit1$slope.knots)
    xx2 <- c(0, input_pvalues[order(input_pvalues[, 2]), 
                              2])
    yy2 <- c(0, seq(1, nmed, by = 1)/nmed)
    gfit2 <- gcmlcm(xx2, yy2, type = "lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots) * gfit2$slope.knots)
    if (alpha1 != 1) 
      Fknots1 <- (Fknots1 - alpha1 * xknots1)/(1 - alpha1)
    else Fknots1 <- rep(0, length(xknots1))
    if (alpha2 != 1) 
      Fknots2 <- (Fknots2 - alpha2 * xknots2)/(1 - alpha2)
    else Fknots2 <- rep(0, length(xknots2))
    orderq1 <- pmax
    orderq2 <- pmax
    gcdf1 <- pmax
    gcdf2 <- pmax
    for (i in 1:length(xknots1)) {
      if (i == 1) {
        gcdf1[orderq1 <= xknots1[i]] <- (Fknots1[i]/xknots1[i]) * 
          orderq1[orderq1 <= xknots1[i]]
      }
      else {
        if (sum(orderq1 > xknots1[i - 1] & orderq1 <= 
                xknots1[i]) > 0) {
          temp <- orderq1[orderq1 > xknots1[i - 1] & 
                            orderq1 <= xknots1[i]]
          gcdf1[orderq1 > xknots1[i - 1] & orderq1 <= 
                  xknots1[i]] <- Fknots1[i - 1] + (Fknots1[i] - 
                                                     Fknots1[i - 1])/(xknots1[i] - xknots1[i - 
                                                                                             1]) * (temp - xknots1[i - 1])
        }
      }
    }
    for (i in 1:length(xknots2)) {
      if (i == 1) {
        gcdf2[orderq2 <= xknots2[i]] <- (Fknots2[i]/xknots2[i]) * 
          orderq2[orderq2 <= xknots2[i]]
      }
      else {
        if (sum(orderq2 > xknots2[i - 1] & orderq2 <= 
                xknots2[i]) > 0) {
          temp <- orderq2[orderq2 > xknots2[i - 1] & 
                            orderq2 <= xknots2[i]]
          gcdf2[orderq2 > xknots2[i - 1] & orderq2 <= 
                  xknots2[i]] <- Fknots2[i - 1] + (Fknots2[i] - 
                                                     Fknots2[i - 1])/(xknots2[i] - xknots2[i - 
                                                                                             1]) * (temp - xknots2[i - 1])
        }
      }
    }
    gcdf1 <- ifelse(gcdf1 > 1, 1, gcdf1)
    gcdf2 <- ifelse(gcdf2 > 1, 1, gcdf2)
    cdf12[, 1] <- gcdf1
    cdf12[, 2] <- gcdf2
    for (i in 1:nmed) {
      fdr11 <- (pmax[i] * cdf12[i, 2] * alpha01)/mean(pmax <= 
                                                        pmax[i])
      fdr12 <- (pmax[i] * cdf12[i, 1] * alpha10)/mean(pmax <= 
                                                        pmax[i])
      fdr2 <- (pmax[i] * pmax[i] * alpha00)/mean(pmax <= 
                                                   pmax[i])
      efdr1[i] <- fdr11 + fdr12 + fdr2
    }
  }
  efdr1.order <- efdr1[order(pmax, decreasing = T)]
  for (i in 2:nmed) {
    efdr1.order[i] <- min(efdr1.order[i], efdr1.order[i - 
                                                        1])
  }
  efdr1 <- efdr1.order[rank(-pmax)]
  return(efdr = efdr1)
}
