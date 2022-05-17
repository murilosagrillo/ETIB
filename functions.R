library(extraDistr)
library(AdequacyModel)
library(GenSA)
library(extraDistr)

#------------------------------------------------------------------
#Exponentiated transmuted-inverted beta (ET-IB) distribution functions#
#------------------------------------------------------------------

#------------------------------------------------------------------
# PDF ET-IB
#------------------------------------------------------------------

detib <- Vectorize(function(y,alpha,beta,lambda,phi,log = FALSE){
  logden <- log(phi *(((1 + lambda) *(y/( 1 + y))^(-1 + alpha) *(-(y/(1 + y)^2) + 1/(1 + y))* (1 - y/( 1 + y))^(-1 + beta))/(beta (alpha, beta)) - ( 2* lambda *(y/( 1 + y))^(-1 + alpha) *(-(y/(1 + y)^2) + 1/(1 + y)) *(1 - y/( 1 + y))^(-1 + beta) *(pbetapr (y, alpha, beta)))/(beta (alpha, beta)))* ((1 + lambda) *(pbetapr (y, alpha, beta)) - lambda* (pbetapr(y, alpha, beta))^2)^(-1 + phi))
  val <- ifelse(log, logden, exp(logden)) 
  return(val)
})

# detib(0.3,0.3,0.4,0.4,0.2)

#------------------------------------------------------------------
# Quantile function ETIB
#------------------------------------------------------------------

qetib <- Vectorize(function(u,alpha,beta,lambda,phi){
  val <- qbetapr((1+lambda-sqrt(((1+lambda)^2)-4*u^(1/phi)*lambda))/(2*lambda),alpha,beta)
  return(val)
})

qetib(0.5,1.5,1.3,0.2,0.4)

#------------------------------------------------------------------
# CDF ET-IB
#------------------------------------------------------------------

petib <- Vectorize(function(q,alpha,beta,lambda,phi,log.p = FALSE){
  cdf <-   ((1 + lambda)* ((pbetapr (q, alpha, beta))) - lambda*((pbetapr (q, alpha, beta)))^2)^phi
  val <- ifelse(log.p, log(cdf), cdf)
  return(val)
})

#petib(0.06833306,1.5,1.3,0.2,0.4)

#------------------------------------------------------------------
# Random number generation ET-IB
#------------------------------------------------------------------

retib <- function(n,alpha,beta,lambda,phi){
  u <- runif(n)
  val <- qbetapr((1+lambda-sqrt(((1+lambda)^2)-4*u^(1/phi)*lambda))/(2*lambda),alpha,beta)
  return(val)
}

#retib(100,0.3,0.9,0.3,0.4)

#------------------------------------------------------------------
# ETIB Estimation ETIB
#------------------------------------------------------------------

"etibfit"<-function(y){
  
  y <- as.vector(y)
  n<-length(y)
  set.seed(137)
  
  fdpETIB<-function(par){
    alpha=par[1]
    beta=par[2] 
    lambda=par[3]
    phi=par[4]
    phi *(((1 + lambda) *(y/( 1 + y))^(-1 + alpha) *(-(y/(1 + y)^2) + 1/(1 + y))* (1 - y/( 1 + y))^(-1 + beta))/(beta (alpha, beta)) - ( 2* lambda *(y/( 1 + y))^(-1 + alpha) *(-(y/(1 + y)^2) + 1/(1 + y)) *(1 - y/( 1 + y))^(-1 + beta) *(pbetapr (y, alpha, beta)))/(beta (alpha, beta)))* ((1 + lambda) *(pbetapr (y, alpha, beta)) - lambda* (pbetapr(y, alpha, beta))^2)^(-1 + phi)
  }
  
  #loglike_theta
  l_theta<- function(par){
    sum(log(fdpETIB(par)))
  }
  
  
  ##INICIO genSA para obter os valores de chutes
  fit.sa <- function(y,fdpETIB) {
    minusllike <- function(y) -sum(log(fdpETIB(c(y[1],y[2],y[3],y[4]))))
    lower <- c(0.001,0.001,-0.99,0.001) #may need some changes here
    upper <- c(250,250,0.99,250)
    out <- GenSA(lower = lower, upper = upper,
                 fn = minusllike, control=list(verbose=TRUE,max.time=2))
    return(out[c("value","par","counts")])
  }
  
  chute<-fit.sa(y,fdpETIB)$par
  ##fim do genSA
  
  res<-optim(chute, l_theta, method="BFGS",#x=x,
             control=list(fnscale=-1, pgtol=1e-20, maxit=200))#, silent=T)
  # if(class(res)=="try-error" || res$conv != 0 ||res$par[2]>.999||res$par[2]<.07) # a classe dos objetos que cont?m o erro,
  
  alpha_emv<-res$par[1]
  beta_emv<-res$par[2]
  lambda_emv<-res$par[3]
  phi_emv<-res$par[4]
  
  fdpETIB_densidade<-function(y){
    phi_emv *(((1 + lambda_emv) *(y/( 1 + y))^(-1 + alpha_emv) *(-(y/(1 + y)^2) + 1/(1 + y))* (1 - y/( 1 + y))^(-1 + beta_emv))/(beta (alpha_emv, beta_emv)) - ( 2* lambda_emv *(y/( 1 + y))^(-1 + alpha_emv) *(-(y/(1 + y)^2) + 1/(1 + y)) *(1 - y/( 1 + y))^(-1 + beta_emv) *(pbetapr (y, alpha_emv, beta_emv)))/(beta (alpha_emv, beta_emv)))* ((1 + lambda_emv) *(pbetapr (y, alpha_emv, beta_emv)) - lambda_emv* (pbetapr(y, alpha_emv, beta_emv))^2)^(-1 + phi_emv)
      }
  
  p<-4
  AIC<-(-2*sum(log( phi_emv *(((1 + lambda_emv) *(y/( 1 + y))^(-1 + alpha_emv) *(-(y/(1 + y)^2) + 1/(1 + y))* (1 - y/( 1 + y))^(-1 + beta_emv))/(beta (alpha_emv, beta_emv)) - ( 2* lambda_emv *(y/( 1 + y))^(-1 + alpha_emv) *(-(y/(1 + y)^2) + 1/(1 + y)) *(1 - y/( 1 + y))^(-1 + beta_emv) *(pbetapr (y, alpha_emv, beta_emv)))/(beta (alpha_emv, beta_emv)))* ((1 + lambda_emv) *(pbetapr(y, alpha_emv, beta_emv)) - lambda_emv* (pbetapr(y, alpha_emv, beta_emv))^2)^(-1 + phi_emv)    ))+2*p)
  BIC<-(-2*sum(log(  phi_emv *(((1 + lambda_emv) *(y/( 1 + y))^(-1 + alpha_emv) *(-(y/(1 + y)^2) + 1/(1 + y))* (1 - y/( 1 + y))^(-1 + beta_emv))/(beta (alpha_emv, beta_emv)) - ( 2* lambda_emv *(y/( 1 + y))^(-1 + alpha_emv) *(-(y/(1 + y)^2) + 1/(1 + y)) *(1 - y/( 1 + y))^(-1 + beta_emv) *(pbetapr (y, alpha_emv, beta_emv)))/(beta (alpha_emv, beta_emv)))* ((1 + lambda_emv) *(pbetapr (y, alpha_emv, beta_emv)) - lambda_emv* (pbetapr(y, alpha_emv, beta_emv))^2)^(-1 +phi_emv)   ))+p*log(length(y)))

  emvs<-c(alpha_emv,beta_emv,lambda_emv,phi_emv)
  dados_ordenados = sort(y)
  cdfETIB<-function(emvs,dados_ordenados){
    ((1 + lambda_emv)* ((pbetapr (dados_ordenados, alpha_emv, beta_emv))) - lambda_emv*((pbetapr (dados_ordenados, alpha_emv, beta_emv)))^2)^phi_emv
  }
  
  cdfETIB_ks_test<-function(y,par){
    ((1 + lambda_emv)* ((pbetapr (dados_ordenados, alpha_emv, beta_emv))) - lambda_emv*((pbetapr (dados_ordenados, alpha_emv, beta_emv)))^2)^phi_emv
  }
  
  v = cdfETIB(emvs, dados_ordenados)
  s = qnorm(v)
  u = pnorm((s - mean(s))/sqrt(var(s)))
  W_temp <- vector()
 # A_temp <- vector()
  for (i in 1:n) {
    W_temp[i] = (u[i] - (2 * i - 1)/(2 * n))^2
#    A_temp[i] = (2 * i - 1) * log(u[i]) + (2 * n + 1 - 2 * i) * log(1 - u[i])
  }
#  A_2 = -n - mean(A_temp)
  W_2 = sum(W_temp) + 1/(12 * n)
  W_star = W_2 * (1 + 0.5/n)
#  A_star = A_2 * (1 + 0.75/n + 2.25/n^2)
  KS = ks.test(x = y, y="cdfETIB_ks_test", par = emvs)
  
  g<-c()
  estim <- emvs
  g$estimativas <- estim
  g$convergencia <- res$conv
  # g$valor<-res$value
  data<-y
  hist(data,nclass = 10, freq = F,xlim=c(0.001,1),ylim=c(0.001,7))
  # plot(density(y),xlim=c(0.001,max(y)))
  curve(fdpETIB_densidade,xlim=c(0.001,1),ylim=c(0.001,.5),lty=3,add=T)
  # g$chute<-chute
  g$AIC<-AIC
  g$BIC<-BIC
  g$W_star<- W_star
#  g$A_star<- A_star
  g$KS<- KS
  return(g)
}
