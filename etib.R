source("functions.R") 

#------------------------------------------------------------------
# ET-IB functions
#------------------------------------------------------------------

#cumulative distribution function of ET-IB distribution with parameters alpha=6.3, beta=36.2, lambda=-0.9 and phi=0.4
detib(y=0.3,alpha=6.3,beta=36.2,lambda=-0.9,phi=0.4)

#probability density function of ET-IB distribution with parameters alpha=6.3, beta=36.2, lambda=-0.9 and phi=0.4
petib(q=0.5,alpha=6.3,beta=36.2,lambda=-0.9,phi=0.4)

#Quantile function of ET-IB distribution with parameters alpha=6.3, beta=36.2, lambda=-0.9 and phi=0.4, for u=1/2 
qetib(u=1/2,alpha=6.3,beta=36.2,lambda=-0.9,phi=0.4)

#Generation of n=1000 pseudo-random numbers of ET-IB distribution with parameters alpha=6.3, beta=36.2, lambda=-0.9 and phi=0.4
pseudo_ETIB<-retib(n=1000,alpha=6.3,beta=36.2,lambda=-0.9,phi=0.4)

#ETIB fit in pseudo_ET-IB data
etibfit(pseudo_ETIB)
