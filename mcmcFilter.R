###irregular grid
library("invgamma", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.3")
#xdatagrid=cbind(xdata,seq(0,96,1))
#filterGrid(param = c(0.336107, 3.265036, 0.147413, 0.264855, -0.065081, 1.188861))
setwd("~/Dropbox/PhD/Meetings and work/August 17/DeBugTheta1_2")

data = scan("drf1.dat")
datagrid = matrix(data, ncol = 2, byrow=TRUE)

filterGrid=function(param=c(0.1, 3, 0.1, 0.1, 0, 1.0),a=c(36,5),C=diag(c(1,0.1)),F=c(1,0),x=xdatagrid)
{
  n=dim(x)[1] # rows = no. obs
  mt=rep(0,2); Vt=matrix(0, ncol=2, nrow=2)
  ll=0
  ll=ll+dnorm(x[1,1],t(F)%*%a,sqrt(t(F)%*%C%*%F+param[6]),log=T)
  #update posterior
  a=a+C%*%F%*%solve(t(F)%*%C%*%F+param[6]^2)%*%(x[1,1]-t(F)%*%a)
  C=C-C%*%F%*%solve(t(F)%*%C%*%F+param[6]^2)%*%t(F)%*%C
  for(i in 1:(n-1))
  {
    tm=x[i,2]
    t=x[(i+1),2]
    mt[1]=m1bt(param[1:5],a,t,tm)
    mt[2]=m2bt(param[1:5],a,t,tm)
    Vt[1,1]=v11bt(param[1:5],c(C[1,1],C[1,2],C[2,2]),t,tm)
    Vt[2,2]=v22bt(param[1:5],c(C[1,1],C[1,2],C[2,2]),t,tm)
    Vt[1,2]=v12bt(param[1:5],c(C[1,1],C[1,2],C[2,2]),t,tm)
    Vt[2,1]=Vt[1,2]
    ll=ll+dnorm(x[i+1,1],t(F)%*%mt,sqrt(t(F)%*%Vt%*%F+param[6]),log=T)
    a=mt+Vt%*%F%*%solve(t(F)%*%Vt%*%F+param[6]^2)%*%(x[i+1,1]-t(F)%*%mt)
    C1 = Vt%*%F%*%solve(t(F)%*%Vt%*%F+param[6]^2)%*%t(F)
    C2 = Vt%*%F%*%solve(t(F)%*%Vt%*%F+param[6]^2)%*%t(F)%*%Vt
    C=Vt-Vt%*%F%*%solve(t(F)%*%Vt%*%F+param[6]^2)%*%t(F)%*%Vt
    #print(Vt)
    # if(is.na(Vt[1,1])){
    #   print(t)
    #   print(tm)
    #   print(i)
    #   break
    # }
    
    #if(i==2953){
     # print(mt)
      #print(Vt)
    #}
    #if(i==2954){
     # print(mt)
      #print(Vt)
      #break
    #}
  }
  # print(ll)
  return(ll)
}

p = c(0.716619, 2.235588, 0.797603, 0.354695, 1.287942, 0.207164)

filterGrid(param = p, x = datagrid)

#logprior(param=c(0.1, 3, 1, 1, 0, 1.0))
#param=c(0.1, 3, 0.1, 0.1, 0, 1.0)
#set.seed(1)
#filterGrid()
v11tb()

logprior=function(param)
{
  dnorm((param[1]),-2,1,log=T)+dnorm(param[2],1.2,0.3,log=T)+dinvgamma(exp(param[3]),3,1, log=TRUE)+dinvgamma(exp(param[4]),3,1, log = TRUE)+dinvgamma(exp(param[6]),3,3, log = TRUE)+dunif(param[5],-pi,pi, log = T)+ sum(exp(param[3:4]))+exp(param[6])
}

#curr=c(log(param[1:4]),param[5],log(param[6]))
#logprior(curr)

mcmcGrid=function(iters=1000,data=datagrid,init=c(0.1, 3, 0.1, 0.1, 0, 1.0),sigtune=diag(rep(0.01,6)),x0=c(36,5),xvar=c(1,0.1))
{
  n=dim(data)[1] # no. data points
  p=length(init) # no. params
  mat=matrix(0,ncol=p+1,nrow=iters)
  mat[1,]=c(init,0)
  #print(mat[1,])
  curr=c(log(init[1:4]),init[5],log(init[6]))
  #print(curr)
  margllcurr=-10000000
  for(i in 2:iters)
  {
    #print("I get here")
    can=rmvn(curr,sigtune)
    if(can[5]<(-pi)){
      can[5] = 2*pi + can[5]
    }
    if(can[5]>pi){
      can[5] = -2*pi + can[5]
    }
    
    margllcan=filterGrid(param=c(exp(can[1:4]),can[5],exp(can[6])),a=x0,C=diag(xvar),F=c(1,0),x=data)
    laprob=margllcan-margllcurr+logprior(can)-logprior(curr)
    rej=0
    if(is.na(margllcan))
    {
      rej=1
    }
    if(rej==0){
      if(log(runif(1))<laprob)
      {
        margllcurr=margllcan
        curr=can
      }
    } 
    mat[i,]=c(exp(curr[1:4]),curr[5],exp(curr[6]),margllcurr)
  }
  return(mat[(2:iters),])
}

#set.seed(1)
outG=mcmcGrid(iters=10000,sigtune=diag(rep(0.01,6)))
#out=mcmcGrid(iters=10000,sigtune=diag(rep(0.01,6)))
outG=mcmcGrid(iters=100000,sigtune = postvar)
mat = matrix(pilot, ncol = 7, byrow=TRUE)
out=cbind(log(outG[,1:4]),outG[,5],log(outG[,6]))
postvar=(1/6)*(2.38^2)*var(out)
write(t(postvar), "postvar.dat", ncolumns = 6)
