
###########################      Poisson sampling     ########################

#####################   IPW estimator    #####################################
ipw <- function(dat)
{ 
    N = length(dat$y)
    wt = dat$D/dat$ex
    the  = sum(wt*dat$y)/N 
    Sig  = sum(dat$y^2*wt^2)/N -  the^2
    c(the, Sig) 
}


###  stablized  IPW estimator    #############################################
###  Kang and Schafer, 2007, Statistics Science 
###        Hajek estimator          ##########################################
sipw <- function(dat)
{   wt = dat$D/dat$ex
    the  = sum(wt*dat$y)/sum(wt)
    gc = dat$y - the 
    Sig  = sum(gc^2*wt^2)/sum(wt)
    c(the, Sig)
}


#############################   Our ELW estimator ############################  
elw  <- function(dat)
{ 
    y   = rep(dat$y, dat$D)
    pw  = rep(dat$ex, dat$D) 
    N   = length(dat$y) 

    n   = sum(dat$D) 
    xi  = n/N + (1-n/N)*pw
    low = min(pw);   
    up  = min(xi)-eeps    

    fun     <- function(alpha){  sum((pw-alpha)/(xi-alpha)) }
    alpha   <- uniroot(fun, interval=c(low, up))$root  
    lambda  =  (N/n-1)/(1-alpha)      
    tmp     =  1+lambda*(pw-alpha) 
    prob.el =  1/(n*tmp) 
   
    theta.elw = sum(y*prob.el)   

    B11 = N*sum(prob.el^2 )
    Bg1=  N*sum(prob.el^2*y)
    Bgg=  N*sum(prob.el^2*y^2)

    Sig = Bgg - theta.elw^2 - (Bg1-theta.elw)^2/(B11-1)

    c(theta.elw,   Sig) 
 }


  
  



#####   Poisson sampling
simu.poisson<-function(nrep, popu, prop, theta)
{ N = length(popu)
  result=NULL
  n0 = round(sum(prop))
  for( i in 1:nrep){  
      d=rbinom(N, rep(1, N), prop)
      if(sum(d)<5)  d[1:n0]=1         
      dat=list(y=popu, ex=prop, D=d) 

      out.ipw=ipw(dat) 
      stat.ipw  = sqrt(N)*(out.ipw[1] - theta)/sqrt(out.ipw[2]) 

      out.sipw=sipw(dat) 
      stat.sipw  = sqrt(N)*(out.sipw[1] - theta)/sqrt(out.sipw[2]) 

      out.elw   = elw(dat) 
      stat.elw  = sqrt(N)*(out.elw[1]-theta) /sqrt(out.elw[2])     
 
      stat = c(stat.ipw, stat.sipw,stat.elw)
      width= 2*1.96*c( sqrt(out.ipw[2]),sqrt(out.sipw[2]), sqrt(out.elw[2]) )/sqrt(N)

      index= ( abs(stat)<1.96 )   
      result = rbind(result, c(index*100, width))
   }    
   round(apply(result, 2, mean), 3) 
}

 
 

#################################################################

start.time=proc.time()   ### record the beginning time
library(sampling)
library(MASS)
eeps=sqrt(.Machine$double.eps)   #### root of machine precision



library(cem)
data(LLvsPSID)
dat.raw  = LLvsPSID
y0 = dat.raw$re78
x0 = dat.raw$re75
ind = (x0>0)
y = y0[ind]
x = x0[ind]
N=sum(ind)

a = 0 ###  or  0

popu =  y/10000 + a
prop0 = x/sum(x) 

 
nrep = 5000     
nn    = c(100, 200, 400)
 

 
  
cov.all = NULL  
theta = mean(popu)

for(i in 1:3){ 
  n = nn[i]
  prop = n*prop0
  prop[prop>0.99] = 0.99
  out = simu.poisson(nrep, popu, prop, theta)  
        cov.all  = rbind(cov.all,  c(n,  out)) 
        print( c(n,  out)) 
 
}

colnames(cov.all)=c('n',  'IPW-cov','SIPW-cov', 'ELW-cov', 'IPW-width', 'SIPW-width', 'ELW-width')
cov.all

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)

 

####################################################################
##   a=2

> cov.all
       n IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
[1,] 100   94.02    89.16   92.12     1.980      0.742     0.598
[2,] 200   94.32    91.08   92.66     1.419      0.573     0.440
[3,] 400   94.28    92.78   94.94     0.989      0.432     0.319
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
用户 系统 流逝 
6.24 0.75 7.12 



####   a=0

> cov.all
       n IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
[1,] 100   94.28    89.84   91.26     0.995      0.748     0.597
[2,] 200   95.08    91.28   92.88     0.695      0.568     0.439
[3,] 400   95.16    92.42   94.46     0.482      0.433     0.321
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
用户 系统 流逝 
6.08 0.75 6.84 

