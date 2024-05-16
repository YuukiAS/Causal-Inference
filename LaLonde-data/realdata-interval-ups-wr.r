
###########################      PPS sampling     ########################

#####################   IPW estimator    #####################################
ipw <- function(dat)
{   alp = mean(dat$D)
    N = length(dat$y)
    wt = dat$D/dat$ex
    the  = sum(wt*dat$y)/N 
    Sig  = alp*sum(dat$y^2*wt^2)/N -  the^2 
    c(the, Sig) 
}

###  stablized  IPW estimator    #############################################
###  Kang and Schafer, 2007, Statistics Science 
###        Hajek estimator          ##########################################
sipw <- function(dat)
{   
   
    the = sum(dat$y*dat$D/dat$ex)/sum(dat$D/dat$ex) 
    alp = mean(dat$D)
    Bgg = sum(dat$y^2*dat$D/dat$ex^2)/sum(dat$D/dat$ex) 
    Bg1 = sum(dat$y*dat$D/dat$ex^2)/sum(dat$D/dat$ex) 
    B11 = sum(dat$D/dat$ex^2)/sum(dat$D/dat$ex) 
  

    Sig  = alp*(Bgg - 2*the*Bg1 + the^2*B11) 
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

    B11 =  N*sum(prob.el^2 )
    Bg1=   N*sum(prob.el^2*y)
    Bgg=   N*sum(prob.el^2*y^2)
    alp = n/N

    Sig1 = alp*Bgg - theta.elw^2 
    Sig2 = (1-alp)^2*(theta.elw*B11- Bg1)^2/((B11*alp-1)*(B11-1)^2)
    Sig3 = (alp*Bg1-theta.elw)^2/(alp*B11-1)
    Sig = Sig1 +Sig2-Sig3

    Sig=max(0, Sig)
    c(theta.elw, Sig) 
 }


#####   Poisson sampling
simu.pps<-function(nrep, popu, prop, n, theta)
{ 
   result=NULL
   N = length(popu)
   prop = prop/sum(prop)
  
   for( i in 1:nrep){  
      index = rmultinom(1, n, prop)
      dat=list(y=popu, ex=n*prop, D=index)       

      out.ipw=ipw(dat) 
      stat.ipw  = sqrt(n)*(out.ipw[1] - theta)/sqrt(out.ipw[2]) 

      out.sipw=sipw(dat) 
      stat.sipw  = sqrt(n)*(out.sipw[1] - theta)/sqrt(out.sipw[2]) 

      out.elw   = elw(dat) 
      stat.elw  = sqrt(n)*(out.elw[1]-theta) /sqrt(out.elw[2])     
 
      stat = c(stat.ipw,stat.sipw, stat.elw) 
      index= ( abs(stat)<1.96 ) 
      width = 2*1.96*sqrt(c(out.ipw[2], out.sipw[2], out.elw[2]))/sqrt(n)
        
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

a = 0 ###  or  2

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
  out = simu.pps(nrep, popu, prop0, n, theta)  
        cov.all  = rbind(cov.all,  c(n,  out)) 
        print( c(n,  out)) 
 
}

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)



##################################################
###   a=2

[1] 100.000  89.940  89.640  90.440   1.112   0.749   0.520
[1] 200.000  93.240  91.320  91.960   0.893   0.572   0.377
[1] 400.000  95.860  91.340  92.060   0.709   0.435   0.270
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
用户 系统 流逝 
7.82 0.69 8.60 


####   a= 0 

[1] 100.000  95.340  89.900  90.540   0.528   0.752   0.521
[1] 200.000  97.340  90.600  91.660   0.416   0.573   0.376
[1] 400.000  98.580  91.580  92.280   0.341   0.430   0.271
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
用户 系统 流逝 
7.29 0.66 7.95 


 