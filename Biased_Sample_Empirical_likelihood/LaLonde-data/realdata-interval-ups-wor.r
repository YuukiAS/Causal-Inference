
###########################      Poisson sampling     ########################


#####################   IPW estimator    #####################################
ipw <- function(dat)
{ 
    N = length(dat$y)
    wt = dat$D/dat$ex
    the  = sum(wt*dat$y)/N 
    Sig  = sum(dat$y^2*wt^2)/N -  sum(dat$y^2*wt)/N 
    c(the, Sig) 
}

###  stablized  IPW estimator    #############################################
###  Kang and Schafer, 2007, Statistics Science 
###        Hajek estimator          ##########################################
sipw <- function(dat)
{   wt = dat$D/dat$ex
    the  = sum(wt*dat$y)/sum(wt)

    Bgg = sum(wt^2*dat$y^2)/sum(wt)
    Bg1 = sum(wt^2*dat$y )/sum(wt) 
    B11 = sum(wt^2)/sum(wt)
    B2  = sum(wt*dat$y^2)/sum(wt)

    Sig  =  Bgg-B2 - (Bg1-the)^2/(B11-1)+(Bg1-the*B11)^2/(B11-1)
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
    alpha   <- uniroot(fun, interval=c(low, up), tol=1e-8)$root  
    lambda  =  (N/n-1)/(1-alpha)      
    tmp     =  1+lambda*(pw-alpha) 
    prob.el =  1/(n*tmp) 
   
    theta.elw = sum(y*prob.el)   

    B11 = N*sum(prob.el^2 )
    Bg1=  N*sum(prob.el^2*y)
    Bgg=  N*sum(prob.el^2*y^2)

    B2 =  sum(y^2*prob.el) 

    Sig = Bgg - B2  - (Bg1-theta.elw)^2/(B11-1)

    c(theta.elw,   Sig) 
 }


  
  



#####   Povital sampling
simu.pivotal<-function(nrep, popu, prop, theta)
{ 
  result=NULL
  N = length(popu)
  n= round(sum(prop)) 
  for( i in 1:nrep){  
      pik= sampling::inclusionprobabilities(prop, n)  # each subject has different prop, result is deterministic
      index =  sampling::UPpivotal(pik)        
      dat=list(y=popu, ex=prop, D=(index>0.5))

      out.ipw=ipw(dat) 
      stat.ipw  = sqrt(N)*(out.ipw[1] - theta)/sqrt(out.ipw[2]) 

      out.sipw=sipw(dat) 
      stat.sipw  = sqrt(N)*(out.sipw[1] - theta)/sqrt(out.sipw[2]) 

      out.elw   = elw(dat) 
      stat.elw  = sqrt(N)*(out.elw[1]-theta) /sqrt(out.elw[2])     
 
      stat = c(stat.ipw, stat.sipw,stat.elw)
      index= ( abs(stat)<1.96 )   

      width= 2*1.96*c( sqrt(out.ipw[2]), sqrt(out.sipw[2]), sqrt(out.elw[2]) )/sqrt(N)


      result = rbind(result, c(100*index, width))
   }    
   round(apply(result, 2, mean),  3) 
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

a = 2 ###  or  0

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
  out = simu.pivotal(nrep, popu, prop, theta)  
        cov.all  = rbind(cov.all,  c(n,  out)) 
        print( c(n,  out)) 
 
}

colnames(cov.all)=c('n',  'IPW-cov','SIPW-cov', 'ELW-cov', 'IPW-width', 'SIPW-width', 'ELW-width')
cov.all

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)


 
####################################################################################################   a= 0

> cov.all
       n IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
[1,] 100   99.98    93.20   95.82     0.985      0.755     0.591
[2,] 200   99.98    93.70   96.70     0.683      0.560     0.424
[3,] 400   99.98    94.44   97.30     0.469      0.418     0.298
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
  �û�   ϵͳ   ���� 
181.52   2.09 187.51 
> 


###   a=2 

> cov.all
       n IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
[1,] 100   99.98    93.22   95.96     1.955      0.746     0.590
[2,] 200   99.98    93.46   96.72     1.408      0.568     0.426
[3,] 400   99.96    94.28   96.92     0.981      0.416     0.298
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
  �û�   ϵͳ   ���� 
179.22   2.14 187.09 

 
