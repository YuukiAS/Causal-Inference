
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
    alpha   <- uniroot(fun, interval=c(low, up), tol=1e-8)$root  
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

N = 3000 
x = runif(N)*2
e = rnorm(N) 
 
prop0 = x/sum(x)
nrep = 5000  
nn= c(250, 500)

rho.all=c(0.2, 0.8)   
cov.all = NULL  


for(i in 1:2){ 
  n = nn[i]
  prop = n*prop0
  for(k in 1:2){
    rho = rho.all[k]
    for(j in 1:4){
        model = j
        reg  = (model==1)*sqrt(3)*rho* x  +  (model==2)*sqrt(3)*rho*(x+x^2) + 
                 (model==3)*(sqrt(3)*rho* x  + 5 ) + (model==4)*(sqrt(3)*rho*(x+x^2) + 5)

        popu  = reg +  sqrt(3)*sqrt(1- rho^2)*abs(e)  
        theta = mean(popu)   
        out = simu.poisson(nrep, popu, prop, theta)  
        cov.all  = rbind(cov.all,  c(n, rho, model, out)) 
        print( c(n, rho, model, out)) 
    } 
  } 
}

colnames(cov.all)=c('n', 'rho', 'model', 'IPW-cov','SIPW-cov', 'ELW-cov', 'IPW-width','SIPW-width', 'ELW-width')
cov.all

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)

  


> cov.all
        n rho model IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
 [1,] 250 0.2     1   93.24    93.82   91.94     0.661      0.376     0.341
 [2,] 250 0.2     2   94.12    92.18   92.04     0.727      0.446     0.373
 [3,] 250 0.2     3   92.40    94.26   92.60     2.512      0.374     0.335
 [4,] 250 0.2     4   92.98    92.18   92.48     2.572      0.453     0.374
 [5,] 250 0.8     1   95.50    88.50   91.04     0.622      0.449     0.320
 [6,] 250 0.8     2   95.58    89.84   93.66     1.006      1.039     0.680
 [7,] 250 0.8     3   93.06    88.68   90.98     2.459      0.449     0.320
 [8,] 250 0.8     4   94.32    91.12   93.50     2.724      1.025     0.679
 [9,] 500 0.2     1   93.28    95.34   93.86     0.476      0.280     0.251
[10,] 500 0.2     2   94.26    93.18   93.52     0.507      0.333     0.275
[11,] 500 0.2     3   93.06    95.16   93.96     1.771      0.278     0.252
[12,] 500 0.2     4   93.94    93.78   93.58     1.811      0.334     0.274
[13,] 500 0.8     1   94.72    90.90   93.52     0.431      0.343     0.232
[14,] 500 0.8     2   96.04    91.84   95.34     0.684      0.776     0.479
[15,] 500 0.8     3   93.80    90.26   93.54     1.752      0.340     0.233
[16,] 500 0.8     4   94.66    91.78   94.92     1.928      0.778     0.479
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
 用户  系统  流逝 
37.31  0.98 38.47 
> 
> 

