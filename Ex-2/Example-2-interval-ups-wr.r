
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
    alpha   <- uniroot(fun, interval=c(low, up), tol=1e-8)$root  
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
        out = simu.pps(nrep, popu, prop, n, theta)  
        cov.all  = rbind(cov.all,  c(n, rho, model, out)) 
        print( c(n, rho, model, out)) 
    } 
  } 
}

colnames(cov.all)=c('n', 'rho', 'model', 'IPW-cov', 'SIPW-cov', 'ELW-cov',  'IPW-width','SIPW-width', 'ELW-width')
cov.all

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)

 
 


> cov.all
        n rho model IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
 [1,] 250 0.2     1   90.76    94.14   91.68     0.540      0.372     0.328
 [2,] 250 0.2     2   91.40    92.54   91.12     0.552      0.438     0.344
 [3,] 250 0.2     3   88.78    94.04   91.94     1.965      0.371     0.328
 [4,] 250 0.2     4   88.98    92.30   91.64     1.973      0.433     0.343
 [5,] 250 0.8     1   94.78    88.50   89.90     0.363      0.436     0.255
 [6,] 250 0.8     2   99.82    90.50   91.76     0.466      1.024     0.485
 [7,] 250 0.8     3   89.68    88.64   90.24     1.830      0.433     0.255
 [8,] 250 0.8     4   90.88    90.06   91.14     1.819      1.021     0.485
 [9,] 500 0.2     1   93.22    94.14   91.50     0.431      0.280     0.244
[10,] 500 0.2     2   94.42    92.96   91.82     0.428      0.330     0.253
[11,] 500 0.2     3   91.36    94.68   92.56     1.518      0.275     0.243
[12,] 500 0.2     4   91.68    93.00   91.84     1.530      0.326     0.253
[13,] 500 0.8     1   98.06    89.64   90.66     0.301      0.331     0.182
[14,] 500 0.8     2   99.98    90.78   91.86     0.410      0.769     0.328
[15,] 500 0.8     3   93.12    90.64   90.84     1.469      0.332     0.182
[16,] 500 0.8     4   93.94    90.56   91.82     1.474      0.769     0.328
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
 用户  系统  流逝 
50.45  1.22 51.94 
> 


> cov.all
        n rho model IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
 [1,] 250 0.2     1   89.98    94.46   91.66     0.569      0.385     0.339
 [2,] 250 0.2     2   90.98    91.60   90.34     0.557      0.449     0.353
 [3,] 250 0.2     3   88.38    94.38   91.64     1.976      0.382     0.338
 [4,] 250 0.2     4   88.30    92.44   91.14     1.969      0.447     0.354
 [5,] 250 0.8     1   94.46    88.38   89.50     0.382      0.436     0.260
 [6,] 250 0.8     2   99.72    89.48   90.58     0.464      1.021     0.486
 [7,] 250 0.8     3   90.26    88.34   88.88     1.930      0.438     0.260
 [8,] 250 0.8     4   90.02    89.32   90.60     1.839      1.021     0.485
 [9,] 500 0.2     1   92.34    94.48   91.34     0.433      0.289     0.253
[10,] 500 0.2     2   93.54    92.96   91.80     0.451      0.333     0.262
[11,] 500 0.2     3   91.60    94.52   91.28     1.560      0.288     0.253
[12,] 500 0.2     4   91.70    93.06   92.28     1.558      0.334     0.263
[13,] 500 0.8     1   97.52    89.46   90.64     0.311      0.333     0.188
[14,] 500 0.8     2  100.00    90.86   91.94     0.419      0.762     0.330
[15,] 500 0.8     3   92.18    89.58   91.32     1.471      0.335     0.187
[16,] 500 0.8     4   94.34    90.56   91.04     1.468      0.775     0.330
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
 用户  系统  流逝 
56.49  0.73 57.90 
> 
>  

