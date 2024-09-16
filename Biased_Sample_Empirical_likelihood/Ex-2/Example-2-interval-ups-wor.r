
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


  
  



#####   Pivotal sampling
simu.pivotal<-function(nrep, popu, prop, theta)
{ 
  result=NULL
  N = length(popu)
  n= round(sum(prop)) 
  for( i in 1:nrep){  
      pik=sampling::inclusionprobabilities(prop, n)
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
        out = simu.pivotal(nrep, popu, prop, theta)  
        cov.all  = rbind(cov.all,  c(n, rho, model, out)) 
        print( c(n, rho, model, out)) 
    } 
  } 
}

colnames(cov.all)=c('n', 'rho', 'model', 'IPW-cov','SIPW-cov', 'ELW-cov', 'IPW-width', 'SIPW-width', 'ELW-width')
cov.all

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)


 
> cov.all
        n rho model IPW-cov SIPW-cov ELW-cov IPW-width SIPW-width ELW-width
 [1,] 250 0.2     1   97.38    94.54   93.14     0.671      0.365     0.339
 [2,] 250 0.2     2   98.98    92.54   93.26     0.746      0.422     0.365
 [3,] 250 0.2     3   98.10    94.56   92.96     2.767      0.370     0.342
 [4,] 250 0.2     4   98.90    92.82   93.18     2.525      0.421     0.364
 [5,] 250 0.8     1  100.00    89.24   94.06     0.634      0.414     0.301
 [6,] 250 0.8     2  100.00    89.96   96.54     1.004      0.975     0.642
 [7,] 250 0.8     3   99.26    89.38   94.40     2.484      0.416     0.302
 [8,] 250 0.8     4  100.00    89.44   96.84     2.726      0.972     0.641
 [9,] 500 0.2     1   96.82    94.34   92.84     0.494      0.270     0.247
[10,] 500 0.2     2   98.66    92.52   93.68     0.528      0.313     0.262
[11,] 500 0.2     3   97.86    94.84   93.40     1.789      0.267     0.246
[12,] 500 0.2     4   98.52    92.48   93.56     1.828      0.311     0.263
[13,] 500 0.8     1   99.98    89.40   94.90     0.438      0.317     0.214
[14,] 500 0.8     2  100.00    91.46   97.48     0.698      0.729     0.436
[15,] 500 0.8     3   98.98    88.80   94.46     1.734      0.310     0.214
[16,] 500 0.8     4   99.88    89.78   97.10     1.936      0.724     0.435
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
  �û�   ϵͳ   ���� 
962.24   2.18 968.30 
