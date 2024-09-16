
#####################   IPW estimator    #####################################
ipw <- function(dat)
{ 
  sum(dat$D*dat$y/dat$ex)/length(dat$y) 
}


###  stablized  IPW estimator    #############################################
###  Kang and Schafer, 2007, Statistics Science 
###        Hajek estimator          ##########################################
sipw <- function(dat)
{  
  sum(dat$D*dat$y/dat$ex)/sum(dat$D/dat$ex) 
}


###  Improved IPW estimator of Zong, Zhu and Zou (2019)   #################### 
###  need to know all the population sampling weights
###              to calculate crit.  
ipw.zzz <- function(dat)
{ 
    N   = length(dat$y)
    index = 1:N
    tmp = sort( dat$ex )
    i = max( which( tmp < 1/(index+1) ) ) 
    crit = tmp[i] 
    
    ex.star = dat$ex  
    ex.star[which(dat$ex <=crit)] = crit
    sum(dat$D*dat$y/ex.star)/N 
}

#####  Trimmed IPW estimator of Crump et al (2009)   #########################    
ipw.chim <- function(dat)
{ 
    tempfun <- function(a) { 
        u = dat$ex*(1- dat$ex)
        v = a*(1-a)  
        z = (u>=v)  
        a + (2*v*sum( z/u )  >  sum(z) )
    }   
    a =  optimize(tempfun, lower=1e-5, upper=0.5, maximum=FALSE)$minimum  

    N   = length(dat$y)
    y   = rep(dat$y,  dat$D)
    pw  = rep(dat$ex, dat$D) 

    ind = ( pw>a & pw< 1-a)
    y.star  = y[ind]
    pw.star = pw[ind]
    sum(y.star/pw.star)/sum(1/pw.star)   
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
    alpha   <- uniroot(fun, interval=c(low, up), tol=eeps)$root
    lambda  =  (N/n-1)/(1-alpha)      
    tmp     =  1+lambda*(pw-alpha) 
    prob.el =  1/(n*tmp) 
   
    sum(y*prob.el)  
 }


 

####  Robust IPW estimator of  Wang and Ma, 2017, JASA   ##################### 
 
ipw.trim<-function(dat)
{ N = length(dat$D)
  n = sum(dat$D) 

  fun1<-function(t, s){  2*n*t^s*mean(dat$ex<=t) - 1 }
  bn1=uniroot(fun1, c(0, 1), s=1)$root 
  bn2=uniroot(fun1, c(0, 1), s=2)$root 

  fun2<-function(t){  n*t^5*mean(dat$ex<=t) - 1  }
  hn=uniroot(fun2, c(0, 1))$root  

  yy = rep(dat$y, dat$D)
  pw = rep(dat$ex,dat$D)
  xx = cbind(rep(1,n),  pw)
  ww = (pw<hn)
  xx1 = cbind(ww,  ww*pw)
  bet = ginv(t(xx1)%*%xx)%*%t(xx1)%*%yy 
  xx.all = cbind(rep(1, N), dat$ex)

  Bnb1 =  - mean( (xx.all%*%bet)*(dat$ex<bn1) )
  ipw.trim1 =  mean(yy*(pw >= bn1)/pw) *n/N            
  ipw.bc1 = ipw.trim1-Bnb1 

  Bnb2 =  - mean( (xx.all%*%bet)*(dat$ex<bn2) )
  ipw.trim2 =  mean(yy*(pw >= bn2)/pw)*n/N
  ipw.bc2 = ipw.trim2 -Bnb2

  c(ipw.bc1, ipw.bc2)
}

 
 
##############################################################################    
simu.one<-function(nrep, gam, N, fun_index, coef)
{
  datt0 = data.gen(gam, 200000, fun_index, coef)
  my   = mean(datt0$y) 
  result=NULL
  for( i in 1:nrep){ 
     # set.seed(i*10+5)
      dat0    = data.gen(gam, N, fun_index, coef) 

      obj=glm(dat0$D~dat0$x, family=binomial(link="logit"))
      ex = obj$fitted.values
      dat = list(y=dat0$y, D=dat0$D, ex=ex, x=dat0$x) 
      est    = c( ipw(dat), sipw(dat), ipw.zzz(dat), 
                  ipw.chim(dat), ipw.trim(dat),  elw(dat))   ### 
      result = rbind(result, est)
   }   
  rmse  = sqrt( apply( (result-my)^2, 2, mean))  
  rmse= sqrt(N)*rmse 
} 

 

###############################################################################
####    Generate data 
data.gen <- function(gam, N, fun_index, coef)
{
   x  = rexp(N, gam)     
   ex = plogis(-x) 
   eta = rnorm(N)    # epsilon
   reg = (fun_index==1)*log(1+x) + (fun_index==2)*((x-1)^2/4+(x-3)^2/4) + 
           (fun_index==3)*(5+ log(1+x)) + (fun_index==4)*((x-1)^2/4+(x-3)^2/4+5)
   y =  reg + eta*coef    
   D   = rbinom(N, 1, ex)
   dat = list(y=y, x=x, D=D)
   return(dat)
}
 
###############################################################################
####    Generate data 
data.gen1 <- function(gam, N, fun_index, coef)
{
   x  = (rnorm(N)+1)*gam     
   ex = plogis(-2+x) 
   eta = rnorm(N)  
   reg = (fun_index==1)*cos(2*pi*x^2) + (fun_index==2)*(x^2) + 
           (fun_index==3)*(5+cos(2*pi*x^2))+ (fun_index==4)*(x^2+5)
   y =  reg + eta*coef    
   D   = rbinom(N, 1, ex)
   dat = list(y=y, x=x, D=D)
   return(dat)
}

 


################################################################################
start.time=proc.time() 
library(sampling)
library(MASS)
eeps=sqrt(.Machine$double.eps)    #### root of machine precision

nrep=5000 
NN=c(500, 2000)
gamm=c(1, 1/2)
cff = c(1, 1/4)  

 
rmse.all = NULL
for(i in 1:2){   ###  3
   N = NN[i]
   for(k in 1:2){   ### 4
      gam = gamm[k] 
      for(a in 1:2){
        coef = cff[a]
        for(j in 1:4){  ## 4
            fun_index = j 
            rmse0 = simu.one(nrep, gam, N, fun_index, coef)
            rmse = c(N, gam, coef, fun_index, rmse0)
             print(round(rmse, 2))
            rmse.all = rbind(rmse.all, rmse)
         }
      }
   }
}

  
result= round(rmse.all,   2)
colnames(result)=c("N", "gam", "coef", "func", "IPW", "SIPW", "ZZZ", "CHIM", "MW1", "MW2", "ELW")
result

end.time=proc.time()
cpu.time = end.time-start.time
print(cpu.time)





################################################################################
################################################################################
###    simulation results  
> result
        N gam coef func     IPW  SIPW   ZZZ  CHIM   MW1   MW2   ELW
rmse  500 1.0 1.00    1    4.00  2.93  2.57  2.93  2.39  2.36  2.45
rmse  500 1.0 1.00    2   23.55  7.77  3.50  7.77  4.04  4.34  3.50
rmse  500 1.0 1.00    3   11.69  2.89  5.70  2.89  4.15  3.91  2.39
rmse  500 1.0 1.00    4   19.25  6.23  6.84  6.23  4.96  5.10  3.55
rmse  500 1.0 0.25    1    3.42  1.63  1.41  1.63  1.03  0.93  0.95
rmse  500 1.0 0.25    2   42.84  7.72  2.91  7.72  3.49  3.83  2.86
rmse  500 1.0 0.25    3   13.75  1.68  5.21  1.68  3.59  3.30  0.97
rmse  500 1.0 0.25    4   16.16  4.35  6.35  4.35  4.36  4.39  2.53
rmse  500 0.5 1.00    1   24.60  5.62  5.86  5.62  4.01  3.81  4.55
rmse  500 0.5 1.00    2  164.77 32.92 28.84 32.92 32.29 33.50 26.46
rmse  500 0.5 1.00    3   60.77  5.74 17.71  5.74  7.56  5.34  4.54
rmse  500 0.5 1.00    4  378.83 34.80 40.33 34.80 33.07 33.59 26.60
rmse  500 0.5 0.25    1   13.15  3.94  4.92  3.94  2.39  2.07  2.19
rmse  500 0.5 0.25    2  106.84 30.86 28.78 30.86 32.32 33.57 26.21
rmse  500 0.5 0.25    3  111.43  3.81 17.54  3.78  6.80  4.39  2.19
rmse  500 0.5 0.25    4  134.39 30.71 40.49 30.71 32.23 33.16 26.06
rmse 2000 1.0 1.00    1    4.94  3.42  2.96  3.42  2.63  2.51  2.67
rmse 2000 1.0 1.00    2   14.68  9.59  4.93  9.59  5.76  6.84  4.72
rmse 2000 1.0 1.00    3   17.44  3.68  6.64  3.68  4.70  4.00  2.59
rmse 2000 1.0 1.00    4   63.04 16.74  8.96 16.74  6.88  7.16  4.89
rmse 2000 1.0 0.25    1    5.25  2.48  1.77  2.48  1.27  1.08  1.07
rmse 2000 1.0 0.25    2   18.69 10.48  4.34 10.48  5.21  6.37  4.11
rmse 2000 1.0 0.25    3   16.06  2.29  6.40  2.29  4.23  3.39  1.06
rmse 2000 1.0 0.25    4   27.11 10.50  8.76 10.50  6.51  6.82  4.07
rmse 2000 0.5 1.00    1   48.22  9.10  8.09  8.91  5.19  4.77  5.61
rmse 2000 0.5 1.00    2  426.28 70.79 49.41 67.50 56.49 62.57 43.49
rmse 2000 0.5 1.00    3   75.96  9.07 24.90  9.07 10.97  6.81  5.64
rmse 2000 0.5 1.00    4  234.17 65.41 65.43 65.41 57.97 62.84 43.61
rmse 2000 0.5 0.25    1   36.09  6.57  7.01  6.50  3.70  3.37  2.54
rmse 2000 0.5 0.25    2  687.50 74.54 49.69 65.31 56.79 62.92 43.67
rmse 2000 0.5 0.25    3   77.50  6.33 24.37  6.33 10.42  5.78  2.58
rmse 2000 0.5 0.25    4 1186.83 75.96 65.42 68.23 57.91 62.84 43.68

> 
> end.time=proc.time()
> cpu.time = end.time-start.time
> print(cpu.time)

1397.28 
> 
> 

