
#####################   IPW estimator    #####################################
ipw <- function(dat)
{ 
  # dat$ex is the propensity score P(D|L)
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
    alpha   <- uniroot(fun, interval=c(low, up), tol=1e-8)$root  
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
  pw = rep(dat$ex, dat$D)
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

 




#####   Poisson sampling
simu.poisson<-function(nrep, popu, prop, theta)
{ 
  result=NULL
  n0 = round(sum(prop))
  for( i in 1:nrep){ 
     # set.seed(i*10+5)
      d=rbinom(N, rep(1, N), prop)
      if(sum(d)<=5) d[1:5]=1  
      
      dat=list(y=popu, ex=prop, D=d)
      est    = c( ipw(dat), sipw(dat),  ipw.zzz(dat),  ipw.trim(dat), elw(dat)) 
      result = rbind(result, est)
   }      
 
   ##############################################################
   ############          boxplot     ############################
  
   txt = paste('Poisson-n=', n, sep='')
   flnm = paste('d:/', txt, '.pdf',   sep='')
   pdf(flnm)
   xnm = c('SIPW', 'ZZZ', 'MW1', 'MW2', 'ELW' )
   boxplot(result[, -1], xaxt='n', xlab='',
               main=paste('Poisson sampling: n=', n,    sep='')
           )
   abline(h=theta)
   axis(1, 1:5, labels=xnm)
   dev.off()
   ##############################################################

   rmse  = round( sqrt( N*apply(  (result-theta)^2, 2, mean)) , 2) 
   colnm = c('IPW','SIPW', 'ZZZ',   'MW1',  'MW2', 'ELW')
   rmse = matrix(rmse, 1)
   colnames(rmse) <-  colnm
   list(result=result, rmse= rmse )
}


#####   Pivotal sampling
simu.pivotal<-function(nrep, popu, prop, theta) 
{ result=NULL
  n= round(sum(prop)) 
  for( i in 1:nrep){  
      set.seed(10+i*757)
      pik=inclusionprobabilities(prop, n)
      index =  UPpivotal(pik)        
      dat=list(y=popu, ex=prop, D=(index>0.5))
      est    = c( ipw(dat), sipw(dat),  ipw.zzz(dat),   ipw.trim(dat), elw(dat)) 
      result = rbind(result, est)
    }  
    
   ##############################################################
   ############          boxplot     ############################
  
   txt = paste('Pivotal-n=', n,  sep='')
   flnm = paste('d:/', txt, '.pdf',   sep='')
   pdf(flnm)
   xnm = c('SIPW', 'ZZZ',   'MW1', 'MW2', 'ELW' )
   boxplot(result[, -1], xaxt='n', xlab='',
               main=paste('Pivotal sampling: n=', n,   sep='')
           )
   abline(h=theta)
   axis(1, 1:5, labels=xnm)
   dev.off()
   ##############################################################

   rmse  = round( sqrt( N*apply((result-theta)^2, 2, mean)), 2) 
   colnm = c('IPW','SIPW', 'ZZZ',   'MW1',  'MW2', 'ELW')
   rmse = matrix(rmse, 1)
   colnames(rmse) <-  colnm 
   print(rmse)
   list(result=result, rmse =rmse )
}


#####   PPS sampling
simu.pps<-function(nrep, popu, prop, n, theta)
{  result=NULL 
   N = length(popu) 
   prop = prop/sum(prop) 
   for( i in 1:nrep){  
      index = rmultinom(1, n, prop) 
      dat=list(y=popu, ex=n*prop, D=index)
      est    = c( ipw(dat), sipw(dat),  ipw.zzz(dat), ipw.trim(dat), elw(dat))  
      result = rbind(result, est)
    }  
 

   ##############################################################
   ############          boxplot     ############################
 
   txt = paste('PPS-n=', n,   sep='')
   flnm = paste('d:/', txt, '.pdf',   sep='')
   pdf(flnm)
   xnm = c('SIPW', 'ZZZ',   'MW1', 'MW2', 'ELW' )
   boxplot(result[, -1], xaxt='n', xlab='',
               main=paste('PPS sampling:  n=', n,     sep='')
           )
   abline(h=theta)
   axis(1, 1:5, labels=xnm)
   dev.off()
   ##############################################################

   rmse  = round( sqrt(N*apply(  (result-theta)^2, 2, mean)) , 2) 
   colnm = c('IPW','SIPW', 'ZZZ',   'MW1',  'MW2', 'ELW')
   rmse = matrix(rmse, 1)
   colnames(rmse) <-  colnm
   list(result=result, rmse=  rmse )
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
ind = (x0> 0)
y = y0[ind]
x = x0[ind]  
N=sum(ind)

a = 2 ###  0 or  2

popu =  y/10000 + a
prop0 = x/sum(x) 

 
nrep = 1000     
nn    = c(100, 200, 400)

theta= mean(popu)   
rmse1=rmse2=rmse3=NULL  

for(i in  1:3){ 
  n = nn[i]
  prop = n*prop0
  prop[prop>0.99] = 0.99
 
  out1 = simu.poisson(nrep, popu, prop, theta)$rmse  
   rmse1 = rbind(rmse1,  c(n,   out1)) 
  print(out1)

  out2 = simu.pivotal(nrep, popu, prop, theta)$rmse  
   rmse2 = rbind(rmse2,  c(n,   out2)) 
  print(out2)


  out3 = simu.pps(nrep, popu, prop, n, theta)$rmse  
  rmse3 = rbind(rmse3,  c(n,   out3)) 
  print(out3)

}


rmse.poisson = round(rmse1, 2)
rmse.pivotal = round(rmse2, 2)
rmse.pps     = round(rmse3, 2)


colnm = c('n','IPW','SIPW', 'ZZZ',   'MW1',  'MW2', 'ELW')
colnames(rmse.poisson) = colnm 
colnames(rmse.pivotal) = colnm 
colnames(rmse.pps) = colnm 

rmse.poisson
rmse.pivotal
rmse.pps

end.time=proc.time()   
cpu.time = end.time-start.time
print(cpu.time)



################################################################################
################################################################################
###    simulation results  

####  a = 0

> rmse.poisson
       n   IPW  SIPW   ZZZ   MW1  MW2  ELW
[1,] 100 13.89 11.16 11.71 11.55 8.38 8.39
[2,] 200  9.35  8.41  8.44  7.85 7.47 6.14
[3,] 400  6.11  6.06  5.66  5.49 5.27 4.07
> rmse.pivotal
       n  IPW SIPW  ZZZ  MW1  MW2  ELW
[1,] 100 8.39 9.99 5.35 6.13 7.64 6.94
[2,] 200 5.07 7.15 3.91 3.62 5.21 4.66
[3,] 400 3.41 5.50 2.82 2.60 2.65 3.31
> rmse.pps
       n  IPW  SIPW  ZZZ  MW1  MW2  ELW
[1,] 100 8.00 11.09 5.69 6.56 8.17 7.62
[2,] 200 5.46  8.70 4.13 4.01 5.63 5.51
[3,] 400 3.68  6.42 3.11 2.86 3.05 3.77
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
 �û�  ϵͳ  ���� 
46.67  1.48 48.23 


## a=2
> rmse.poisson
       n   IPW  SIPW   ZZZ   MW1   MW2  ELW
[1,] 100 27.60 10.28 23.17 18.99  8.99 8.18
[2,] 200 19.27  8.41 16.33 13.89 11.03 6.14
[3,] 400 13.40  6.06 11.54 10.35  8.85 4.07
> rmse.pivotal
       n   IPW SIPW   ZZZ  MW1  MW2  ELW
[1,] 100 20.43 9.99 10.17 7.93 8.47 6.93
[2,] 200 12.17 7.15  7.63 5.30 6.55 4.66
[3,] 400  8.50 5.50  5.83 4.32 3.48 3.31
> rmse.pps
       n   IPW  SIPW   ZZZ  MW1  MW2  ELW
[1,] 100 19.95 11.09 11.61 8.45 9.13 7.62
[2,] 200 13.86  8.70  8.19 5.94 7.06 5.51
[3,] 400  9.26  6.42  6.56 4.92 3.98 3.77
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
 �û�  ϵͳ  ���� 
52.39  1.10 53.57 
> 

###  a=10

> rmse.poisson
       n    IPW  SIPW   ZZZ   MW1   MW2  ELW
[1,] 100 106.46 11.18 71.34 49.63 14.98 8.44
[2,] 200  65.79  8.41 51.58 40.25 27.02 6.14
[3,] 400  46.55  6.06 37.66 31.83 24.46 4.07
> rmse.pivotal
       n   IPW SIPW   ZZZ   MW1   MW2  ELW
[1,] 100 76.91 9.99 37.15 20.71 13.11 6.94
[2,] 200 46.69 7.15 27.69 16.82 14.20 4.66
[3,] 400 33.65 5.50 21.70 14.56  9.38 3.31
> rmse.pps
       n   IPW  SIPW   ZZZ   MW1   MW2  ELW
[1,] 100 77.01 11.09 43.66 22.14 14.33 7.62
[2,] 200 56.61  8.70 31.36 19.33 15.39 5.51
[3,] 400 37.61  6.42 25.12 17.80 10.85 3.77
> 
> end.time=proc.time()
> cpu.time = end.time-start.time
> print(cpu.time)
 �û�  ϵͳ  ���� 
46.06  1.56 47.77 
