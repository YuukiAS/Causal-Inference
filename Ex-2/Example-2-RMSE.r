
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


#####   Poisson sampling
simu.poisson<-function(nrep, popu, prop, other)
{ 
  result=NULL
  n0 = round(sum(prop))
  for( i in 1:nrep){ 
     # set.seed(i*10+5)
      d=rbinom(N, rep(1, N), prop)
      if(sum(d)==0) d[1:n0]=1        ######防止poisson抽样中没有单元被抽中的情况
      
      dat=list(y=popu, ex=prop, D=d)
      est    = c( ipw(dat), sipw(dat),  ipw.zzz(dat),   elw(dat)) 
      result = rbind(result, est)
   }  
   theta = other[1]
   rho   = other[2]    
   model  = other[3]    


 
   ##############################################################
   ############          boxplot     ############################
   rho1 = 10*rho   
   txt = paste('Poisson-n=', n, '-rho=', rho1, '-Model=', model, sep='')
   flnm = paste('d:/', txt, '.pdf',   sep='')
   pdf(flnm) 
   boxplot(result[, -1], xaxt='n', xlab='',
               main=paste('Poisson sampling: n=', n, ',  rho=', rho, ', Model=', model,  sep='')
           )
   abline(h=theta)
   xnm = c('IPW','SIPW', 'ZZZ',   'ELW')
   axis(1, 1:3, labels=xnm[-1])
   dev.off()
   ##############################################################

   rmse  = round( sqrt( N*apply(  (result-theta)^2, 2, mean)) , 2) 
   colnm = c('IPW','SIPW', 'ZZZ',   'ELW')
   rmse = matrix(rmse, 1)
   colnames(rmse) <-  colnm
   list(result=result, rmse= rmse )
}


#####   Pivotal sampling
simu.pivotal<-function(nrep, popu, prop, other) 
{ result=NULL
  n= round(sum(prop)) 
  for( i in 1:nrep){  
      set.seed(10+i*757)
      pik=inclusionprobabilities(prop, n)
      index =  UPpivotal(pik)        
      dat=list(y=popu, ex=prop, D=(index>0.5))
      est    = c( ipw(dat), sipw(dat),  ipw.zzz(dat),   elw(dat)) 
      result = rbind(result, est)
    }  

   theta = other[1]
   rho   = other[2]    
   model  = other[3]     
   ##############################################################
   ############          boxplot     ############################
   rho1 = 10*rho  
   txt = paste('Pivotal-n=', n, '-rho=', rho1, '-Model=', model, sep='')
   flnm = paste('d:/', txt, '.pdf',   sep='')
   pdf(flnm) 
   boxplot(result[, -1], xaxt='n', xlab='',
               main=paste('Pivotal sampling: n=', n, ', rho=', rho, ', Model=', model,  sep='')
           )
   abline(h=theta) 
   xnm = c('IPW','SIPW', 'ZZZ',   'ELW')
   axis(1, 1:3, labels=xnm[-1])
   dev.off()
   ##############################################################

   rmse  = round( sqrt( N*apply((result-theta)^2, 2, mean)), 2) 
   colnm = c('IPW','SIPW', 'ZZZ',   'ELW')
   rmse = matrix(rmse, 1)
   colnames(rmse) <-  colnm 
   print(rmse)
   list(result=result, rmse =rmse )
}


#####   PPS sampling
simu.pps<-function(nrep, popu, prop, n, other)
{  result=NULL 
   N = length(popu) 
   prop = prop/sum(prop) 
   for( i in 1:nrep){  
      index = rmultinom(1, n, prop) 
      dat=list(y=popu, ex=n*prop, D=index)
      est    = c( ipw(dat), sipw(dat),  ipw.zzz(dat),  elw(dat))  
      result = rbind(result, est)
    }  
   theta = other[1]
   rho   = other[2]    
   model  = other[3]   
   ##############################################################
   ############          boxplot     ############################
   rho1 = 10*rho  
   txt = paste('PPS-n=', n, '-rho=', rho1, '-Model=', model, sep='')
   flnm = paste('d:/', txt, '.pdf',   sep='')
   pdf(flnm) 
   boxplot(result[, -1], xaxt='n', xlab='',
               main=paste('PPS sampling:  n=', n,  ', rho=', rho, ', Model=', model,  sep='')
           )
   abline(h=theta) 
   xnm = c('IPW','SIPW', 'ZZZ',   'ELW')
   axis(1, 1:3, labels=xnm[-1])
   dev.off()
   ##############################################################

   rmse  = round( sqrt(N*apply(  (result-theta)^2, 2, mean)) , 2) 
   colnm = c('IPW','SIPW', 'ZZZ',  'ELW')
   rmse = matrix(rmse, 1)
   colnames(rmse) <-  colnm
   list(result=result, rmse=  rmse )
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
nn    = c(250, 500)

rho.all=c(0.2, 0.8)   
rmse1=rmse2=rmse3=NULL  

for(i in  1:2){ 
  n = nn[i]
  prop = n*prop0
  for(k in 1:2){
     rho = rho.all[k]
     for(j in 1:4){
        model = j
        reg  = (model==1)*sqrt(3)*rho* x  +  (model==2)*sqrt(3)*rho*(x+x^2) + 
                 (model==3)*(sqrt(3)*rho* x  + 5 ) + (model==4)*(sqrt(3)*rho*(x+x^2) + 5)

        popu =  reg +  sqrt(3)*sqrt(1- rho^2)*abs(e)  
        theta= mean(popu)  
        other = c(theta, rho, model) 

        out1 = simu.poisson(nrep, popu, prop, other)$rmse  
        rmse1 = rbind(rmse1,  c(n, rho, model, out1)) 
        print(out1)

        out2 =simu.pivotal(nrep, popu, prop, other)$rmse 
        rmse2 = rbind(rmse2,  c(n, rho, model, out2))   
        print(out2)

        out3 = simu.pps(nrep, popu, prop, n, other)$rmse 
        rmse3 = rbind(rmse3,  c(n, rho, model, out3))  
        print(out3)
    } 
  }
}


rmse.poisson = round(rmse1, 2)
rmse.pivotal = round(rmse2, 2)
rmse.pps     = round(rmse3, 2)


colnm = c('n', 'rho', 'model','IPW','SIPW', 'ZZZ',   'ELW')
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


> rmse.poisson
        n rho model   IPW  SIPW   ZZZ  ELW
 [1,] 250 0.2     1 12.62  5.95  8.63 5.45
 [2,] 250 0.2     2 12.06  6.93  9.51 5.91
 [3,] 250 0.2     3 42.33  5.78 31.78 5.48
 [4,] 250 0.2     4 42.93  6.93 32.59 5.91
 [5,] 250 0.8     1  9.45  7.33  8.27 4.93
 [6,] 250 0.8     2 14.29 16.84 13.55 9.92
 [7,] 250 0.8     3 40.46  7.33 31.53 4.93
 [8,] 250 0.8     4 43.52 16.84 35.42 9.92
 [9,] 500 0.2     1  9.04  4.34  6.28 3.93
[10,] 500 0.2     2 11.38  5.36  6.95 4.16
[11,] 500 0.2     3 36.21  4.59 23.56 3.89
[12,] 500 0.2     4 36.55  5.36 24.06 4.16
[13,] 500 0.8     1  8.07  5.66  5.85 3.39
[14,] 500 0.8     2 10.79 12.84  9.24 6.56
[15,] 500 0.8     3 33.58  5.66 23.09 3.39
[16,] 500 0.8     4 35.35 12.84 25.53 6.56
> rmse.pivotal
        n rho model    IPW  SIPW   ZZZ  ELW
 [1,] 250 0.2     1  10.73  5.87  6.60 5.32
 [2,] 250 0.2     2  10.58  6.94  6.38 5.56
 [3,] 250 0.2     3  38.14  5.87 22.46 5.32
 [4,] 250 0.2     4  37.92  6.94 22.12 5.56
 [5,] 250 0.8     1   6.57  7.26  4.04 4.23
 [6,] 250 0.8     2   6.27 16.63  3.63 7.62
 [7,] 250 0.8     3  34.47  7.26 20.38 4.23
 [8,] 250 0.8     4  33.65 16.63 19.07 7.62
 [9,] 500 0.2     1  40.67  4.38  5.03 3.78
[10,] 500 0.2     2  40.64  5.36  4.89 3.92
[11,] 500 0.2     3 179.21  4.38 17.48 3.78
[12,] 500 0.2     4 179.18  5.36 17.27 3.92
[13,] 500 0.8     1  24.90  5.88  3.07 2.91
[14,] 500 0.8     2  24.84 13.29  2.79 5.03
[15,] 500 0.8     3 163.52  5.88 15.86 2.91
[16,] 500 0.8     4 163.41 13.29 15.03 5.03
> rmse.pps
        n rho model    IPW  SIPW   ZZZ  ELW
 [1,] 250 0.2     1  80.60  6.02  6.80 5.32
 [2,] 250 0.2     2  80.58  7.23  6.56 5.55
 [3,] 250 0.2     3 355.73  6.02 23.13 5.32
 [4,] 250 0.2     4 355.70  7.23 22.78 5.55
 [5,] 250 0.8     1  49.36  7.75  4.16 4.26
 [6,] 250 0.8     2  49.30 17.74  3.72 7.84
 [7,] 250 0.8     3 324.56  7.75 20.95 4.26
 [8,] 250 0.8     4 324.45 17.74 19.60 7.84
 [9,] 500 0.2     1   7.41  4.39  5.11 3.87
[10,] 500 0.2     2   7.28  5.31  4.94 4.03
[11,] 500 0.2     3  27.46  4.39 18.01 3.87
[12,] 500 0.2     4  27.29  5.31 17.76 4.03
[13,] 500 0.8     1   4.53  5.74  3.13 3.00
[14,] 500 0.8     2   4.27 13.14  2.79 5.28
[15,] 500 0.8     3  24.97  5.74 16.36 3.00
[16,] 500 0.8     4  24.32 13.14 15.40 5.28
> 
> end.time=proc.time()   
> cpu.time = end.time-start.time
> print(cpu.time)
   用户    系统    流逝 
1121.18    5.45 1136.06 
> 
> 

