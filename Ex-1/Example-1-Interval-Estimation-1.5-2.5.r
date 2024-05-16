 

#####################   IPW estimator    #####################################
ipw <- function(dat)
{  
   N=length(dat$y) 
    wt = dat$D/dat$ex
    the  = sum(wt*dat$y)/N
    gc = dat$y*wt - the 
    Sig  = sum(gc^2)/N
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


 

####  Robust IPW estimator of    Ma and Wang, 2019, JASA   ##################### 
 
ipw.mw<-function(dat)
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


  if (sum(dat$D==1 & dat$ex >= hn) <= 2) {
      Bnb1 = Bnb2 = 0
    } else {
       Bnb1 =  - mean( (xx.all%*%bet)*(dat$ex<bn1) )
       Bnb2 =  - mean( (xx.all%*%bet)*(dat$ex<bn2) ) 
    }

  ipw.trim1 =  mean(yy*(pw >= bn1)/pw) *n/N            
  theta.mw1 = ipw.trim1-Bnb1   
  tmp1  = sum( ( yy*(pw >= bn1)/pw - theta.mw1)^2)
  Sn.mw1 = (tmp1 +   (N-n)*theta.mw1^2)/(N-1) 

  ipw.trim2 =  mean(yy*(pw >= bn2)/pw)*n/N
  theta.mw2 = ipw.trim2 -Bnb2
  tmp2  = sum( ( yy*(pw >= bn2)/pw - theta.mw2)^2)
  Sn.mw2 = (tmp2 +   (N-n)*theta.mw2^2)/(N-1)  

  c(theta.mw1, theta.mw2, Sn.mw1, Sn.mw2)
}

##############################################################################
quan.low<-function(x){  quantile(x, 0.025) } 
quan.up<-function(x){  quantile(x, 0.975) } 


##############################################################################
crit.cal<-function(dat, theta.est, B)
{
   N = length(dat$ex)
   M = round(N/log(N))  
   stat.all = NULL

   for(i in 1:B){
       ind = sample.int(N, M, replace=FALSE)
       y.s  = dat$y[ind]
       ex.s = dat$ex[ind]
       D.s  = dat$D[ind]
       dat.s= list(y=y.s, ex=ex.s, D=D.s) 

       out.ipw=ipw(dat.s)
       stat.ipw  = sqrt(M)*(out.ipw[1]-theta.est[1])/sqrt(out.ipw[2])  

       out.sipw=sipw(dat.s)
       stat.sipw  = sqrt(M)*(out.sipw[1]-theta.est[2])/sqrt(out.sipw[2])  
 
       out.mw = ipw.mw(dat.s) 
       stat.mw1 = sqrt(M)*(out.mw[1]-theta.est[3])/sqrt(out.mw[3])
       stat.mw2 = sqrt(M)*(out.mw[2]-theta.est[4])/sqrt(out.mw[4]) 

       out.elw = elw(dat.s)    
       stat.elw  = sqrt(M)*(out.elw[1]-theta.est[5])/sqrt(out.elw[2])  

       stat.all=rbind(stat.all, c(stat.ipw, stat.sipw, stat.mw1, stat.mw2, stat.elw))
   }


  stat.ipw  = stat.all[, 1]
  mean.ipw  = mean(stat.ipw)
  qu.ipw = quantile(abs(stat.ipw-mean.ipw), 0.95)
  low.ipw  = mean.ipw  - qu.ipw
  up.ipw   = mean.ipw  + qu.ipw

  stat.sipw  = stat.all[, 2]
  mean.sipw  = mean(stat.sipw)
  qu.sipw = quantile(abs(stat.sipw-mean.sipw), 0.95)
  low.sipw  = mean.sipw  - qu.sipw
  up.sipw   = mean.sipw  + qu.sipw


  stat.elw  = stat.all[, 5]
  mean.elw  = mean(stat.elw)
  qu.elw = quantile(abs(stat.elw-mean.elw), 0.95)
  low.elw   = mean.elw  - qu.elw
  up.elw    = mean.elw  + qu.elw

   

  crit.low = c(low.ipw, low.sipw, apply(stat.all[, 3:4], 2, quan.low), low.elw )
  crit.up  = c(up.ipw, up.sipw, apply(stat.all[, 3:4], 2, quan.up),  up.elw) 
  c(crit.low, crit.up) 
}

 


################################################################################
simu.cov<-function(nrep, gam, N, fun_index, coef, my, B)
{
   result=NULL
   
   for( i in 1:nrep){
     # set.seed(i*10+5)
      dat    = data.gen(gam, N, fun_index, coef) 


      out.ipw=ipw(dat)
      theta.ipw = out.ipw[1]
      stat.ipw  = sqrt(N)*(theta.ipw - my)/sqrt(out.ipw[2])        
      st.ipw    = sqrt(N)*abs(theta.ipw - my)/sqrt(out.ipw[2]) 

      out.sipw=sipw(dat)
      theta.sipw = out.sipw[1]
      stat.sipw  = sqrt(N)*(theta.sipw - my)/sqrt(out.sipw[2])        
      st.sipw    = sqrt(N)*abs(theta.sipw - my)/sqrt(out.sipw[2]) 

      out.sipw=sipw(dat)
      theta.sipw = out.sipw[1]
      stat.sipw  = sqrt(N)*(theta.sipw - my)/sqrt(out.sipw[2])        
      st.sipw    = sqrt(N)*abs(theta.sipw - my)/sqrt(out.sipw[2]) 

      out.mw = ipw.mw(dat)
      theta.mw1 = out.mw[1]
      theta.mw2 = out.mw[2]  

      stat.mw1  = sqrt(N)*(theta.mw1 - my)/sqrt(out.mw[3])
      stat.mw2  = sqrt(N)*(theta.mw2 - my)/sqrt(out.mw[4]) 

      out.elw   = elw(dat) 
      theta.elw = out.elw[1]  
      stat.elw  = sqrt(N)*(theta.elw-my) /sqrt(out.elw[2])  

      theta.est = c(theta.ipw, theta.sipw,theta.mw1, theta.mw2, theta.elw)
      stat = c(stat.ipw, stat.sipw, stat.mw1, stat.mw2, stat.elw ) 
  
      crit  = crit.cal(dat, theta.est, B)
      index = (stat-crit[1:5] > 0 )*(stat-crit[6:10] < 0 )
      
      st.elw = sqrt(N)*abs(theta.elw  - my)/sqrt(out.elw[2])
      index=c(st.ipw<1.96, index[1],st.sipw<1.96, index[2:4],  st.elw<1.96, index[5])      

      wd.ipw  = 2*1.96*sqrt(out.ipw[2])/sqrt(N)
      wd.sipw = 2*1.96*sqrt(out.sipw[2])/sqrt(N)
      wd.elw  = 2*1.96*sqrt(out.elw[2]) /sqrt(N)
      tmp = (crit[6:10]-crit[1:5])*1.96*sqrt(c(out.ipw[2], out.sipw[2], out.mw[3],out.mw[4], out.elw[2])) /sqrt(N)

      width = c(wd.ipw, tmp[1], wd.sipw, tmp[2:4], wd.elw, tmp[5])
#      width =  c(wd1, (crit[5:8]-crit[1:4])*tmp, wd2) 

      result = rbind(result, c(index*100, width))
      print( round(c(gam, coef, fun_index, i,  stat, crit), 2))
    }
  
   cvrg = round( apply(result, 2, mean), 4)


   out = c(N, gam,  coef, fun_index,  cvrg )
   out = matrix(out, 1)
   colnames(out)=c('N', 'gamma','c', 'Model', 'IPW.an','IPW.re', 'SIPW.an','SIPW.re', 'MW1', 'MW2', 
                 'ELW.an', 'ELW.re', 'IPW.an','IPW.re', 'SIPW.an','SIPW.re', 'MW1', 'MW2', 'ELW.an', 'ELW.re')
   return(out)
}
 



###############################################################################
####    Generate data 
data.gen<- function(gam, N, fun_index, coef)
{
   ex  = (runif(N))^(1/(gam-1))
   eta = (rchisq(N, 4)-4)/sqrt(8)  
   reg = (fun_index==1)*cos(2*pi*ex) + (fun_index==2)*(1-ex) + 
           (fun_index==3)*(5+cos(2*pi*ex))+ (fun_index==4)*(1-ex+5)
   y =  reg + eta*coef    
   D   = rbinom(N, 1, ex)
   dat = list(y=y, ex=ex, D=D)
   return(dat)
}


################################################################################
start.time=proc.time() 
library(MASS)
eeps=sqrt(.Machine$double.eps)    #### root of machine precision

nrep=5000
NN=c(2000)
gamm=c( 1.5,  2.5)
cff = c(1,  0.1) 
B = 1000

cov.all = NULL
for(i in 1:1){
   N = NN[i]
   for(k in 1:2){
      gam = gamm[k] 
      for(a in 1:2){
        coef = cff[a]
        for(j in 1:4){  
            fun_index = j 
            dat0 = data.gen(gam, 200000, fun_index, coef)
            my   = mean(dat0$y) 
            cvrg = simu.cov(nrep, gam, N, fun_index, coef, my, B)
            cov.all = rbind(cov.all, cvrg)
            print(cov.all)
         }
      }
   }
}
 
cov.all

end.time=proc.time()
cpu.time = end.time-start.time
print(cpu.time)

write.table(cov.all, 'e:/cov-all-ex1.txt')
 

 
###############################################################################

> cov.all
         N gamma   c Model IPW.an IPW.re SIPW.an SIPW.re   MW1   MW2 ELW.an ELW.re IPW.an    IPW.re SIPW.an SIPW.re    MW1    MW2 ELW.an ELW.re
 [1,] 2000   1.5 1.0     1  76.14  84.52   78.42   88.86 89.44 92.86  82.58  91.48 0.5561    3.0290  0.4363  1.4553 1.0606 0.7803 0.3297 1.0441
 [2,] 2000   1.5 1.0     2  78.24  87.24   85.08   91.14 89.68 94.52  81.62  91.04 1.5935 6397.7332  0.3669  1.1593 1.1658 0.8313 0.3124 1.0475
 [3,] 2000   1.5 1.0     3  76.76  86.40   77.84   88.94 82.74 81.44  82.02  90.94 2.9483   18.4471  0.4469  1.5755 3.3253 1.5070 0.3338 1.0618
 [4,] 2000   1.5 1.0     4  78.74  86.32   85.36   91.60 82.32 83.32  82.54  92.38 3.2826  121.8590  0.3700  1.1601 3.2768 1.5515 0.3124 1.0465
 [5,] 2000   1.5 0.1     1  76.94  88.86   79.58   91.54 79.76 68.98  92.86  95.72 0.4914    6.3610  0.3040  1.0316 0.5435 0.2612 0.1260 0.2874
 [6,] 2000   1.5 0.1     2  76.22  87.68   76.54   89.20 81.92 82.40  87.62  91.72 0.4777    5.0524  0.1302  0.4369 0.5522 0.2243 0.0497 0.1223
 [7,] 2000   1.5 0.1     3  77.40  86.18   78.00   90.44 80.50 78.90  91.84  94.84 2.9527   22.3009  0.3048  1.0362 3.1775 1.2930 0.1259 0.2876
 [8,] 2000   1.5 0.1     4  77.88  85.96   76.22   89.00 81.76 79.12  85.38  89.28 2.9585   31.3594  0.1310  0.4412 3.1254 1.3613 0.0499 0.1227
 [9,] 2000   2.5 1.0     1  94.02  91.40   93.60   92.46 91.92 93.14  93.32  93.20 0.1782    0.3415  0.1778  0.3513 0.3231 0.3232 0.1694 0.3437
[10,] 2000   2.5 1.0     2  93.52  93.48   93.72   93.54 92.38 93.40  93.72  93.60 0.1694    0.3466  0.1507  0.3039 0.3248 0.3212 0.1417 0.2870
[11,] 2000   2.5 1.0     3  94.14  93.68   93.86   92.74 89.40 88.16  93.56  93.40 0.6303    1.2613  0.1775  0.3487 1.0294 0.8376 0.1695 0.3441
[12,] 2000   2.5 1.0     4  94.26  93.34   93.80   94.06 90.24 88.96  93.46  93.94 0.6821    1.3462  0.1501  0.2991 1.1180 0.9486 0.1416 0.2869
[13,] 2000   2.5 0.1     1  95.18  92.88   94.94   93.92 91.18 90.90  94.04  94.72 0.1075    0.2036  0.1078  0.2139 0.1769 0.1575 0.0989 0.2080
[14,] 2000   2.5 0.1     2  93.76  94.42   94.02   94.90 89.40 88.28  94.50  94.04 0.0953    0.2080  0.0549  0.1192 0.1592 0.1272 0.0354 0.0711
[15,] 2000   2.5 0.1     3  94.42  93.92   94.64   93.48 89.36 87.96  94.02  94.80 0.6168    1.2319  0.1083  0.2144 0.9867 0.7841 0.0990 0.2083
[16,] 2000   2.5 0.1     4  94.08  92.96   94.34   95.08 89.90 88.70  94.38  94.32 0.6646    1.3062  0.0546  0.1177 1.0813 0.9017 0.0353 0.0709
> 
> end.time=proc.time()
> cpu.time = end.time-start.time
> print(cpu.time)
    用户     系统     流逝 
57247.02     7.88 57299.64 

