 

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
{  
   wt = dat$D/dat$ex
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
{ 
  N = length(dat$D)
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
crit.cal<-function(dat, theta.est, B)
{
   quan.low<-function(x){  quantile(x, 0.025) } 
   quan.up<-function(x){  quantile(x, 0.975) } 

   N = length(dat$ex)
   M = round(N/log(N))  
   stat.all = NULL

   for(i in 1:B) {
       ind = sample.int(N, M, replace=FALSE)
       y.s  = dat$y[ind]
       ex.s = dat$ex[ind]
       D.s  = dat$D[ind]
       dat.s= list(y=y.s, ex=ex.s, D=D.s) 

       # `stat` below means Wald statistic
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
  # * Here we use `quantile` to derive the CI
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
gamm=c( 1.3, 1.9)
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

write.table(cov.all, 'e:/cov-all-ex1-39.txt')
 

 
#################################################################################

> cov.all
         N gamma   c Model IPW.an IPW.re SIPW.an SIPW.re   MW1   MW2 ELW.an ELW.re IPW.an    IPW.re SIPW.an SIPW.re    MW1    MW2 ELW.an ELW.re
 [1,] 2000   1.3 1.0     1  60.64  78.06   65.26   85.20 90.28 89.20  69.68  89.12 1.1174  110.7574  0.5597  2.5656 2.0096 1.1866 0.4489 2.1523
 [2,] 2000   1.3 1.0     2  59.68  78.70   76.16   88.30 92.34 94.90  67.86  88.96 0.9791   29.4286  0.4848  2.0223 2.1248 1.2648 0.4360 2.2280
 [3,] 2000   1.3 1.0     3  60.82  78.88   67.84   86.44 77.70 80.60  70.14  88.92 7.2252 5683.8083  0.5586  2.5268 4.8063 1.6642 0.4537 2.1835
 [4,] 2000   1.3 1.0     4  60.16  77.20   75.94   88.90 77.10 82.76  68.14  89.02 6.3866  840.5085  0.4799  1.9510 4.6730 1.6682 0.4387 2.2516
 [5,] 2000   1.3 0.1     1  60.58  82.76   64.90   86.94 74.00 33.14  90.06  93.78 0.9837  195.4618  0.3535  1.4533 0.7718 0.2765 0.1200 0.2983
 [6,] 2000   1.3 0.1     2  59.08  80.20   61.04   84.48 76.04 78.84  75.22  85.90 0.8264   24.0514  0.1466  0.6074 0.7987 0.2278 0.0567 0.1832
 [7,] 2000   1.3 0.1     3  60.24  79.32   63.32   87.66 74.66 72.54  85.80  90.12 5.0666   91.1688  0.3523  1.4512 4.4242 1.2242 0.1198 0.2989
 [8,] 2000   1.3 0.1     4  60.94  77.58   60.78   84.36 75.70 76.06  68.26  77.74 5.7126  216.1039  0.1482  0.6166 4.3178 1.2742 0.0570 0.1860
 [9,] 2000   1.9 1.0     1  89.62  90.56   89.48   92.48 89.74 93.04  90.96  93.86 0.2894    0.9503  0.2726  0.6674 0.5173 0.4679 0.2252 0.5252
[10,] 2000   1.9 1.0     2  91.24  93.82   92.98   94.70 90.64 93.12  91.84  93.76 0.2890    0.7772  0.2269  0.5364 0.5542 0.4905 0.1974 0.4675
[11,] 2000   1.9 1.0     3  90.24  92.88   91.40   93.74 86.40 84.48  91.98  94.48 1.2878    3.6956  0.2773  0.6691 1.8230 1.1731 0.2267 0.5294
[12,] 2000   1.9 1.0     4  90.28  91.82   91.38   93.06 86.88 85.16  90.64  92.80 1.3011    3.3728  0.2256  0.5247 1.8512 1.2607 0.1977 0.4691
[13,] 2000   1.9 0.1     1  88.76  91.16   89.10   93.40 86.30 87.24  92.46  95.20 0.2039    0.5658  0.1897  0.4886 0.2878 0.2013 0.1204 0.2699
[14,] 2000   1.9 0.1     2  89.02  93.34   88.96   93.38 85.58 85.92  92.44  93.84 0.2038    0.6686  0.0892  0.2355 0.2955 0.1815 0.0427 0.0909
[15,] 2000   1.9 0.1     3  89.86  92.64   89.12   93.48 85.24 82.12  92.66  95.18 1.2342    3.3981  0.1884  0.4822 1.7600 1.0800 0.1205 0.2703
[16,] 2000   1.9 0.1     4  90.24  91.76   90.10   94.24 86.30 85.26  93.28  93.42 1.2864    3.1604  0.0897  0.2346 1.7884 1.1774 0.0428 0.0910
> 
> end.time=proc.time()
> cpu.time = end.time-start.time
> print(cpu.time)
    �û�     ϵͳ     ���� 
56317.39     7.06 56366.16 



