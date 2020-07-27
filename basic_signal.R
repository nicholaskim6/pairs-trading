

closeAllConnections()
rm(list=ls())
library(rkdb)
h1 <- open_connection('orf474', 6000)
h2 <- open_connection('orf474', 7000)

date = '2017.10.13'
sym1="GEZ8"
sym2="GEH9"

q1=paste0("([] time:0D00:00:00+10000000000*til 360*15)")
q2=paste0("select mid:(last bid+last ask)%2 by time from quote where
          date=",date,",sym=`",sym1)
q3=paste0("select mid:(last bid+last ask)%2 by time from quote where
          date=",date,",sym=`",sym2)
q4=paste0("aj[`time;",q1,";",q2,"]")
q5=paste0("aj[`time;",q1,";",q3,"]")

dt1=execute(h1,q4)
dt2=execute(h1,q5)c

n=nrow(dt1)

mean1=c(dt1$mid[1])
mean2=c(dt2$mid[1])
std1=c(0)
std2=c(0)
cov=c(dt1$mid[1]*dt2$mid[1])
corr=c(NA)

for (i in 2:nrow(dt1)){
  mean1[i]=(mean1[i-1]*(i-1)+dt1$mid[i])/i
  mean2[i]=(mean2[i-1]*(i-1)+dt2$mid[i])/i
  
  std1[i]=sqrt((std1[i-1]**2*(i-1)+(i-1)*mean1[i-1]**2+dt1$mid[i]**2
                -i*mean1[i]**2 )/i)
  std2[i]=sqrt((std2[i-1]**2*(i-1)+(i-1)*mean2[i-1]**2+dt2$mid[i]**2
                -i*mean2[i]**2 )/i)
  
  cov[i]=cov[i-1]+dt1$mid[i]*dt2$mid[i]
  corr[i]=(cov[i]/i-mean1[i]*mean2[i])/(std1[i]*std2[i])
}

start_index=360*7+1
dt1=dt1[start_index:n,]
dt2=dt2[start_index:n,]
mean1=mean1[start_index:n]
mean2=mean2[start_index:n]
std1=std1[start_index:n]
std2=std2[start_index:n]
cov=cov[start_index:n]
corr=corr[start_index:n]


V=cbind(c(1,1),c(1,-1))/sqrt(2)
x_y=rbind((dt1$mid-mean1)/std1,(dt2$mid-mean2)/std2)

x_y_star=matrix(ncol = 2,nrow=nrow(dt1))
for (i in 1:nrow(dt1)){
  xi_eta=t(V)%*%c((dt1$mid[i]-mean1[i])/std1[i],(dt2$mid[i]-mean2[i])/std2[i])
  xi_eta[2]=0
  x_y_star[i,]=diag(c(std1[i],std2[i]))%*%V%*%xi_eta+c(mean1[i],mean2[i])
}

signal_x=x_y_star[,1]-dt1$mid
signal_y=x_y_star[,2]-dt2$mid

ticksize=execute(h1, "select minpxincr by inst from instinfo where
                 inst=`GE")$minpxincr[[1]]

lag=10
x_diff=diff(dt1$mid,lag=lag)/ticksize
y_diff=diff(dt2$mid,lag=lag)/ticksize
x_signal=head(signal_x,-lag)/ticksize
y_signal=head(signal_y,-lag)/ticksize

cor(x_diff,x_signal)
cor(y_diff,y_signal)

lm1=lm(x_diff~x_signal)
summary(lm1)

lm2=lm(y_diff~y_signal)
summary(lm2)
  
plot(x_signal,x_diff, pch=5, main = "GEZ8/GEH9 \n 2017-10-13 / 7:00-15:00 / 10 sec / 120 sec fwd", col="darkgreen",xlim=c(-0.6,0.6),ylim=c(-2,2),xlab="Signal in Ticks", ylab="Forward return in ticks")
points(y_signal,y_diff,col='orange', pch=1)
legend(-0.55, 1.80, pch=c(5,1), legend=c("GEZ8: 0.488 ± 0.0600", "GEH9: 0.758 ± 0.057"),
       col=c("darkgreen", "orange"), cex=0.8)
abline(lm1,col="darkgreen")
abline(lm2,col="orange")
