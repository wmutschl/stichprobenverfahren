###########################################
#### Aufgabe Einschlusswahrscheinlichkeiten
###########################################
# Matrix mit Einschlusswahrscheinlichkeiten
Iks <- function(x,y) as.numeric(is.element(x,y))	
N <- 4
n <- 2

# a)
S <- combn(1:N,n)
M <- choose(N,n)	
ps <- rep(1/M,M)
ind <- apply(S,2,function(z) Iks(1:N,z)); ind
pi_k = colSums(t(ind)*ps);
round(pi_k,2)
sum(pi_k)

#c)
M <- 3
S <- cbind(c(1,3),c(1,4),c(2,4))
ps <- c(0.1,0.6,0.3)
ind <- apply(S,2,function(z) Iks(1:N,z)); ind
pi_k <- colSums(t(ind)*ps)
round(pi_k,2)
sum(pi_k)

#inefficient solution
pi_kl <- matrix(NA,N,N)
for (k in 1:N){
  for (l in 1:N) {
    pi_kl[k,l] <- sum(apply(S,2,function(z) Iks(k,z)*Iks(l,z))*ps)
  }
}

#d)
Delta_kl <- matrix(NA,N,N)
for (k in 1:N){
  for (l in 1:N) {
    Delta_kl[k,l] <- pi_kl[k,l] - pi_kl[k,k]*pi_kl[l,l]
  }
}

##################################################################
#### Aufgabe Schätzung mithilfe von Einschlusswahrscheinlichkeiten
##################################################################
# Matrix mit Einschlusswahrscheinlichkeiten
Iks <- function(x,y) as.numeric(is.element(x,y))
Y <- c(1,2,5,12,30)
N <- length(Y)
n <- 3
ps <- 1:M/sum(1:M); round(ps,3)
#a)
M <- choose(N,n);M
#b)
S <- combn(N,n);S
ind <- apply(S,2,function(z) Iks(1:N,z)); ind
#c)
pi_k <- colSums(t(ind)*ps);round(pi_k,3)
pi_kl <- matrix(NA,N,N)
for (k in 1:N){
  for (l in 1:N){
    pi_kl[k,l] <- sum(apply(S,2,function(z) Iks(k,z)*Iks(l,z))*ps)
  }
}
round(pi_kl,3)
#d)
Delta_kl <- matrix(NA,N,N)
for (k in 1:N){
  for (l in 1:N) {
    Delta_kl[k,l] <- pi_kl[k,l] - pi_kl[k,k]*pi_kl[l,l]
  }
}
round(Delta_kl,2)
#e)
ybar.hat <- 1/N*apply(S,2,function(z) sum(Y[z]/pi_k[z]))
round(ybar.hat,2)
#f)
mean(Y) #wahrer Wert
sum(ybar.hat*ps) #unverzerrter Schätzer ergibt wahren Wert
#g)
# Funktion die Varianz des Horvitz-Thompson Schätzers für jede Stichprobe schätzt
vhatHT <- function(s){
  n <- length(s)
  sl <- rep(NA,n)
  sk <- sl
  for (j1 in 1:n){
    k <- s[j1]
    for (j2 in 1:n) {
      l <- s[j2]
      sl[j2] <- 1/pi_kl[k,l]*(pi_kl[k,l]/(pi_k[k]*pi_k[l])-1)*Y[k]*Y[l]
    }
    sk[j1] <- sum(sl)
  }
  sum(sk)/N^2
}
vHT <- apply(S,2,vhatHT)
round(vHT,2)
# Erster Wert ist negativ! Dies kann passieren beim Varianz Schätzer von Horvitz-Thompson

# Alternativ Yates-Grundi Schätzer
vhatYG <- function(s){
  n <- length(s)
  sl <- rep(NA,n)
  sk <- sl
  for (j1 in 1:n){
    k <- s[j1]
    for (j2 in 1:n) {
      l <- s[j2]
      sl[j2] <- Delta_kl[k,l]/pi_kl[k,l]*(Y[k]/pi_k[k]-Y[l]/pi_k[l])^2
    }
    sk[j1] <- sum(sl)
  }
  sum(sk)*(-1)/(2*N^2)
}
vYG <- apply(S,2,vhatYG)
round(vYG,2)

#h)
sum((ybar.hat-mean(Y))^2*ps) # wahrer Wert der Varianz des Schätzers
sum(vHT*ps)
sum(vYG*ps)


########
psid <- read.csv2("psid.csv")
N <- nrow(psid)
n <- 20
Y <- psid$wage
f <- n/N
Ybar <- mean(Y);Ybar
V.Ybar <- (1-f)/n*var(Y); V.Ybar

f_var <- function(y,N){
  n <- length(y)
  f <- n/N
  return((1-f)/n*var(y))
}

B <- 100000
m <- rep(NA,B)
v <- rep(NA,B)
for (b in 1:B) {
  y <- sample(Y,n)
  m[b] <- mean(y)
  v[b] <- f_var(y,N)
}

vm <- v/10^6 
V.Ybarm <- V.Ybar/10^6

plot(density(m),lwd=2,xlab='Schätzwert',main='')
arrows(Ybar,5e-06,Ybar,0,length=0.1,angle=25)
text(Ybar,7.5e-06,expression(bar(y)[U]))

plot(density(vm),lwd=2,xlab='Schätzwert in Millionen',main='')
arrows(V.Ybarm+1000,0.002,V.Ybarm,0,length=0.1,angle=25)
text(V.Ybarm+1000,0.0026,expression(V(bar(y)[U])))

psid <- read.csv2('psid.csv')
N <- nrow(psid)
B <- 10000
n <- 30
f <- n/N
e <- rep(NA,B)
v <- rep(NA,B)
Y <- psid$sector==7
for (i in 1:B){
  y <- sample(Y,n)
  e[i] <- mean(y)
  v[i] <- (1-f)/n*var(y)
}
plot(density(e))
# annähernd normal verteilt
plot(density(v))
# linksschief, definitiv nicht normalverteilt



library(pps)
B <- 10000; N <- 5; n <- 3
e_sample <- matrix(NA,B,n)
e_sampford <- matrix(NA,B,n)
p <- 4:8/sum(4:8);p

for (i in 1:B){
  e_sample[i,] <- sample (1:N,n,prob=p)
  e_sampford[i,] <- sampford(p,n)
}
pi_emp_sample <- rep(NA,N)
pi_emp_sampford <- rep(NA,N)

for (i in 1:N){
  pi_emp_sample[i] <- sum(apply(e_sample,1, function(z) i%in%z))
  pi_emp_sampford[i] <- sum(apply(e_sampford,1, function(z) i%in%z))
}

rbind(p*n,round(pi_emp_sample/B,3),round(pi_emp_sampford/B,3))




library(samplingbook)
data(influenza)
summary(influenza)

# 1) Usage of pps.sampling
set.seed(123)
pps <- pps.sampling(z=influenza$population,n=20,method='sampford')
pps
sample <- influenza[pps$sample,]
sample











# 2) Usage of htestimate
set.seed(123)
pps <- pps.sampling(z=influenza$population,n=20,method='midzuno')
sample <- influenza[pps$sample,]
# htestimate()
N <- nrow(influenza)
# exact variance estimate
PI <- pps$PI
htestimate(sample$cases, N=N, PI=PI, method='ht')
htestimate(sample$cases, N=N, PI=PI, method='yg')
# approximate variance estimate
pk <- pps$pik[pps$sample]
htestimate(sample$cases, N=N, pk=pk, method='hh')
pik <- pps$pik
htestimate(sample$cases, N=N, pk=pk, pik=pik, method='ha')
# without pik just approximative calculation of Hajek method
htestimate(sample$cases, N=N, pk=pk, method='ha') 
# calculate confidence interval based on normal distribution for number of cases
est.ht <- htestimate(sample$cases, N=N, PI=PI, method='ht')
est.ht$mean*N  
lower <- est.ht$mean*N - qnorm(0.975)*N*est.ht$se
upper <- est.ht$mean*N + qnorm(0.975)*N*est.ht$se
c(lower,upper) 
# true number of influenza cases
sum(influenza$cases)


Flaeche <- c(1232,327,1346,1285,428,871,1042,1262,497,1016,651,1170,2630,515,895,1055,2110,979,671,120,541,1331,842,162,206)
Reisflaeche <- c(688,231,768,898,417,697,785,1190,338,745,392,1055,2400,330,810,1026,1666,929,565,101,516,1036,568,137,107)
N <- 892
n <- 25
X.dot <- 568565
sum(Reisflaeche/Flaeche)
y.sum <- X.dot * sum(Reisflaeche/Flaeche)/n
y.sum
var(Reisflaeche/Flaeche)
var.y.sum <- X.dot^2 * var(Reisflaeche/Flaeche)/n
var.y.sum
sqrt(var.y.sum)

lower <- y.sum - sqrt(var.y.sum)*qnorm(0.975)
upper <- y.sum + sqrt(var.y.sum)*qnorm(0.975)
cbind(lower, upper)

## VL
psid <- read.csv2("psid.csv")
Y <- psid$wage
X <- psid$eduyears
cor(Y,X)
Yt <- Y/1000
plot(X,Yt,xlim=c(5,18),xlab = "Ausbildungsjahre",ylab = "Jahreslohn in 1000 USD")
abline(lm(Yt~X),lwd=2)

plot(X,Yt,xlim=c(5,18),ylim=c(0,200),xlab = "Ausbildungsjahre",ylab = "Jahreslohn in 1000 USD")
abline(lm(Yt~X),lwd=2)


##Klausur
Y <- c(1, 2, 4, 3, 5, 7, 6, 8, 9)
Z <- c(1, 1, 1, 2, 2, 2, 3, 3, 3)
N <- length(Y)

mY <- mean(Y)
vY <- var(Y)
n <- 6
M_SI <- choose(N,n)	
vSI <- (1-6/9)/6*7.5

#2
Nh <- 3
M_h <- choose(3,2)^3
nh <- c(2,2,2)
fh <- nh/Nh
Nh <- tapply(Y,Z,length);Nh
Wh <- Nh/N
mYh <- tapply(Y,Z,mean)
vYh <- tapply(Y,Z,var)

1/N^2*sum(Nh^2*(1-fh)/nh*vYh)

# Between strata
Vext <- sum((mYh-mY)^2*Nh)/(N-1);Vext
# Within strata
Vint <- sum(vYh*(Nh-1))/(N-1);Vint
# Fraction explained (%)
round(Vext/vY*100,1)
cor(Y,Z)
Vext/vY
o <- order(Y)
Z[o]


#3
NI <- 3
nI <- 2
fI <- nI/NI

tY <- sum(Y)
si <- combn(3,2);si
htY <- NI/nI*apply(si,2,
                   function(z) sum(Y[Z%in%z]))
htY
mean(htY)
var(htY)*(NI-1)/NI
ti <- tapply(Y,Z,sum);ti
vhtY <- NI^2*(1-fI)/nI*var(ti);vhtY

v_C <- 1/N^2*NI^2*(1-fI)/nI*var(ti)

#4
y1 <- 3; y2 <- 10; y3 <- 7
p1 <- 0.06; p2 <- 0.2;p3<-0.1 
m <-3
pi1 <- 1-(1-p1)^m
pi2 <- 1-(1-p2)^m
pi3 <- 1-(1-p3)^m
t_hh <- 1/3*(y1/p1+y2/p2+y3/p3)

p1 = 0.06
p2=0.2;p3=0.1
y1=3
y2=10
y3=7
pi1 = 1-(1-p1)^3
pi2 = 1-(1-p2)^3
pi3 = 1-(1-p3)^3
pi_12 = pi1+pi2-(1-(1-p1-p2)^3)
pi_13 = pi1+pi3-(1-(1-p1-p3)^3)
pi_23 = pi2+pi3-(1-(1-p2-p3)^3)
pi_12^(-1)*(pi_12/(pi1*pi2)-1)*y1*y2 + pi_13^(-1)*(pi_13/(pi1*pi3)-1)*y1*y3 + pi_23^(-1)*(pi_23/(pi2*pi3)-1)*y2*y3