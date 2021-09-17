## EN.553.632 FA20 Bayesian Statistics
## Final Project: Inverse ECG
## Changxin Lai (clai29@jhmi.edu)
## 12/20/2020

library('MASS')
library('progress')
library('GPfit')

# load data
A <- as.matrix(read.csv('A.csv', header = FALSE))
x <- read.csv('x.csv', header = FALSE)
y <- read.csv('y.csv', header = FALSE)
u <- read.csv('vertices.csv', header = FALSE)

N.x <- nrow(x)
N.y <- nrow(y)

x.mean <- mean(as.matrix(x))
x.var <- var(as.vector(as.matrix(x,ncol=1)))

y.mean <- mean(as.matrix(y))
y.var <- var(as.vector(as.matrix(y,ncol=1)))

x.train <- x[,1:500]
y.train <- y[,1:500]

x.test <- x[,501:2499]
y.test <- y[,501:2499]

# normalize spatial locations
u[,1] <- u[,1]-min(u[,1])
u[,2] <- u[,2]-min(u[,2])
u[,3] <- u[,3]-min(u[,3])

u_scale <- 1.01*max(c(max(u[,1]),max(u[,2]),max(u[,3])))

u <- u/u_scale

# fit Gaussian process

## sample a subset of training data
set.seed(123)
sample_ratio <- 0.2
n.sample <- floor(N.x*sample_ratio)
sample_index <- sample(seq(1, N.x), n.sample)

## estimate Gaussian process parameters
beta <- numeric()
sig2 <- numeric()
for (t in seq(50, 500, 50)){
  x.train.t <- (x-x.mean)[sample_index,t]
  u.train <- u[sample_index,]
  
  GPmodel = GP_fit(u.train,x.train.t,maxit = 200)
  
  beta <- rbind(beta, GPmodel$beta)
  sig2 <- rbind(sig2, GPmodel$sig2)
}

beta_mean <- colMeans(beta)
sig2_mean <- colMeans(sig2)

# calculate covariance matrix for x
logR <- -(10^beta_mean[1])*as.matrix(dist(u[,1]))-
  (10^beta_mean[2])*as.matrix(dist(u[,2]))-
  (10^beta_mean[3])*as.matrix(dist(u[,3]))

SIGMA.x <- sig2_mean*exp(logR)

# estimate covariance for noise
alpha <- 1
beta <- 1

noise.train <- as.matrix(y.train)-as.matrix(A)%*%as.matrix(x.train)
noise.train.flat <- as.vector(noise.train)

alpha_n <- alpha + length(noise.train.flat)/2
beta_n <- as.numeric(beta + t(noise.train.flat)%*%noise.train.flat/2)

# estimate x on test set
test_length <- 1000

## estimate x using covariance from Gaussian process
x.low.post <- numeric()
x.mean.post <- numeric()
x.high.post <- numeric()
pb <- progress_bar$new(
  format = "  progress [:bar] :current/:total  eta: :eta",
  total = test_length, clear = FALSE, width= 60)
sigma2.n <- beta_n/(alpha_n-1)
for (t in 501:(501+test_length-1)){
  y.t <- y[,t]
  
  SIGMA.x.post <- t(A)%*%A/sigma2.n+solve(SIGMA.x) 
  MU.x.post <- solve(SIGMA.x.post)%*%(t(A)%*%as.matrix(y.t))/sigma2.n
  
  x.low.post.t <- qnorm(0.025,MU.x.post,diag(SIGMA.x.post))
  x.mean.post.t <- MU.x.post
  x.high.post.t <- qnorm(0.9755,MU.x.post,diag(SIGMA.x.post))
  
  x.low.post <- cbind(x.low.post,x.low.post.t)
  x.mean.post <- cbind(x.mean.post,x.mean.post.t)
  x.high.post <- cbind(x.high.post,x.high.post.t)
  
  pb$tick()
  Sys.sleep(1/100)
}
## calculate correlation coefficients
norm.x <- sqrt(rowSums(as.matrix(x.test[,1:test_length])^2))
norm.x.hat <- sqrt(rowSums(x.mean.post^2))
CC <- rowSums(as.matrix(x.test)[,1:test_length]*x.mean.post)/norm.x/norm.x.hat

## estimate x using empirical covariance
x.low.post.emp <- numeric()
x.mean.post.emp <- numeric()
x.high.post.emp <- numeric()
pb <- progress_bar$new(
  format = "  progress [:bar] :current/:total  eta: :eta",
  total = test_length, clear = FALSE, width= 60)
SIGMA.x.empirical <- 1/10*as.matrix(x.train[,seq(50, 500, 50)]-rowMeans(x.train[,seq(50, 500, 50)]))%*%
  t(as.matrix(x.train[,seq(50, 500, 50)]-rowMeans(x.train[,seq(50, 500, 50)])))
for (t in 501:(501+test_length-1)){
  y.t <- y[,t]
  
  MU.x.post <- SIGMA.x.empirical%*%t(A)%*%solve(A%*%SIGMA.x.empirical%*%t(A)+sigma2.n*diag(dim(A)[1]))%*%as.matrix(y.t)
  
  x.mean.post.t <- MU.x.post
  x.mean.post.emp <- cbind(x.mean.post.emp,x.mean.post.t)

  pb$tick()
  Sys.sleep(1/100)
}
## calculate correlation coefficients
norm.x <- sqrt(rowSums(as.matrix(x.test[,1:test_length])^2))
norm.x.hat <- sqrt(rowSums(x.mean.post.emp^2))
CC.emp <- rowSums(as.matrix(x.test)[,1:test_length]*x.mean.post.emp)/norm.x/norm.x.hat

# plot CC distribution
plot(density(CC),lwd=2,col='red',xlim=c(-0.3,1),
     xlab="Correlataion Coefficient",main=' ')
lines(density(CC.emp),lwd=2,col='blue')
legend('topleft',legend=c("CC using empirical variance",
                          "CC using Gaussian process variance"),
       col=c('blue','red'),lty=1)

# plot signals
standardize <- function(x){
  x.mean=mean(x)
  x.var=var(x)
  return((x-x.mean)/sqrt(x.var))
}
location = 201
plot(seq(1:500),standardize(as.numeric(x.test[location,1:500])),'l',
     ylim=c(-3,3),
     xlab='Time',ylab='Potential (Standarized)',col='black')
lines(seq(1:500),standardize(x.mean.post.emp[location,1:500]),'l',col='blue')
lines(seq(1:500),standardize(x.mean.post[location,1:500]),'l',col='red')
legend("topright",legend=c("True signal",
                       "Reconstructed signal using empirical variance",
                       "Reconstructed signal using Gaussian process variance"),
       col=c('black','blue','red'),lty=1,cex=0.8)



###############################################################################

# estimation using MCMC
run_MCMC = FALSE

if (run_MCMC){
  test_length = 1000
  x.mean.post.MCMC <- numeric()
  pb <- progress_bar$new(
    format = "  progress [:bar] :current/:total  eta: :eta",
    total = test_length, clear = FALSE, width= 60)
  n.MC <- 500
  for (t in 501:(501+test_length-1)){
    y.t <- y[,t]
    
    x.sample <- numeric()
    sigma2.n.sample <- numeric()
    sigma2.n.i <- 1/rgamma(1, alpha_n, beta_n)
    for (i in 1:n.MC){
      SIGMA.x.post <- t(A)%*%A/sigma2.n.i+solve(SIGMA.x) 
      MU.x.post <- solve(SIGMA.x.post)%*%(t(A)%*%as.matrix(y.t))/sigma2.n.i
      
      x.sample.i <- mvrnorm(1, MU.x.post, SIGMA.x.post)
      x.sample <- rbind(x.sample, x.sample.i)
      
      alpha_n.post <- alpha_n + N.y/2
      beta_n.post <- beta_n + t(y.t-A%*%x.sample.i)%*%(y.t-A%*%x.sample.i)/2
      
      sigma2.n.i <- rgamma(1, alpha_n.post, beta_n.post)
      sigma2.n.sample <- rbind(sigma2.n.sample, sigma2.n.i)
    }
    x.mean.post.MCMC <- rbind(x.mean.post.MCMC, rowMeans(x.sample))
    pb$tick()
    Sys.sleep(1/100)
  }
  
  norm.x <- sqrt(rowSums(as.matrix(x.test[,1:test_length])^2))
  norm.x.hat <- sqrt(rowSums(x.mean.post.MCMC^2))
  CC.MCMC <- rowSums(as.matrix(x.test)[,1:test_length]*x.mean.post.MCMC)/norm.x/norm.x.hat
  
  plot(density(CC.MCMC),lwd=2,
       xlab="Correlataion Coefficient",main='MCMC results')
  plot(seq(1:test_length),standardize(as.numeric(x.mean.post.MCMC[location,])),'l',
       ylim=c(-3,3),
       xlab='Time',ylab='Potential (Standarized)',col='black',
       main='MCMC results')
}


