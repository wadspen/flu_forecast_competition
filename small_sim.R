#simple simulation of asg to assess our MCMC fitting
library(mixtools)
library(boot)
library(rstan)
library(dplyr)


ASG <- function(x,beta1,beta2,eta,mu,sig12,sig22) {
  
  (beta1 + (eta-beta1)*exp(-(x-mu)^2/(2*sig12)))*(x<mu) +
    (beta2 + (eta-beta2)*exp(-(x-mu)^2/(2*sig22)))*(x>=mu)
}

theta <- c(-5.012532,-5.290287,-3.070383,26.612454,24.012684,128.339209)

n <- c(101088, 113676, 105918, 108208, 106342, 109962, 108314, 99628, 
  102408, 101594, 100834, 99680, 92962, 88028, 93908, 93558, 94218, 
  91056, 92622, 93226, 90628, 82466, 80802, 88742, 86772, 85768, 
  85150, 86964, 86950, 84248, 83294, 80728, 83990, 102576, 81926, 
  84590, 85912, 86244, 87242, 45982, 49003, 47872, 47207, 45833, 
  47714, 46932, 42671, 48312, 45579, 46044, 44987, 41553)


sigbet <- .2
siget <- .1
sigmu <- 2.5
sigsig1 <- 3
sigsig2 <- 25
sigs <- c(sigbet,sigbet,siget,sigmu,sigsig1,sigsig2)

test <- rnorm(5000,sig22,sigsig2)

rasg <- function(x,theta,sigs) {
  pars <- rmvnorm(1,theta,diag(sigs^2))
  y <- inv.logit(ASG(x,pars[1],pars[2],pars[3],pars[4],pars[5],pars[6]))
  return(y)
}

small_sim <- data.frame(x=1:52) %>% 
  group_by(x) %>% 
  summarise(
    y=rasg(x,theta=theta,sigs=sigs)
  ) 

small_sim$n <- n

small_sim <- small_sim %>% 
  group_by(x,y,n) %>%
  summarise(
    c = rbinom(1,n,y)
  )

# theta[5:6] <- theta[5:6]^2

model_prior = "
functions {
   real asg(row_vector param, row_vector betas, real x) {
    real beta1;
    real beta2;
    real eta;
    real mu;
    real sig21;
    real sig22;
    real ASG;
    
    beta1 = param[1];
    beta2 = param[2];
    eta = param[3];
    mu = param[4];
    sig21 = param[5];
    sig22 = param[6];
    
    ASG = ((beta1 + (eta-beta1)*exp(-((x-mu)^2)/(2*sig21^2)))*(x < mu) +
       (beta2 + (eta-beta2)*exp(-((x-mu)^2)/(2*sig22^2)))*(mu <= x));
    return ASG;
  }
}

data {
int<lower=0> N;                                 // number of points
int<lower=0> n[N];                              // patients
int<lower=0> c[N];                              // positive cases
real<lower=0> x[N];                             // week
int<lower=0> n_params;                          // number of parameters ASG


vector[n_params] m0;
vector[n_params,n_params] C0;
row_vector[6] betas;
}
parameters {

row_vector[real real real<lower=0> real<lower=13> real<lower=0> real<lower=0>] theta;

}

model {

theta ~ multi_normal(m0,C0);

for (i in 1:N) c[i] ~ binomial_logit(n[i],asg(theta, betas, x[i])); 

}

//generated quantities {
  //  vector[N] postpred_pr;
    //for (i in 1:N)
      //postpred_pr[i] = binomial_rng(n[i],asg(theta, betas,
       //                                             x[i]));
  //}
"


model_prior = "
functions {
   real asg(real beta1, real beta2, real eta, real mu, real sig1,
   real sig2, real x) {
    
    real ASG;
    
    ASG = ((beta1 + (eta-beta1)*exp(-((x-mu)^2)/(2*sig1^2)))*(x < mu) +
       (beta2 + (eta-beta2)*exp(-((x-mu)^2)/(2*sig2^2)))*(mu <= x));
    return ASG;
  }
}

data {
int<lower=0> N;                                 // number of points
int<lower=0> n[N];                              // patients
int<lower=0> c[N];                              // positive cases
real<lower=0> x[N];                             // week
int<lower=0> n_params;                          // number of parameters ASG


vector[n_params] m0;
vector[n_params] C0;
row_vector[6] betas;
}
parameters {

real beta1;
real beta2;
real eta;
real mu;
real sig1;
real sig2;

}

model {

beta1 ~ normal(m0[1],C0[1]);
beta2 ~ normal(m0[2],C0[2]);
eta ~ normal(m0[3],C0[3]);
mu ~ normal(m0[4],.01);
sig1 ~ normal(m0[5],C0[5]);
sig2 ~ normal(m0[6],C0[6]);

for (i in 1:N) c[i] ~ binomial_logit(n[i],asg(beta1,beta2,eta,mu,sig1,sig2,x[i])); 

}

//generated quantities {
  //  vector[N] postpred_pr;
    //for (i in 1:N)
      //postpred_pr[i] = binomial_rng(n[i],asg(theta, betas,
       //                                             x[i]));
  //}
"


m =stan_model(model_code=model_prior)
# m = stan_model(model_code = model_prior)

# small_simt <- small_sim
small_sim <- small_simt[1:20,]
dat = list(N = nrow(small_sim),
           n = small_sim$n,
           c = small_sim$c,
           x = small_sim$x,
           n_params = 6,
           m0 = theta,
           C0 = sigs,
           betas = theta)
rmle = sampling(m, dat, chains=1,iter=10000,warmup=5000)


par <- rstan::extract(rmle)#,pars=c('theta'))
thetas <- do.call('cbind',par[1:6])
traceplot(rmle)#,pars=c('theta'))

predps <- apply(thetas,MARGIN=2,FUN=mean)
w <- seq(0,53,length.out=1001)
y <- inv.logit(ASG(w,predps[1],predps[2],predps[3],predps[4],
                   predps[5]^2,predps[6]^2))

small_sim %>% 
  ggplot() +
  geom_point(aes(x=x,y=c/n)) +
  geom_line(data=data.frame(w,y),aes(x=w,y=y),colour='red')
# plot(y~w,type='l')

thetadf <- data.frame(thetas)
colnames(thetadf) <- c('beta1','beta2','eta','mu','sig1','sig2')
pairs(thetadf)

library(gplots)
hist2d(x=thetadf$mu,y=thetadf$sig2)

library(ggplot2)
thetadf %>%
  ggplot(aes(x=mu,y=sig2)) +
  geom_hex() +
  theme_bw()










