#simple simulation of asg to assess our MCMC fitting
library(mixtools)


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


model_prior = "
functions {
   real asg(row_vector param, real eta, real x) {
    real beta1;
    real beta2;
    //real eta;
    real mu;
    real sig21;
    real sig22;
    real ASG;
    
    beta1 = param[1];
    beta2 = param[2];
    //eta = param[3];
    mu = param[4];
    //print(mu);
    sig21 = param[5];
    sig22 = param[6];
    
    ASG = (beta1 + (eta-beta1)*exp(-(x-mu)^2/(2*sig21^2)))*(x <mu) +
       (beta2 + (eta-beta2)*exp(-(x-mu)^2/(2*sig22^2)))*(x >= mu);
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
matrix[n_params,n_params] C0;
real eta;
}
parameters {

row_vector[n_params] theta;

}

model {

theta ~ multi_normal(m0,C0);

for (i in 1:N) c[i] ~ binomial(n[i],inv_logit(asg(theta, eta, x[i]))); 

}

generated quantities {
    vector[N] postpred_pr;
    for (i in 1:N)
      postpred_pr[i] = binomial_rng(n[i],inv_logit(asg(theta, eta,
                                                    x[i])));
  }
"
m =stan_model(model_code=model_prior)
# m = stan_model(model_code = model_prior)


dat = list(N = nrow(small_sim),
           n = small_sim$n,
           c = small_sim$c,
           x = small_sim$x,
           n_params = 6,
           m0 = theta,
           C0 = diag(sigs^2),
           eta = theta[3])
rmle = sampling(m, dat, chains=1,iter=10000,warmup=4000)


par <- rstan::extract(rmle,pars=c('theta'))
thetas <- par$theta
traceplot(rmle,pars=c('theta'))

predps <- apply(thetas,MARGIN=2,FUN=mean)
w <- seq(0,53,length.out=1001)
y <- inv.logit(ASG(w,predps[1],predps[2],predps[3],predps[4],
                   predps[5]^2,predps[6]^2))

small_sim %>% 
  ggplot() +
  geom_point(aes(x=x,y=c/n)) +
  geom_line(data=data.frame(w,y),aes(x=w,y=y))
plot(y~w,type='l')






