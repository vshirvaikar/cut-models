# StatML Module 1 Project - Bayesian Modelling and Computation
# Efthymios Costa, Henry Antonio Palasciano, and Vik Shirvaikar

# Based on code from Jacob, P. E., O’Leary, J., and Atchade, Y. F., 
# "Unbiased Markov chain Monte Carlo with couplings", 2017
# and data from Maucort-Boulch, D., Franceschi, S., Plummer, M.,
# "International Correlation between Human Papillomavirus Prevalence 
# and Cervical Cancer Incidence", 2008

set.seed(2021)
library(MASS)

# Module 1 uses HPV data to form a Beta distribution for phi
J <- 13
Z_hpv <- c(32, 0, 25, 46, 51, 127, 51, 38, 4, 4, 41, 173, 164)
Z_pop <- c(911, 634, 570, 3260, 802, 1490, 830, 711, 518, 755, 678, 1601, 812)
phi_alpha <- 1 + Z_hpv
phi_beta <- 1 + Z_pop - Z_hpv

module1_sample <- function(nsamples){
  phis <- matrix(nrow = nsamples, ncol = J)
  for (j in 1:J){
    phis[,j] <- rbeta(nsamples, shape1 = phi_alpha[j], shape2 = phi_beta[j])
  }
  return(phis)
}

# Module 2 uses cervical cancer data and phi from module 1 to estimate theta
Y_cerv <- c(17.5, 18.6, 38.2, 11.2, 59.3, 85.7, 43.6, 61.5,
            50.5, 17.6, 128.4, 231.5, 175.4)
Y_pop <- c(100000, 100000, 100000, 100000, 100000, 100000, 100000,
           100000, 100000, 100000, 100000, 100000, 100000)
Y_pop_norm <- log(10**(-3) * Y_pop)

hyper2 <- list(theta_mean_prior = 0, theta_sd_prior = sqrt(1000))
dprior2 <- function(theta, hyper2){
  return(sum(dnorm(theta, mean = hyper2$theta_mean_prior, 
                   sd = hyper2$theta_sd_prior, log = TRUE)))
}

module2_loglike <- function(phi, theta, ncases, Y_pop_norm){
  score <- mu <- logmu <- 0
  for (j in 1:J){
    logmu = theta[1] + phi[j] * theta[2] + Y_pop_norm[j]
    mu = exp(logmu)
    score = score + Y_cerv[j] * logmu - mu
  }
  return(score)
}

get_kernels <- function(phi, sigma_proposal, init_mean, init_sigma){
  target <- function(x){
    return(module2_loglike(phi, x, Y_cerv, Y_pop_norm) + dprior2(x, hyper2))
  }

  single_kernel <- function(state){
    chain_state <- state$chain_state
    current_pdf <- state$current_pdf
    proposal_value <- chain_state + mvrnorm(1, c(0, 0), chol(sigma_proposal))
    proposal_pdf <- target(proposal_value)
    if (log(runif(1)) < (proposal_pdf - current_pdf)){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf))
    }
  }
  
  rinit <- function(){
    chain_state <- mvrnorm(1, init_mean, init_sigma)
    current_pdf <- target(chain_state)
    return(list(chain_state = chain_state, current_pdf = current_pdf))
  }
  return(list(target = target, single_kernel = single_kernel, rinit = rinit))
}

# Sample phis and compute posterior mean under first model
nsamples <- 1000
phis <- module1_sample(nsamples)
phi_hat <- colMeans(phis)

## Run standard MCMC for inference on theta given phi_hat
sigma_proposal <- diag(0.1, 2, 2)
init_mean <- rep(0, 2)
init_sigma <- diag(1, 2, 2)
pb <- get_kernels(phi_hat, sigma_proposal, init_mean, init_sigma)
niterations <- 5e3
chain <- matrix(0, nrow = niterations, ncol = 2)

state <- pb$rinit()
for (iteration in 1:niterations){
  state <- pb$single_kernel(state)
  chain[iteration,] <- state$chain_state
}

matplot(chain, type = "l")

hist(chain[,2], probability = TRUE, breaks=seq(min(chain[,2]), 20, by=.5), 
     ylim=c(0,1), xlim=c(min(chain[,2]),15),
     main=bquote("Cut posterior approximation for" ~ theta[2]),
     xlab=expression(theta[2]), ylab='Density')
d <- density(chain[,2], bw=0.25)
lines(d,col='red',lwd=2, xlim=c(min(chain[,2]),20))

hist(chain[,1], probability = TRUE, breaks=seq(-2, 0, by=.1), 
     ylim=c(0,7), xlim=c(-2, 0),
     main=bquote("Cut posterior approximation for" ~ theta[1]),
     xlab=expression(theta[1]), ylab='Density')
d <- density(chain[,1], bw=0.05)
lines(d,col='red',lwd=2, xlim=c(-2, 0))
