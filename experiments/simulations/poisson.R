library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

model = "poisson"
n = 500 # Number of clusters
k = 7 # Observations per cluster
sigma = 1.5 # standard deviation of random effects


# Models of the form
# g(mu_{ij}) = beta0 + beta1 * x_{ij} + b_i
# with b_i ~ N(0, sigma^2)
# g(x) = log(x)


### Poisson model 1 - zero inflated
beta0 = -2.5
beta1 = -2

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n) # x_{ij} = (j - 4)/10 from Tan (2021)
Z = vector(mode = "list", length = n) # Will just be vector of 1's
etahat = vector(mode = "list", length = n)

x_ij = (c(1:k) - 4)/10

set.seed(987)
for (i in 1:n){
  bi = rnorm(1, 0, sigma)
  eta = beta0 + beta1*x_ij + bi
  mu = exp(eta)
  y[[i]] = rpois(k, mu)
  X[[i]] = matrix(c(rep(1, k), x_ij), nrow = k)
  Z[[i]] = matrix(rep(1, k), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
xvec = rep(x_ij, n)
id = rep(1:n, each = k)

mean(yvec == 0) # 83.8%
mean(yvec <= 5) # 99.5%

pois1df <- data.frame(y = yvec, x = xvec, id = id)

pois1glm <- glm(y ~ x, data = pois1df, family = poisson)
pois1pred <- predict(pois1glm, type = "response")

pois1prior <- rvb::KNprior(model, pois1pred, Z)

# Run RVB

pois1_RVB1 <- rvb::Alg_RVB1(y, X, Z, pois1prior, etahat, model)
pois1_RVB2 <- rvb::Alg_RVB2(y, X, Z, pois1prior, etahat, model)

# INLA
pois1_inlaprior <- list(prec = list(param = c(pois1prior$nu/2, pois1prior$Sinv/2)))
pois1_INLA <- inla(y ~ x + f(id, model = "iid", hyper = pois1_inlaprior),
                   data = pois1df, family = "poisson")
summary(pois1_INLA)

# Stan
pois1_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = pois1prior$nu, S = (pois1prior$S %>% as.matrix),
                     binom = 0, n_binom = rep(1,n))

pois1_stanfit <- stan(file = "experiments/stan/glmm.stan",
                      data = pois1_standt, chains = 4, iter = 25000,
                      warmup = 12500, cores = 4)

pois1_stansum <- summary(pois1_stanfit)
pois1_stansum$summary[1:2]

traceplot(pois1_stanfit, pars = c("beta[1]", "beta[2]"))

rvb::summary_table(pois1_RVB1, pois1_RVB2, pois1_INLA, pois1_stanfit, n, p, r)

##############################################################################

### Poisson model 2
beta0 = 1.5
beta1 = 0.5

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n)
Z = vector(mode = "list", length = n)
etahat = vector(mode = "list", length = n)

set.seed(765)
for (i in 1:n){
  bi = rnorm(1, 0, sigma)
  eta = beta0 + beta1*x_ij + bi
  mu = exp(eta)
  y[[i]] = rpois(k, mu)
  X[[i]] = matrix(c(rep(1, k), x_ij), nrow = k)
  Z[[i]] = matrix(rep(1, k), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
xvec = rep(x_ij, n)
id = rep(1:n, each = k)

mean(yvec == 0) # 14.6%
mean(yvec == 1) # 12.3%
mean(yvec == 2) # 9.62%
mean(yvec == 3) # 8.23%
mean(yvec == 4) # 5.83%
mean(yvec >= 5) # 49.4%

pois2df <- data.frame(y = yvec, x = xvec, id = id)

pois2glm <- glm(y ~ x, data = pois2df, family = poisson)
pois2pred <- predict(pois2glm, type = "response")

pois2prior <- rvb::KNprior(model, pois2pred, Z)

# Run RVB

pois2_RVB1 <- rvb::Alg_RVB1(y, X, Z, pois2prior, etahat, model)
pois2_RVB2 <- rvb::Alg_RVB2(y, X, Z, pois2prior, etahat, model)

# INLA
pois2_inlaprior <- list(prec = list(param = c(pois2prior$nu/2, pois2prior$Sinv/2)))
pois2_INLA <- inla(y ~ x + f(id, model = "iid", hyper = pois2_inlaprior),
                   data = pois2df, family = "poisson")
summary(pois2_INLA)

# Stan
pois2_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = pois2prior$nu, S = (pois2prior$S %>% as.matrix),
                     binom = 0, n_binom = rep(1,n))

pois2_stanfit <- stan(file = "experiments/stan/glmm.stan",
                      data = pois2_standt, chains = 4, iter = 25000,
                      warmup = 12500, cores = 4)

pois2_stansum <- summary(pois2_stanfit)
pois2_stansum$summary[1:2]

rvb::summary_table(pois2_RVB1, pois2_RVB2, pois2_INLA, pois2_stanfit, n, p, r)

traceplot(pois2_stanfit, pars = c("beta[1]", "beta[2]"))
