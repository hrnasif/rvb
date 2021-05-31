library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

setwd("experiments/simulations/")


model = "binomial"
n = 500 # Number of clusters
k = 7 # Observations per cluster
sigma = 1.5 # standard deviation of random effects


# Models of the form
# g(mu_{ij}) = beta0 + beta1 * x_{ij} + b_i
# with b_i ~ N(0, sigma^2)
# g(x) = logit(x)


### Bernoulli model 1 - concentrated at boundaries
beta0 = -2.5
beta1 = 4.5

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n) # x_{ij} ~ Bernoulli(0.5)
Z = vector(mode = "list", length = n) # Will just be vector of 1's
etahat = vector(mode = "list", length = n)

set.seed(546)

for (i in 1:n){
  bi = rnorm(1, 0, sigma)
  x_ij = rbinom(k, size = 1, prob = 0.5)
  eta = beta0 + beta1*x_ij + bi
  mu = plogis(eta) # expit
  y[[i]] = rbinom(k, size = 1, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij), nrow = k)
  Z[[i]] = matrix(rep(1, k), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
xvec = sapply(1:n, function(x){X[[x]][,2]}) %>% as.vector()
id = rep(1:n, each = k)

bern1df <- data.frame(y = yvec, x = xvec, id = id)

bern1glm <- glm(y ~ x, data = bern1df, family = binomial)
bern1pred <- predict(bern1glm, type = "response")

bern1prior <- rvb::KNprior(model, bern1pred, Z)

# Run RVB

bern1_RVB1 <- rvb::Alg_RVB1(y, X, Z, bern1prior, etahat, model)
bern1_RVB2 <- rvb::Alg_RVB2(y, X, Z, bern1prior, etahat, model)

# INLA
bern1_inlaprior <- list(prec = list(param = c(bern1prior$nu/2, bern1prior$Sinv/2)))
bern1_INLA <- inla(y ~ x + f(id, model = "iid", hyper = bern1_inlaprior),
                   data = bern1df, family = "binomial")
summary(bern1_INLA)

# Stan
bern1_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = bern1prior$nu, S = (bern1prior$S %>% as.matrix),
                     binom = 1, n_binom = rep(1,n))

bern1_stanfit <- stan(file = "../stan/glmm.stan",
                      data = bern1_standt, chains = 4, iter = 10000,
                      warmup = 5000, cores = 4)

bern1_stansum <- summary(bern1_stanfit)
bern1_stansum$summary[1:2]


##############################################################################

### bernson model 2
beta0 = 0
beta1 = 1

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n)
Z = vector(mode = "list", length = n)
etahat = vector(mode = "list", length = n)

x_ij = (c(1:k) - 4)/10

set.seed(623)
for (i in 1:n){
  bi = rnorm(1, 0, sigma)
  eta = beta0 + beta1*x_ij + bi
  mu = plogis(eta)
  y[[i]] = rbinom(k, size = 1, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij), nrow = k)
  Z[[i]] = matrix(rep(1, k), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
xvec = rep(x_ij, n)
id = rep(1:n, each = k)

bern2df <- data.frame(y = yvec, x = xvec, id = id)

bern2glm <- glm(y ~ x, data = bern2df, family = binomial)
bern2pred <- predict(bern2glm, type = "response")

bern2prior <- rvb::KNprior(model, bern2pred, Z)

# Run RVB

bern2_RVB1 <- rvb::Alg_RVB1(y, X, Z, bern2prior, etahat, model)
bern2_RVB2 <- rvb::Alg_RVB2(y, X, Z, bern2prior, etahat, model)

# INLA
bern2_inlaprior <- list(prec = list(param = c(bern2prior$nu/2, bern2prior$Sinv/2)))
bern2_INLA <- inla(y ~ x + f(id, model = "iid", hyper = bern2_inlaprior),
                   data = bern2df, family = "binomial")
summary(bern2_INLA)

# Stan
bern2_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = bern2prior$nu, S = (bern2prior$S %>% as.matrix),
                     binom = 1, n_binom = rep(1,n))

bern2_stanfit <- stan(file = "../stan/glmm.stan",
                      data = bern2_standt, chains = 4, iter = 10000,
                      warmup = 5000, cores = 4)

bern2_stansum <- summary(bern2_stanfit)
bern2_stansum$summary[1:2]
