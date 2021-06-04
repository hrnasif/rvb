library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

setwd("experiments/simulations/")


model = "binomial"
n = 200 # Number of clusters. Decreased from sims in paper
k = 7 # Observations per cluster
sigma1 = 1.5 # standard deviation of random intercept
sigma2 = 1 # standard deviation for random slope


# Models of the form
# g(mu_{ij}) = beta0 + beta1 * x_{ij} + b_i
# with b_i ~ N(0, sigma^2)
# g(x) = logit(x)


### Bernoulli model 1 - concentrated at boundaries
beta0 = -2.5
beta1 = 4.5

# new time effect
beta2 = -2
x_ijt = 0.2*(1:k)

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n) # x_{ij} ~ Bernoulli(0.5)
Z = vector(mode = "list", length = n) # intercept and slope
etahat = vector(mode = "list", length = n)

set.seed(409)

for (i in 1:n){
  bi1 = rnorm(1, 0, sigma1)
  bi2 = rnorm(1, 0, sigma2)
  x_ij = rbinom(k, size = 1, prob = 0.5)
  eta = beta0 + beta1*x_ij + (beta2 + bi2)*x_ijt + bi1
  mu = plogis(eta) # expit
  y[[i]] = rbinom(k, size = 1, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij, x_ijt), nrow = k)
  Z[[i]] = matrix(c(rep(1, k), x_ijt), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
xvec1 = sapply(1:n, function(x){X[[x]][,2]}) %>% as.vector()
xvec2 = sapply(1:n, function(x){X[[x]][,3]}) %>% as.vector()
id = rep(1:n, each = k)

bern1exdf <- data.frame(y = yvec, x1 = xvec1, x2 = xvec2, id = id)

bern1exglm <- glm(y ~ x1 + x2, data = bern1exdf, family = binomial)
bern1expred <- predict(bern1exglm, type = "response")

bern1exprior <- rvb::KNprior(model, bern1expred, Z)

# Run RVB

bern1ex_RVB1 <- rvb::Alg_RVB1(y, X, Z, bern1exprior, etahat, model)
bern1ex_RVB2 <- rvb::Alg_RVB2(y, X, Z, bern1exprior, etahat, model)

# INLA
bern1exdf$id2 <- bern1exdf$id + n
bern1ex_inlaprior <- list(theta1 = list(param = c(bern1exprior$nu,
                                            bern1exprior$Sinv[1,1],
                                            bern1exprior$Sinv[2,2],
                                            bern1exprior$Sinv[1,2])))

bern1ex_INLA <- inla(y ~ x1 + x2 + f(id, model = "iid2d", hyper = bern1ex_inlaprior, n = 2*n) +
                       f(id2, x2, copy = "id"), data = bern1exdf, family = "binomial")
summary(bern1ex_INLA)

# Stan
bern1ex_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = bern1exprior$nu, S = (bern1exprior$S %>% as.matrix),
                     binom = 1, n_binom = rep(1,n))

bern1ex_stanfit <- stan(file = "../stan/glmm.stan",
                      data = bern1ex_standt, chains = 4, iter = 10000,
                      warmup = 8000, cores = 4)

bern1ex_stansum <- summary(bern1ex_stanfit)
bern1ex_stansum$summary[1:3]


##############################################################################

### bernoulli model 2
beta0 = 0
beta1 = 2

# new time effect
beta2 = -0.5
sigma2 = 1
x_ijt = 0.2*(1:k)

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n)
Z = vector(mode = "list", length = n)
etahat = vector(mode = "list", length = n)

#x_ij = (c(1:k) - 4)/10

set.seed(111)
for (i in 1:n){
  bi1 = rnorm(1, 0, sigma1)
  bi2 = rnorm(1, 0, sigma2)
  x_ij = rbinom(k, 1, 0.5)
  eta = beta0 + beta1*x_ij + (bi2 + beta2)*x_ijt + bi1
  mu = plogis(eta)
  y[[i]] = rbinom(k, size = 1, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij, x_ijt), nrow = k)
  Z[[i]] = matrix(c(rep(1, k), x_ijt), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
xvec1 = sapply(1:n, function(x){X[[x]][,2]}) %>% as.vector()
xvec2 = rep(x_ijt, n)
id = rep(1:n, each = k)

bern2exdf <- data.frame(y = yvec, x1 = xvec1, x2 = xvec2, id = id)

bern2exglm <- glm(y ~ x1 + x2, data = bern2exdf, family = binomial)
bern2expred <- predict(bern2exglm, type = "response")

bern2exprior <- rvb::KNprior(model, bern2expred, Z)

# Run RVB

bern2ex_RVB1 <- rvb::Alg_RVB1(y, X, Z, bern2exprior, etahat, model)
bern2ex_RVB2 <- rvb::Alg_RVB2(y, X, Z, bern2exprior, etahat, model)

# INLA
bern2exdf$id2 <- bern2exdf$id + n
bern2ex_inlaprior <- list(theta1 = list(param = c(bern2exprior$nu,
                                                  bern2exprior$Sinv[1,1],
                                                  bern2exprior$Sinv[2,2],
                                                  bern2exprior$Sinv[1,2])))

bern2ex_INLA <- inla(y ~ x1 + x2 + f(id, model = "iid2d", hyper = bern2ex_inlaprior, n = 2*n) +
                       f(id2, x2, copy = "id"), data = bern2exdf, family = "binomial")
summary(bern2ex_INLA)

# Stan
bern2ex_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                     x = X, z = Z, sdbeta0 = 10,
                     nu = bern2exprior$nu, S = (bern2exprior$S %>% as.matrix),
                     binom = 1, n_binom = rep(1,n))

bern2ex_stanfit <- stan(file = "../stan/glmm.stan",
                      data = bern2ex_standt, chains = 4, iter = 10000,
                      warmup = 8000, cores = 4)

bern2ex_stansum <- summary(bern2ex_stanfit)
bern2ex_stansum$summary[1:3]
