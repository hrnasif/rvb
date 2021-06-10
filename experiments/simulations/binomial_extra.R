library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

model = "binomial"
n = 200 # Number of clusters. Decreased from sims in paper
k = 7 # Observations per cluster
m = 20 # Binomial sample size
sigma1 = 1.5 # standard deviation of random intercept
sigma2 = 1 # standard deviation for random slope


# Models of the form
# g(mu_{ij}) = beta0 + beta1 * x_{ij} + b_i
# with b_i ~ N(0, sigma^2)
# g(x) = logit(x)


### binomial model 1 - concentrated at boundaries
beta0 = -2.5
beta1 = 4.5

# new time effect
beta2 = -2
x_ijt = 0.2*(1:k)

# Create simulation data
y = vector(mode = "list", length = n)
X = vector(mode = "list", length = n) # x_{ij} ~ binomial(0.5)
Z = vector(mode = "list", length = n) # intercept and slope
etahat = vector(mode = "list", length = n)

set.seed(4681)

for (i in 1:n){
  bi1 = rnorm(1, 0, sigma1)
  bi2 = rnorm(1, 0, sigma2)
  x_ij = rbinom(k, size = 1, prob = 0.5)
  eta = beta0 + beta1*x_ij + (beta2 + bi2)*x_ijt + bi1
  mu = plogis(eta) # expit
  y[[i]] = rbinom(k, size = m, prob = mu)
  X[[i]] = matrix(c(rep(1, k), x_ij, x_ijt), nrow = k)
  Z[[i]] = matrix(c(rep(1, k), x_ijt), nrow = k)
  etahat[[i]] = digamma(y[[i]] + 0.5) + digamma(1 - y[[i]] + 0.5)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])


# Obtain prior
yvec = unlist(y)
x1vec = sapply(1:n, function(x){X[[x]][,2]}) %>% as.vector()
x2vec = sapply(1:n, function(x){X[[x]][,3]}) %>% as.vector()
id = rep(1:n, each = k)

binom1exdf <- data.frame(y = yvec, x1 = xvec1, x2 = xvec2, id = id)

# transforming to long data for glm
m_total <- m*n*k
m_array <- rep(m, n*k)
startindex <- c(1, cumsum(m_array)[-(n*k)] + 1)
endindex <- cumsum(m_array)

y_long <- numeric(m_total)
id_long <- numeric(m_total)
x1_long <- numeric(m_total)
x2_long <- numeric(m_total)

for (i in 1:n){
  for (j in 1:k){
    index <- (i-1)*k + j
    rows <- startindex[index]:endindex[index]
    y_long[rows] <- c(rep(1, yvec[index]), rep(0, m - yvec[index]))
    id_long[rows] <- rep(i, m)
    x1_long[rows] <- rep(x1vec[index], m)
    x2_long[rows] <- rep(x2vec[index], m)
  }
}

binom1exdf_long <- data.frame(y_long = y_long, x1_long = x1_long, x2_long = x2_long,
                              id_long = id_long)

binom1exglm <- glm(y_long ~ x1_long + x2_long, data = binom1exdf_long,
                   family = binomial)
binom1expred_long <- predict(binom1exglm, type = "response")
binom1expred <- binom1expred_long[startindex]

binom1exprior <- rvb::KNprior(model, binom1expred, Z)

# Run RVB

binom1ex_RVB1 <- rvb::Alg_RVB1(y, X, Z, binom1exprior, etahat, model, m)
binom1ex_RVB2 <- rvb::Alg_RVB2(y, X, Z, binom1exprior, etahat, model, m)

# INLA
binom1exdf_long$id_long2 <- binom1exdf_long$id_long + n
binom1ex_inlaprior <- list(theta1 = list(param = c(binom1exprior$nu,
                                                  binom1exprior$Sinv[1,1],
                                                  binom1exprior$Sinv[2,2],
                                                  binom1exprior$Sinv[1,2])))

binom1ex_INLA <- inla(y_long ~ x1_long + x2_long +
                        f(id_long, model = "iid2d", hyper = binom1ex_inlaprior, n = 2*n) +
                        f(id_long2, x2_long, copy = "id_long"),
                      data = binom1exdf_long, family = "binomial")
summary(binom1ex_INLA)

# Stan
binom1ex_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                       x = X, z = Z, sdbeta0 = 10,
                       nu = binom1exprior$nu, S = (binom1exprior$S %>% as.matrix),
                       binom = 1, n_binom = rep(m,n))

binom1ex_stanfit <- stan(file = "experiments/stan/glmm.stan",
                        data = binom1ex_standt, chains = 4, iter = 10000,
                        warmup = 8000, cores = 4)

binom1ex_stansum <- summary(binom1ex_stanfit)
binom1ex_stansum$summary[1:3, c(1,3)]


##############################################################################

### binomial model 2
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

binom2exdf <- data.frame(y = yvec, x1 = xvec1, x2 = xvec2, id = id)

binom2exglm <- glm(y ~ x1 + x2, data = binom2exdf, family = binomial)
binom2expred <- predict(binom2exglm, type = "response")

binom2exprior <- rvb::KNprior(model, binom2expred, Z)

# Run RVB

binom2ex_RVB1 <- rvb::Alg_RVB1(y, X, Z, binom2exprior, etahat, model)
binom2ex_RVB2 <- rvb::Alg_RVB2(y, X, Z, binom2exprior, etahat, model)

# INLA
binom2exdf$id2 <- binom2exdf$id + n
binom2ex_inlaprior <- list(theta1 = list(param = c(binom2exprior$nu,
                                                  binom2exprior$Sinv[1,1],
                                                  binom2exprior$Sinv[2,2],
                                                  binom2exprior$Sinv[1,2])))

binom2ex_INLA <- inla(y ~ x1 + x2 + f(id, model = "iid2d", hyper = binom2ex_inlaprior, n = 2*n) +
                       f(id2, x2, copy = "id"), data = binom2exdf, family = "binomial")
summary(binom2ex_INLA)

# Stan
binom2ex_standt <- list(M = n*k, N = n, K = k, P = p, R = r, y = yvec,
                       x = X, z = Z, sdbeta0 = 10,
                       nu = binom2exprior$nu, S = (binom2exprior$S %>% as.matrix),
                       binom = 1, n_binom = rep(m,n))

binom2ex_stanfit <- stan(file = "experiments/stan/glmm.stan",
                        data = binom2ex_standt, chains = 4, iter = 10000,
                        warmup = 8000, cores = 4)

binom2ex_stansum <- summary(binom2ex_stanfit)
binom2ex_stansum$summary[1:3, c(1,3)]
