library(rvb)
library(INLA)
library(rstan)
library(tidyverse)

#setwd("experiments/data analysis/")
################################################################################


### Data preparation
model <- "poisson"
epil <- MASS::epil

# See Section 8.2 of Tan (2021) for descriptions
sz <- epil$y
trt <- as.numeric(epil$trt) - 1
lb4 <- log(epil$base/4)
lb4trt <- lb4*trt
clage <- center(log(epil$age))
V4 <- epil$V4
Visit <- c(-0.3, -0.1, 0.1, 0.3)
id <- epil$subject
unique_id <- unique(id)

n = length(unique_id) # Number of subjects
vni = rep(4, n) # Number of observations per subject
y <- vector(mode = "list", length = n) # List of responses
etahat <- vector(mode = "list", length = n)

# Making lists of response values and initializing etahat
for (i in 1:n){
  rows = which(id == unique_id[i])
  y[[i]] = sz[rows]
  etahat[[i]] = digamma(y[[i]] + 0.5)
}
yvec = unlist(y)
y_stats <- c(min(yvec), max(yvec), mean(yvec), median(yvec)) # 0, 102, 8.25, 4

df = data.frame(sz = sz, lb4 = lb4, trt = trt, clage = clage, Visit = rep(Visit, n),
                V4 = V4, id = id)

SS = 50000 # Sample size

################################################################################

#### Random Intercept Model

X = vector(mode = "list", length = n)
Z = vector(mode = "list", length = n)

# Lists of covariate data
for (i in 1:n){
  rows = which(id == unique_id[i])
  Z[[i]] = matrix(1, nrow = 4) # random intercept
  X[[i]] = matrix(c(rep(1,4), lb4[rows], trt[rows], (lb4[rows]*trt[rows]),
                    clage[rows], V4[rows]), nrow = 4) # Design matrix
}

p = ncol(X[[1]]) # dimension of fixed effects
r = ncol(Z[[1]]) # dimension of random effects
G = p + 0.5*r*(r + 1) # Total number of global variables

set.seed(123)

glmfit <- glm(sz ~ lb4 + trt + I(lb4*trt) + clage + V4, data = df, family = poisson(link = "log"))
summary(glmfit)
pred = predict(glmfit, type = "response")

Wprior1 <- rvb::KNprior(model, pred, Z) # Getting prior

# Run RVB1 and RVB2
RVB1_int <- Alg_RVB1(y, X, Z, Wprior1, etahat, model)
RVB2_int <- Alg_RVB2(y, X, Z, Wprior1, etahat, model)

# Compare with INLA
prec.prior <- list(prec = list(param = c(Wprior1$nu/2, Wprior1$Sinv/2)))
INLA1 <- inla(sz ~ lb4 + trt + I(lb4*trt) + clage + V4 + f(id, model = "iid",
                                                           hyper = prec.prior),
              data = df, control.predictor = list(compute = T),
              family = "poisson")
summary(INLA1)

# Compare with Stan
stan_data_1 <- list(M = n*4, N = n, K = 4, P = p, R = r, y = yvec,
                  x = X, z = Z, sdbeta0 = 10,
                  nu = Wprior1$nu, S = (Wprior1$S %>% as.matrix()),
                  binom = 0, n_binom = rep(1,n))

stan1 <- stan(file = "../stan/glmm.stan",
              data = stan_data_1, chains = 4, iter = 25000, warmup = 12500,
              cores = 4)

stan1sum <- summary(stan1)

traceplot(stan1, pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]",
                          "beta[6]"))



bs <- rvb::postb_RVB1_b1(RVB1, y, X, Z, etahat, SS, L = n, model, m = 1)
mean1_rvb1 <- rowMeans(bs)
sd1_rvb1 <- matrixStats::rowSds(bs)

boxplot(mean1_rvb1, sd1_rvb1)

###############################################################################

### Random Intercept + Slope

X = vector(mode = "list", length = n)
Z = vector(mode = "list", length = n)

# Same as earlier
for (i in 1:n){
  rows = which(id == unique_id[i])
  Z[[i]] = matrix(c(rep(1,4), Visit), nrow = 4)
  X[[i]] = matrix(c(rep(1,4), lb4[rows], trt[rows], (lb4[rows]*trt[rows]),
                    clage[rows], Visit), nrow = 4)
}

p = ncol(X[[1]])
r = ncol(Z[[1]])
G = p + 0.5*r*(r + 1)

glmfit <- glm(sz ~ lb4 + trt + I(lb4*trt) + clage + Visit,
              data = df, family = poisson(link = "log"))
summary(glmfit)
pred = predict(glmfit, type = "response")

Wprior2 <- rvb::KNprior(model, pred, Z) # New prior

RVB1_slope <- rvb::Alg_RVB1(y, X, Z, Wprior2, etahat, model, m = 1)
RVB2_slope <- rvb::Alg_RVB2(y, X, Z, Wprior2, etahat, model, m = 1)


# Compare with INLA
df$id2 <- df$id + n
prec.prior2 <- list(theta1 = list(param = c(Wprior2$nu,
                                            Wprior2$Sinv[1,1],
                                            Wprior2$Sinv[2,2],
                                            Wprior2$Sinv[1,2])))

INLA2 <- inla(sz~lb4 + trt + I(lb4*trt) + clage + Visit +
                f(id,model = "iid2d", hyper = prec.prior2, n = 2*n) +
                f(id2, Visit, copy = "id"),
              data = df, family = "poisson", control.predictor = list(compute = T))
summary(INLA2)

# Compare with Stan
stan_data_2 <- list(M = n*4, N = n, K = 4, P = p, R = r, y = yvec,
                    x = X, z = Z, sdbeta0 = 10,
                    nu = Wprior2$nu, S = (Wprior2$S %>% as.matrix()),
                    binom = 0, n_binom = rep(1,n))

stan2 <- stan(file = "../stan/glmm.stan",
              data = stan_data_2, chains = 4, iter = 25000, warmup = 12500,
              cores = 4)

stan2sum <- summary(stan2)

traceplot(stan2, pars = c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]",
                          "beta[6]"))

