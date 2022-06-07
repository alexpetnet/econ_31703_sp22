### sasha petrov - - pset 2 for bonhomme's topics in ecm (spring 2022)

library(MASS)

lasso.objective <- function(b, data, lambda) {
  .5 * (1 / nrow(data)) * sum((data[, 1] - data[, -1] %*% b) ^ 2) +
    lambda * sum(abs(b))
}

lasso.update <- function(b, data, lambda) {
  for (k in 1:length(b)) {
    ytilde <- data[, 1] - data[, -c(1, k + 1)] %*% b[-k]
    cov <- ytilde[, 1] %*% data[, k + 1]
    var <- data[, k + 1] %*% data[, k + 1]
    b[k] <- ifelse(abs(cov / dim(data)[1]) < lambda, 0, 
                   cov / var - sign(cov) * lambda )
  }
  b
}

# standardise function

sdata <- function(data) {
  data[, 1] <- data[, 1] - mean(data[, 1])
  data[, -1] <- apply(data[, -1], 2, function(x){(x - mean(x)) / sqrt( mean((x - mean(x))^2) ) })
  data
}


lasso <- function(data, lambda, b_initial = rep(0, p),
                  eps = 1e-06, max = 1000) {
  # standardise
  data <- sdata(data)
  
  b <- matrix(b_initial, nrow = 2, ncol = length(b_initial), byrow = TRUE)
  o <- c()
  
  crit <- 0
  i <- 1
  while (i <= max) {
    i <- i + 1
    b[2, ] <- b[1, ]
    b[1, ] <- lasso.update(b[2, ], data, lambda)
    o <- c(o, lasso.objective(b[1, ], data, lambda))
    if (max(abs(b[1, ] - b[2, ])) < eps) {
      crit = 1
      break
    }
  }
  
  ans <- list(b[1, ], o,
              ifelse(crit, 'eps', 'max'))
  
  return(ans)
  
}

# sampling from the pset 1

sample_sn <- function(p, n, N) {
  array(rnorm(n * (p + 1) * N), dim = c(n, p + 1, N))
}

sample_d <- function(srn, sigma_chol) {
  s <- apply(srn, 3, function(x) {x %*% sigma_chol}, simplify = FALSE)
  lapply(s, function(x) {cbind(x[, 1] - x[, 2] + x[, dim(srn)[2]], x[, 1:(dim(srn)[2] - 1)])})
}

p <- 90
n <- 100

N <- 10000 # number of samples

rns <- sample_sn(p, n, N)
mc <- sample_d(rns, diag(p + 1))

# calculating ols coefs

ols <- function(sample, nreg) {
  ginv(sample[, 2:(nreg+1)]) %*% sample[, 1]
}

bols <- sapply(mc, function(x){ols(x, p)[1:2]})
blasso <- lapply(mc, function(x){lasso(x, 20)})

ridge <- function(data, lambda) {
  x <- data[, -1]
  y <- data[, 1]
  solve(t(x) %*% x + lambda * diag(dim(x)[2])) %*% (t(x) %*% y)
}

bridge <- sapply(mc, function(x){ridge(x, 20)})

## are there any non-zero coefs in lasso? 
## i don't think so for \lambda = 20 lol

zero <- sapply(blasso, function(x) {TRUE %in% (x[[1]] != 0)})
print(paste('how many times do we get non-zero coefs?', sum(zero)))

## storing and comparing coef estimates

comp <- matrix(0, 2, 3)

comp[, 1] <- apply(bols, 1, mean)
comp[, 2] <- apply(do.call(cbind, blasso[1, ]), 1, mean)[1:2]
comp[, 3] <- apply(bridge[1:2, ], 1, mean)

### 1f - - plotting coefs against penalisation intensity

s <- sample(1:N, 1)

# looking for the lowest lambda that yields all 0s

covt <- function(data, k, b = rep(0, p)) {
  ytilde <- data[, 1] - data[, -c(1, k + 1)] %*% b[-k]
  ytilde[, 1] %*% data[, k + 1] / dim(data)[1]
}

ind <- 1:p
covtilde <- covt(sdata(mc[[s]]), ind)
lmax <- max(abs(covtilde))

grid <- seq(.01, 1, .01)

coefs_1f <- matrix(grid, p, length(grid), byrow = TRUE)

coefs_1f <- apply(coefs_1f, 2, function(x){lasso(mc[[s]], x[1] * lmax)[[1]]})

plot(0, 0, xlim = c(0, 1), ylim = c(-1, 1), type = 'n',
     xlab = expression(lambda), ylab = expression(hat(beta)))
cl <- rainbow(5)
for (i in 1:5) {
  lines(grid, coefs_1f[i, ], col = cl[i])
}
legend("topright", legend = 1:5, col=cl, pch=1)



### exercise 2


# keep in mind that the size of the dataset is T + 1 !!!
# but the number of data points is still T !!! eheh

mspe <- function(data, lambda, tau = c(100, 10), h = 1) {
  # i'm reversing the order of observations bc i rely later on
  # the first obs (row) being the earliest period
  Y <- data[nrow(data):1, 1]
  X <- data[nrow(data):1, -1]
  
  B <- dim(data)[1] - 1 - tau[1] - h + 1 # the length of Yhat
  
  Yhat = c() # in the end will be of length B
  
  i <- 0
  while (length(Yhat) < B) {
    start1 <- 1 + i * tau[2]
    start2 <- tau[1] + 1 + i * tau[2]
    
    xt <- X[start1:(start1 + tau[1] - 1), ] # should never face the problem of the right bound exceeding T + 1, right?
    ytph <- Y[(start1 + h):(start1 + h + tau[1] - 1)]
    d <- cbind(ytph, xt)
    bl <- lasso(d, lambda, b_initial = rep(0, dim(xt)[2]))[[1]]
    
    rb <- min(dim(data)[1] - h, start2 + tau[2] - 1)
    yh <- X[start2:rb, ] %*% bl
    Yhat <- c(Yhat, yh)
    
    i <- i + 1
  }
  
  mean( (Y[(length(Y) - B + 1):length(Y)] - Yhat)^2 )
}

## sampling

T <- 200
no <- T + 1
nr <- 50
N2 <- 20 # the number of samples

rho <- .9
sigma <- diag(nr) * (1 - rho ^ 2)
beta <- c(5, rep(0, nr - 1))

eta <- rnorm(no)
epsilon <- matrix(rnorm(no * nr), no, nr)
epsilon[-no, ] <- epsilon[-no, ] %*% chol(sigma)

X <- epsilon
for (i in 1:(no - 1)) {
  X[no - i, ] <- rho * X[no - i + 1, ] + X[no - i, ]
}

Y <- X %*% beta + eta

dataset <- cbind(Y, X)

sigma_b <- diag(nr) * (1 - rho ^ 2)

s2 <- function(T, nr, sigma, beta) {
  sigma_0 <- sigma * (1 - rho ^ 2) ^ (-1)
  
  no <- T + 1
  
  eta <- rnorm(no)
  epsilon <- matrix(rnorm(no * nr), no, nr)
  # THE LAST ROW IS THE FIRST PERIOD DON'T ASK WHY
  epsilon[-no, ] <- epsilon[-no, ] %*% chol(sigma)
  epsilon[no, ] <- epsilon[no, ] %*% chol(sigma_0)
  
  X <- epsilon
  for (i in 1:(no - 1)) {
    X[no - i, ] <- rho * X[no - i + 1, ] + X[no - i, ]
  }
  
  Y <- X %*% beta + eta
  
  cbind(Y, X)
}

mc2 <- vector('list', N2)
mc2 <- lapply(mc2, function(x) {s2(T, nr, rho, beta)})


lambdas <- seq(.01, 1.5, length.out = 20)

cross.validation <- function(data, lambda_seq, tau = c(100, 10), h = 1) {
  sapply(lambda_seq, function(x) {mspe(data, x, tau, h)})
}

mspe_seq <- sapply(mc2, function(x) {cross.validation(x, lambdas)})

av_mspe <- apply(mspe_seq, 1, mean)

plot(lambdas, av_mspe, ylim = c(7.8, 12),
     xlab = expression(lambda), ylab = expression(bar(MSPE)))

lambda_cv <- lambdas[which.min(av_mspe)]

lasso.forecast <- function(data, lambda, h = 1) {
  yth <- data[(nrow(data) - h):1, 1]
  xt <- data[nrow(data):(1 + h), -1]
  d <- cbind(yth, xt)
  lasso(d, lambda, b_initial = rep(0, ncol(xt)))[[1]]
}

end_lasso <- sapply(mc2, function(x) {lasso.forecast(x, lambda_cv)})

table_2b <- apply(end_lasso, 1, mean)

## 2c - introducing some dependence

tilde_sigma <- matrix(.8 * (1 - rho ^ 2), 10, 10)
diag(tilde_sigma) <- (1 - rho ^ 2)
sigma_c <- kronecker(diag(5), tilde_sigma)

mc2c <- vector('list', N2)
mc2c <- lapply(mc2c, function(x) {s2(T, nr, sigma_c, beta)})

mspe_seq_c <- sapply(mc2c, function(x) {cross.validation(x, lambdas)})

av_mspe_c <- apply(mspe_seq_c, 1, mean)

plot(lambdas, av_mspe_c, ylim = c(7.8, 11),
     xlab = expression(lambda), ylab = expression(bar(MSPE)))

lambda_cv_c <- lambdas[which.min(av_mspe_c)]


end_lasso_c <- sapply(mc2c, function(x) {lasso.forecast(x, lambda_cv_c)})

table_2c <- apply(end_lasso_c, 1, mean)


## exporting all results

save(bols, blasso, bridge, comp, grid, coefs_1f, lambdas, av_mspe,
     av_mspe_c, table_2b, table_2c,
     file = 'output.RData')
