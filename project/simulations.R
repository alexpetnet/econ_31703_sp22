library(matlib)
library(MASS)

n <- 10

X <- runif(n)
U <- rnorm(n)
Y <- X^2 + U

data <- matrix(c(Y, X), n, 2)

plot(X, Y)
plot(X, U)

X2 <- runif(1000) ^ 2
Y2 <- X2^2 + U
plot(X2, Y2)

X3 <- runif(1000) ^ 3
Y3 <- X3^2 + U
plot(X3, Y3)

model.mc <- function()
  
kl <- function(p, q) {
  d <- log(p / q)
  mean(p * d)
}

yitzhaki <- function(t, draws) {
  v <- var(draws)
  e <- mean(draws[draws >= t]) - mean(draws[draws <= t])
  p <- mean(draws >= t) * mean(draws <= t)
  (1 / v) * e * p
} # is it okay to do inclusion for both terms ?

beta <- apply( matrix(X), 1, function(x){yitzhaki(x, X)} )
X.d <- dbeta(X, .5, 1)

kl.beta <- function(X, n) {
  beta <- apply( matrix(X), 1, function(x){yitzhaki(x, X)} )
  X.d <- dbeta(X, 1 / n, 1)
  kl(X.d, beta)
}

powers <- c(1:10)
draws.powers <- apply( matrix(powers), 1, function(x) {runif(10000) ^ x} )

kl.powers <- apply( draws.powers, 2, function(x) {kl.beta(x, )} )

kl.powers <- rep(0, 10)
for (i in 1:10) {
  kl.powers[i] <- kl.beta(draws.powers[ , i], i)
}

pdf('kl.pdf')
plot(powers, kl.powers, ylim = c(0, 300000))
dev.off()

# a function that plots yitzhaki weights and density
plot.yd <- function(X, n) {
  beta <- apply( matrix(X), 1, function(x){yitzhaki(x, X)} )
  X.d <- dbeta(X, 1 / n, 1)
  
  m <- matrix(c(X, beta, X.d), length(X), 3)
  m <- m[order(m[ , 1]), ]
  
  plot(m[, 1], m[, 3], ylim = c(0, 2))
  lines(m[, 1], m[, 2])
}

# a function that does ridge rkhs -- vapnik kernel
ridge <- function(data, lambda) {
  Y <- data[, 1]
  X <- data[, 2]
  N <- nrow(data)
  
  id <- matrix(1, N, N)
  xy <- X %*% t(X)
  x2y2 <- X^2 %*% t(X^2)
  K <- id + 2 * xy + x2y2
  
  alpha <- solve(K + lambda * id, Y)
  fhat <- K %*% alpha
  
  list( matrix(c(alpha, fhat), N, 2), K )
}

ridge.der <- function(data, alpha, lambda) {
  Y <- data[, 1]
  X <- data[, 2]
  N <- nrow(data)
  
  y <- t(matrix(X, N, N))
  xy2 <- X %*% t(X^2)
  Kder <- 2 * y + 2 * xy2
  Kder %*% alpha
}

ex <- ridge(data, 0)
ex <- cbind(ex, Y)

# a function that plots average derivative estimate as a function of lambda
plot.ad.lambda <- function(data, lambdas) {
  Y <- data[, 1]
  X <- data[, 2]
  N <- nrow(data)
  
  #alpha <- ridge(data, lambdas)
  #alphas <- apply(matrix(lambdas), 1, function(x){ridge(data, x)})
  
  #av.der <- mean(ridge.der(data, alpha, lambdas))
  avd <- apply(matrix(lambdas), 1, 
               function(x) 
                 {1 - mean(ridge.der(data, ridge(data, x)[, 1], x))})
  
  avd
} 

ls <- seq(1, 100, length.out = 10)

pdf('ad_lambdas.pdf')
plot(ls, plot.ad.lambda(data, ls))
dev.off()

plot.level.lambda <- function(data, lambdas) {
  Y <- data[, 1]
  X <- data[, 2]
  N <- nrow(data)
  
  #alpha <- ridge(data, lambdas)
  #alphas <- apply(matrix(lambdas), 1, function(x){ridge(data, x)})
  
  #av.der <- mean(ridge.der(data, alpha, lambdas))
  mse <- apply(matrix(lambdas), 1, 
               function(x) 
               {mean( (ridge(data, x)[, 2] - X ^ 2) ^ 2 )})
  
  mse
}

pdf('mse_lambdas.pdf')
plot(ls, plot.level.lambda(data, ls))
dev.off()