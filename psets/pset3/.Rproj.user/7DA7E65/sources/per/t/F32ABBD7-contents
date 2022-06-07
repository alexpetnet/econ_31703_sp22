library(xtable)

### topics in ecm (bonhomme, spring 22) -- pset3 -- sasha petrov ===============

kmeans.objective <- function(K, data, clustering) {
  X <- data[, -1]
  
  cl <- rep(NA, K)
  
  for (i in 1:K) {
    Xk <- X[clustering == i, , drop = FALSE]
    m <- apply(Xk, 2, mean)
    Xkd <- apply(Xk, 1, function(x) {sum( (x - m)^2 )})
    cl[i] <- sum(Xkd)
  }
  
  sum(cl)
}

kmeans.mean.update <- function(K, data, clustering) {
  X <- data[, -1, drop = FALSE]
  p <- ncol(X)
  
  M <- matrix(NA, K, p)
  
  for (i in 1:K) {
    Xk <- X[clustering == i, , drop = FALSE]
    M[i, ] <- apply(Xk, 2, mean)
  }
  
  M
}

choose.cluster <- function(X, means) {
  diff <- t(apply(means, 1, function(x) {x - X}))
  norm <- apply(diff, 1, function(x) {sum(x ^ 2)})
  which.min(norm)
}

kmeans.clustering.update <- function(K, data, means) {
  X <- data[, -1]
  N <- nrow(X)
  
  clustering <- rep(NA, N)
  
  apply(X, 1, function(x) {choose.cluster(x, means)})
}

# a function to compute the update in the means
means.update <- function(means) {
  diff <- means[ , , 1] - means[ , , 2]
  norm <- apply(diff, 1, function(x) {sum(x ^ 2)})
  max(norm)
}

kmeans <- function(K, data, initial.clustering, max = 1000, eps = 1e-6) {
  means <- array(0, dim = c(K, ncol(data[, -1]), 2))
  means[ , , 2] <- kmeans.mean.update(K, data, initial.clustering)
  
  obj <- kmeans.objective(K, data, initial.clustering)
  
  iter <- 1
  
  while((iter <= max) & (means.update(means) > eps)) {
    iter <- iter + 1
    
    stop <- 0 # change to 1 only if max iterations is hit
    
    means[ , , 1] <- means[ , , 2]
    cluster <- kmeans.clustering.update(K, data, means[ , , 2])
    means[ , , 2] <- kmeans.mean.update(K, data, cluster)
    for (j in 1:K) {
      if (NaN %in% means[j, , 2]) {
        means[j, , 2] <- means[j, , 1]
      }
    }
    
    obj <- c(obj, kmeans.objective(K, data, cluster))
    
    if (iter == max) {
      stop <- 1
    }
  }
  
  ans <- list(clustering = cluster,
              means = means[ , , 2],
              objective = obj,
              status = stop)
  
  ans
}

### mc simulations =============================================================

# sample generation function

sample.mc <- function(p = 2, N = 200, 
                   mu = matrix(0, p + 1, 3), 
                   sigma = rep( list( diag(p + 1) ), 3) ) {
  u <- runif(N)
  theta <- rep(1, N) + (u >= .2) + (u >= .5)
  
  s.chol <- lapply(sigma, chol)
  
  s <- t( apply(matrix(theta), 1, function(x) {mu[ , x] + 
      s.chol[[x]] %*% rnorm(p + 1)}) )
  
  matrix(c(s, theta), 200, p + 2)
}

mu <- matrix(c(0, 0, 0, 2, 1, 0, 1, 0, 1), 3)
s <- rep( list( matrix(.1, 3, 3) + diag(c(.9, .1, .1)) ), 3 )

d <- sample.mc(mu = mu, sigma = s)

plot.cluster <- function(K, data, cluster, name) {
  X <- data[ , -1]
  
  pdf(name)
  cl <- rainbow(K)
  plot(X, col = cl[cluster],
       xlab = expression(X[1]), ylab = expression(X[2]))
  l <- apply(matrix(b), 1, function(x) 
    {parse( text = paste0('theta[', as.character(x), ']'))})
  legend("topright", legend = 1:K, 
         col=cl, pch=1)
  dev.off()
}

kmeans(3, d[ , -4], sample(3, 200, replace = TRUE))

r <- rep(list(NA), 10)
r <- lapply(r, function(x) 
  {kmeans(3, d[ , -4], sample(3, 200, replace = TRUE))})
r.o <- sapply(r, function(x) {x$objective[length(x$objective)]})

plot.cluster(3, d[ , -4], d[ , 4], 'true_clustering.pdf')
plot.cluster(3, d[ , -4], r[[7]]$clustering, 'best_k_means.pdf')


## (f)

n.mc <- 1000

samples <- rep(list(NA), n.mc)
samples <- lapply(samples, function(x) {sample.mc(mu = mu, sigma = s)})

k.tilde.cluster <- function(d, kt, mu) {
  r <- kmeans(kt, d, sample(kt, 200, replace = TRUE))
  d <- cbind( d, mu[1, ][d[ , 4]] )
  m <- kmeans.mean.update(kt, d[ , 4:5], r$clustering)
  d <- cbind( d, m[r$clustering] )
  
  bias <- mean( ( d[ , 5] - d[ , 6] ) ^ 2 )
  
  d <- cbind( d, d[ , 1] - d[ , 5] )
  v <- kmeans.mean.update(kt, d[ , 6:7], r$clustering)
  d <- cbind( d, v[r$clustering] )
  
  var <- mean(d[ , 8] ^ 2)
  
  c(bias, var)
}

# a function that returns a list of estimation outcomes given k's to try
k.trial <- function(d, k.tilde = 2:5) {
  o <- as.list(k.tilde)
  lapply(o, function(x) {k.tilde.cluster(d, x, mu)})
}

k.mean.var <- function(kt, mu, sam) {
  t <- sapply(sam, function(x) {k.tilde.cluster(x, kt, mu)})
  apply(t, 1, mean)
}

# the results table
table.f <- apply( matrix(2:5), 1, function(x) {k.mean.var(x, mu, samples)} )