library('MASS')

# give a sample given the number of regressors and the size


sigma <- diag(p + 1)

sample_sn <- function(p, n, N) {
  array(rnorm(n * (p + 1) * N), dim = c(n, p + 1, N))
}

# transform so that there's dependence
sample <- function(srn, sigma_chol) {
  s <- apply(srn, 3, function(x) {x %*% sigma_chol}, simplify = FALSE)
  lapply(s, function(x) 
    {cbind(x[, 1] - x[, 2] + x[, dim(srn)[2]], x[, 1:(dim(srn)[2] - 1)])})
}

sample_r <- function(full_sample, sigma, p) {
  c <- chol(sigma)
  sp <- full_sample[, c(1:p, dim(full_sample)[2])]
  spc <- sp %*% c
  cbind(spc[, 1] - spc[, 2] + spc[, p + 1], spc[, 1:p])
}


p <- 90
n <- 100

N <- 10000 # number of samples

rns <- sample_sn(p, n, N)
mc <- sample(rns, diag(p + 1))

# b -- bivariate ols

# calculating ols coefs

ols <- function(sample, nreg) {
  ginv(sample[, 2:(nreg+1)]) %*% sample[, 1]
}

mv <- function(ptilde, samples) {
  bols <- sapply(samples, function(x){ols(x, ptilde)[1]})
  c(mean(bols), var(bols))
}

table_b <- mv(1, mc)

print(table_b)

# c -- being jittery about how many regressors to use

p_tries <- c(5, 10, 50, 85, 90)
tries <- length(p_tries) # number of 'guesses' for 
                        # how many regressors to include
table_c <- matrix(0, 2, tries)

for (i in 1:tries) {
  table_c[, i] <- mv(p_tries[i], mc)
}

print(table_c)



# d -- predicting a variable with regressors that are orthogonal to it lol

p_tries_d <- c(5, 10, 50, 85, 90)
tries_d <- length(p_tries_d)

table_d <- matrix(0, 2, tries_d)

min_eig_pred_er <- function(ptilde, samples) {
  eig <- sapply(samples, function(x)
    {min(eigen(1 / dim(x)[1] * t(x[, 3:(ptilde+1)]) %*% 
                 x[, 3:(ptilde+1)])$values )})
  #proj <- lapply(samples, function(x){lm(x[, 2] ~ x[, 3:(ptilde+1)] -1)})
  #sq_pred <- sapply(proj, function(x) mean(x$fitted.values^2))
  
  proj <- lapply(samples, function(x)
    {x[, 3:(ptilde+1)] %*% ols(x[, -1], ptilde - 1)})
  sq_pred <- sapply(proj, function(x) mean(x^2))
  
  c(mean(eig), mean(sq_pred))
}

for (i in 1:tries_d) {
  table_d[, i] <- min_eig_pred_er(p_tries_d[i], mc)
}

print(table_d)



# e -- introducing dependence 

N_tries <- c(100, 200, 500, 1000)
tries_e <- length(N_tries)

table_e <- array(0, dim = c(3, tries_e, 3))

N_samples <- 1000

rho_e = c(0, 0.5, 0.9)
p_e <- 90

# a function that gives out a covariance matrix given \rho
sig <- function(rho, p) {
  sig <- matrix(rho, p + 1, p + 1)
  sig[(p+1), ] <- rep(0, p + 1)
  sig[, (p+1)] <- rep(0, p + 1)
  diag(sig) <- rep(1, p + 1)
  sig
}

sigma_e <- lapply(as.list(rho_e), sig, p = p_e) # store covariance matrices 
#    not to make them anew each time they're called
sigma_chol_e <- lapply(sigma_e, chol)


# a function that returns mean and variance of beta_1 and 
# the mean lowest eigenvalue
e_func <- function(sample, p) {
  beta_bols <- sapply(sample, function(x){ols(x, p)[1]})
  
  eig <- sapply(sample, function(x)
    {min(eigen(1 / dim(x)[1] * t(x[, 2:(p+1)]) %*% x[, 2:(p+1)])$values )})
  
  c(mean(beta_bols), var(beta_bols), mean(eig))
}

eig_rho <- vector('list', length = length(rho_e))
eig_big_rho <- vector('list', length = length(rho_e))

for (j in 1:length(rho_e)) {
  
  mc_e <- vector('list', length = tries_e)
  for (i in 1:tries_e) {
    mc_e[[i]] <- sample(srn[1:N_tries[i], 1:(p_e+1), ], sigma_chol_e[[j]])
    
    table_e[, i, j] <- e_func(mc_e[[i]], p_e)
    
    
    # 1 step ahead -- to avoid remaking samples again for the next question
    if (N_tries[i] == 1000) {
      all_eig <- lapply(mc_e[[i]], function(x) 
        {sort(eigen(1 / dim(x)[1] * t(x[, 2:(p_e+1)]) %*% 
                      x[, 2:(p_e+1)])$values)})
      almost_all_eig <- lapply(all_eig, function(x) {x[-p_e]})
      biggest_eig <- lapply(all_eig, function(x) {x[p_e]})
      eig_rho[[j]] <- colMeans(do.call(rbind, almost_all_eig))
      eig_big_rho[[j]] <- mean(unlist(biggest_eig))
      
      rm(all_eig, almost_all_eig, biggest_eig)
    }
    
  }
  
  rm(mc_e)
}



print(table_e)

# f - how do other eigenvalues behave


eig_rho <- vector('list', length = length(rho_e))
eig_big_rho <- vector('list', length = length(rho_e))
for (i in 1:length(rho_e)) {
  all_eig <- lapply(mc_e[[i]], function(x) 
    {sort(eigen(1 / dim(x)[1] * t(x[, 2:(p+1)]) %*% x[, 2:(p+1)])$values)})
  almost_all_eig <- lapply(all_eig, function(x) {x[-p]})
  biggest_eig <- lapply(all_eig, function(x) {x[p]})
  eig_rho[[i]] <- colMeans(do.call(rbind, almost_all_eig))
  eig_big_rho[[i]] <- mean(unlist(biggest_eig))
}

x_rho <- c(sapply(rho_e, function(x) rep(x, p_e - 1)))
y_rho <- unlist(eig_rho)
plot(x_rho, y_rho)

print(eig_big_rho)


# g

# i guess we need to resample with a higher number of regressors ? 
#     because 20 * log N is above 90 for any N we're considering here ...

p_g <- 900 # should be enough lol

sigma_g <- lapply(as.list(rho_e), sig, p = p) # store covariance matrices 
                                # not to make them anew each time they're called

srn <- sample_rn(max(N_tries), p_g, N_samples) # this will blow up your memory!!


table_g <- array(0, dim = c(6, tries_e, 3))

coefs_eig <- vector('list', length = length(rho_e))
coefs_eig <- lapply(coefs_eig, function(x) {rep(list(matrix(0, 4, N_samples)), 
                                                tries_e)})



for (j in 1:length(rho_e)) {
  for (s in 1:N_samples) {
    sample <- sample_r(srn[, , s], sigma_g[[j]], p_g)
    
    for (i in 1:tries_e) {
      ss <- 0.9 * N_tries[i]
      slog <- round(20 * log(N_tries[i]))
      coefs_eig[[j]][[i]][1, s] <- ols(sample[1:N_tries[i], ], ss)[1]
      coefs_eig[[j]][[i]][2, s] <- min(eigen(1 / N_tries[i] * 
                                        t(sample[1:N_tries[i], 2:(ss+1)]) %*% 
                                        sample[1:N_tries[i], 2:(ss+1)])$values )
      coefs_eig[[j]][[i]][3, s] <- ols(sample[1:N_tries[i], ], slog)[1]
      coefs_eig[[j]][[i]][4, s] <- min(eigen(1 / N_tries[i] * 
                                        t(sample[1:N_tries[i], 2:(slog+1)]) %*% 
                                      sample[1:N_tries[i], 2:(slog+1)])$values )
    }
    
  }

}


for (j in 1:length(rho_e)) {
  for (i in 1:tries_e) {
    table_g[, i, j] <- c(mean(coefs_eig[[j]][[i]][1, ]), 
                         var(coefs_eig[[j]][[i]][1, ]),
                         mean(coefs_eig[[j]][[i]][2, ]),
                         mean(coefs_eig[[j]][[i]][3, ]),
                         var(coefs_eig[[j]][[i]][3, ]),
                         mean(coefs_eig[[j]][[i]][4, ]))
  }
}


# save all of the tables for reporting

save(table_b, table_c, table_d, table_e, eig_rho, eig_big_rho, table_g, 
     file = 'output.RData')