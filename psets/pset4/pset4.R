### topics in ecm (bonhomme, spring 22) -- pset4 -- sasha petrov ===============

mse_gain <- function(y, m) {
  new <- sum( ( y[1:m] - mean(y[1:m]) ) ^ 2 ) + 
    sum( ( y[-c(1:m)] - mean(y[-c(1:m)]) ) ^ 2 )
  sum( (y - mean(y)) ^ 2 ) - new 
}

tree.grow <- function(data) {
  N <- nrow(data)
  p <- ncol(data) - 1
  
  mse_matrix <- matrix(NA, N, p)
  
  for ( i in 1:p ) {
    s <- data[order(data[ , p + 1]), ]
    mse_matrix[ , i] <- apply(matrix(1:N), 1, 
                              function(x) {mse_gain(s[ , 1], x)},
                              simplify = TRUE)
  }
  
  split <- which(mse_matrix == max(mse_matrix), arr.ind = TRUE)[1, ]
  k <- split[1]
  m <- split[2]
  
  group <- s[ , k + 1] > sort(s[ , k + 1])[m]
  
  c(split, mse_matrix[m, k], group)
}

tree.update <- function(data, split, assign, 
                        min.size.for.split = 5, max.depth = 10) {
  ntl <- nrow(split) + 1
  
  terminal.leafs <- matrix(NA, ntl, 3)
  terminal.leafs[ , 1] <- unique(assign[ , 1])
  terminal.leafs[ , 2] <- apply(terminal.leafs[ , 1, drop = FALSE], 1,
                                function(x) {sum( assign[ , 1] == x )})
  terminal.leafs[ , 3] <- apply(terminal.leafs[ , 1, drop = FALSE], 1,
                                function(x) {assign[assign[ , 1] == x, 2][1]})
  
  if ( all( (terminal.leafs[ , 2] < min.size.for.split) |
           (terminal.leafs[ , 3] > max.depth) ) ) {
    return("gameover")
  }
}