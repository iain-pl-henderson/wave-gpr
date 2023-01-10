selectPtsToPredict <- function(X, x0,R,c){
  # X : matrix of size n x 4  (x, y, z, t)
  #browser()
  n <- nrow(X)
  X[,1:3] <- X[,1:3] - matrix(x0, nrow = nrow(X), ncol = 3, byrow = TRUE)

  return( which( abs(sqrt(X[,1]^2 + X[,2]^2 + X[,3]^2) - c*X[,4])/R < 1 ) )
}

# S <- abs(sqrt(X[,1]^2 + X[,2]^2 + X[,3]^2) - c*X[,4])/R
# return(which(S < 1))