# Functions that compute gradient and hessian
#   for fitting logistic regression model at a particular solution.

log_reg_obj_function <- function(beta,X,y) {
  # negative log likelihood
  sample.size <- nrow(X)
 -(1/sample.size)*sum(y*(X%*%beta) - log(1 + exp(X%*%beta)))
}

log_reg_grad_function <- function(beta,X,y) {
  sample_size <- nrow(X)
  p <- as.matrix(predict(X, beta),ncol = 1)
  return(-as.vector((1/sample_size)*t(X)%*%(y - p)) + beta)
}

log_reg_hess_function <- function(beta,X,y) {
  sample_size <- nrow(X)
  p <- predict(X, beta)
  W <- diag(x = p*(1 - p) )
  return((1/sample_size)*t(X)%*%W%*%X + diag(x = rep(1,length(beta))) )
}

predict <- function(X,beta) {
  # helper function
  # X is n x p matrix
  # returns vector of probabilities
  return(as.vector(1/(1 + exp(- (X%*%beta)))))
}