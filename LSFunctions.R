# Functions to compute the gradient and objective function
# for the standard least squares function.

ls_grad_function <- function(beta,X,y) {
  # Returns the sample gradient for a least squares objective (1/2) || X%*%beta - y ||^2
  return( (1/nrow(X))*(t(X)%*%X%*%beta - t(X)%*%y ) )
}

ls_obj_function <- function(beta,X,y) {
  return((1/nrow(X))*sum( (X%*%beta - y)^2 ))
}
