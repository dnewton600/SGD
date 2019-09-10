# Code for DRA/IRA

# At each iteration:
#   1. Get pre-specified number of samples from the oracle, define ERF
#   2. Optimize the ERF until sample gradient is below some threshold

# DRA vs. IRA
# For DRA, we keep all previous samples
# For IRA, we obtain new samples


# --- Function for using a dataset --- #

RA <- function(initial_x,
               oracle_dataset,
               RA_type = 'DRA',
               gradient_function, # function to calculate sample gradient
               objective_function, # function to calculate sample objective
               solver,
               q = 1.1, # sample size increase factor
               starting_n = 10,
               tol_0 = 1,
               num_iters = 10) {    
  
  # initialize x and sample size
  x <- initial_x
  n <- starting_n/q
  
  # initialize the local dataset
  dataset <- list()
  dataset$X <- matrix(0,nrow = 0, ncol = length(x))
  dataset$y <- matrix(0, nrow = 0, ncol = 1)
  
  # record keeping for plotting
  x_record <- matrix(0, nrow = length(x), ncol = num_iters + 1)
  x_record[,1] <- x
  work_record <- rep(0,num_iters + 1)
  
  # run RA
  for(ii in 1:num_iters) {

    # 1. Choose sample size
    n <- ceiling(q*n)
    
    # 2. Update sample-path problem (i.e. query new observations)
    if(RA_type == "DRA") {
      work_record[ii + 1] <- work_record[ii] + (n-length(dataset$y))
      dataset <- queryOracle(dataset, n - length(dataset$y),oracle_dataset)
      
    } else if(RA_type == "IRA") {
      work_record[ii + 1] <- work_record[ii] + n
      dataset <- list()
      dataset$X <- matrix(0,nrow = 0, ncol = length(x))
      dataset$y <- matrix(0, nrow = 0, ncol = 1)
      dataset <- queryOracle(dataset,n,oracle_dataset) # in other file
    }
    
    # 3. Optimize the SPP using the solver() up to tolerance tol
    tol <- tol_0/sqrt(n)
    
    x <- solver(starting.x = x, 
                tol = tol, 
                obj.function = function(x) objective_function(x,dataset$X, dataset$y), 
                grad.function = function(x) gradient_function(x,dataset$X, dataset$y) )
    
    # update records for plotting
    x_record[,ii+1] <- x
  }
  
  return(list(x_record, work_record))
}







# --- Polyak on synthetic LS --- #

Polyak_ls <- function(initial_x,
                      true_beta,
                      gradient_function,
                      enlarge_dataset,
                      initial_step_size = 1,
                      num_iters) {
  
  x <- initial_x 
  x_avg <- x
  alpha <- initial_step_size
  error_record <- rep(0,num_iters + 1)
  error_record[1] <- sqrt(sum( (x - true_beta)^2 ))
  
  for(ii in 1:num_iters) {
    dataset <- list()
    dataset$X <- matrix(0,nrow = 0, ncol = length(x))
    dataset$y < matrix(0, nrow = 0, ncol = 1)
    dataset <- enlarge_dataset(dataset,1,true_beta)
    x <- x - (alpha/sqrt(ii))*gradient_function(x,dataset$X,dataset$y)
    x_avg <- (x_avg*ii + x)/(ii + 1)
    error_record[ii + 1] <- sqrt(sum( (x_avg - true_beta)^2 ))
  }
  
  return(list(error_record, 1:(num_iters+1) ))
}
