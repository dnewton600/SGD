# This code implements several SGD variants for logistic regression.
# It can be applied to any dataset that has the following fields:

# dataset$x : n x p data matrix
# dataset$y : vector of y observations
# dataset$n : sample size
# dataset$counter : used to keep track of sampling. Set to 1 initially.

g_n_fast <- function(sample_size, beta, dataset, lambda) {
# Function to compute gradients
# sample_size: size of subsample used to compute gradient
# beta: the vector of paramaters at which to compute the gradient
# lambda: the L2 regularization parameter
  
  # First, decide whether the counter will exceed the dataset size
  remainder <- (dataset$counter + sample_size - 1) - dataset$n
  if(remainder <= 0) {
    indices <- dataset$counter:(dataset$counter + sample_size - 1)
  }
  else {
    indices <- c( dataset$counter:dataset$n, 1:(remainder - 1) )
  }
  
  # Get X and y subsamples
  X_ss <- dataset$x[indices,,drop = FALSE]
  y_ss <- dataset$y[indices]
  
  # vector of predicted probabilities
  p <- predict(X_ss, beta)
  
  return(-as.vector((1/sample_size)*t(X_ss)%*%(y_ss - p)) + lambda*beta)
}

predict <- function(X,beta) {
  # X: n x p matrix
  # beta: vector of parameters
  # returns vector of probabilities
  return(as.vector(1/(1 + exp(- (X%*%beta)))))
}

#### Functions for Newton CG method ####
newton_direction <- function(beta, dataset, n_hess, grad, cg_iters, n_grad, lambda) {
  # Uses truncated CG to find approximate newton step
  remainder <- (dataset$counter + n_grad - 1) - dataset$n
  if(remainder <= 0) {
    indices <- dataset$counter:(dataset$counter + n_grad - 1)
  } else {
    indices <- c( dataset$counter:dataset$n, 1:(remainder) )
  }
  
  X_ss <- as.matrix(dataset$x[sample(indices,n_hess),,drop = FALSE])
  
  # matrix of weights
  p <- predict(X_ss, beta)
  W <- diag(p*(1-p))
  
  # now we perform CG
  return(CG(Ax = Hp, X_ss, W, b = -grad, lambda, num_iters = cg_iters))
}

Hp <- function(X,W,p,lambda) {
  # Returns Hessian-vector product
  # X is the n x p data matrix
  # W is the outputed W from newton_direction()
  # p is a p + 1 dimensional vector
  # 
  as.vector( ( (1/nrow(X) )*(t(X)%*%W%*%X) + diag(rep(lambda, length(p))) )  %*%p)
}

CG <- function(Ax, X, W, b, lambda, num_iters) {
  # Implements the conjugate gradient method. 
  # Ax is a function that takes x and returns Ax
  # b is the target vector
  # number of iterations should be <= length of b

  x <- rep(0, length(b))
  r <- b - Ax(X,W,x,lambda)
  d <- r
  
  for(ii in 1:num_iters) { 
    AxXWd <- Ax(X,W,d,lambda)
    alpha <- sum(r^2)/(sum(d*AxXWd))
    x <- x + alpha*d
    r_prev <- r
    r <- r - alpha*AxXWd
    beta <- sum(r^2)/sum(r_prev^2)
    d <- r + beta*d
  }
  return(x)
}

#########

#### Algorithms ####
logistic_SGD <- function(num_iters, dataset, start = 1, log_like_mat = FALSE, record_every = 10, step_size_start = 1, init_dig = 0, lambda = 1) {
  
  # initialize solution; set counter to start
  beta <- rep(init_dig, ncol(dataset$x))
  dataset$counter <- start
  
  # Create matrix to keep track of the -log likelihood as the algorithm runs.
  # This will be used for plotting the path.
  log_like <- matrix(0, ncol = 2, nrow = num_iters/record_every + 1)
  index_counter <- 0
  log_like[1,1] <- index_counter
  log_like[1,2] <- full_likelihood(beta, dataset, lambda)
  index_counter <- index_counter + 1
  
  # Run the algorithm for num_iters 
  for(ii in 1:num_iters) {
    if(ii %% 100 == 1) print(ii)
    alpha <- step_size_start/ii
    beta <- beta - alpha*g_n_fast(1, beta, dataset, lambda)
    dataset$counter <- (dataset$counter %% dataset$n) + 1
    
    # Record the log-likelihood as often as is requested by the user
    if( (log_like_mat == TRUE) && (ii %% record_every == 0) )  {
      index_counter <- index_counter + 1
      log_like[index_counter,1] <- ii
      log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
    }
  }
  
  return(list(beta_hat = beta,log_like = log_like))
}

logistic_batch_SGD <- function(num_iters, dataset, sample_size, start = 1, log_like_mat = FALSE, record_every = 1, step_size_start = 1, init_dig = 0, lambda = 1) {
  
  beta <- rep(init_dig,ncol(dataset$x))
  dataset$counter <- start
  log_like <- matrix(0, ncol = 2, nrow = num_iters/record_every + 1)
  index_counter <- 0
  log_like[1,1] <- index_counter
  log_like[1,2] <- full_likelihood(beta, dataset, lambda)
  index_counter <- index_counter + 1
  
  for(ii in 1:num_iters) {
    if(ii %% 100 == 1) print(ii)
    alpha <- step_size_start/ii
    beta <- beta - alpha*g_n_fast(sample_size, beta, dataset, lambda)
    dataset$counter <- (dataset$counter + sample_size) %% dataset$n
    
    if( (log_like_mat == TRUE) && (ii %% record_every == 0) )  {
      index_counter <- index_counter + 1
      log_like[index_counter,1] <- ii*sample_size
      log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
    }

  }
  return(list(beta_hat = beta,log_like = log_like))
}

momentum_SGD <- function(num_iters, dataset, start = 1, log_like_mat = TRUE, momentum = 1, record_every = 10, alpha = 1, init_dig = 0, lambda = 1) {
  beta <- rep(init_dig, ncol(dataset$x))
  beta_prevs <- list(beta,beta)
  dataset$counter <- start
  
  log_like <- matrix(0, ncol = 2, nrow = num_iters/record_every + 1)
  index_counter <- 0
  log_like[1,1] <- index_counter
  log_like[1,2] <- full_likelihood(beta, dataset, lambda)
  index_counter <- index_counter + 1
  
  for(ii in 1:num_iters) {
    if(ii %% 100 == 1) print(ii)
    beta <- beta - (alpha/ii)*g_n_fast(1, beta, dataset, lambda) + momentum*(beta - beta_prevs[[2]]) 
    beta_prevs[[2]] <- beta_prevs[[1]]
    beta_prevs[[1]] <- beta
    
    dataset$counter <- (dataset$counter %% dataset$n) + 1
    
    if( (log_like_mat == TRUE) && (ii %% record_every == 0) )  {
      index_counter <- index_counter + 1
      log_like[index_counter,1] <- ii
      log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
    }
  }
  
  return(list(beta_hat = beta,log_like = log_like))
}

newton_CG <- function(num_iters, dataset, sampling, theta = 1, start = 1, log_like_mat = FALSE, n_grad, prop_hess, cg_iters, init_dig = 0, init_beta = NULL, alpha = 1, lambda) {
    # num iters is number of iterations of the outer loop of the newton algorithm
    # dataset is a dataset with fields x, y, counter, and n
    # sampling is 'standard', 'norm', or 'ip'
    # start is the starting index of the dataset from which we are sampling
    # log_like_mat is a boolean that records whether or not to track progress of the algorithm (for plotting)
    # record_every describes how often the likelihood is evaluated (for plotting)
    # If sampling == 'standard,' n_grad is the number of samples used to estimate the gradient
    # prop_hess is the proportion of samples of the gradient to be used in estimating the Hessian
    # cg_iters is the number of conjugate gradient iterations used to approximate the newton system of equations
  
    #bookkeeping <- matrix(0, nrow = 2000, ncol = 3)
   # colnames(bookkeeping) <- c("counter","var_est","norm_sq")
  
    if(!is.null(init_beta)) {
      beta <- init_beta
    }
    else {
      beta <- rep(init_dig,ncol(dataset$x))
    }
    dataset$counter <- start
    
    grad_n_lower <- 100
    grad_n_upper <- min(2000, dataset$n)
    max_n_hess <- 150
    
    log_like <- matrix(0, ncol = 2, nrow = num_iters + 1)
    index_counter <- 0
    work_counter <- 0
    log_like[1,1] <- index_counter
    log_like[1,2] <- full_likelihood(beta, dataset, lambda)
    index_counter <- index_counter + 1
    
    for(ii in 1:num_iters) {
      cat("\n iteration number:", ii,"\n")
      cat("counter:", dataset$counter, "\n")
      
    if(sampling == 'norm') {
      grad_1 <- g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter + 1) %% dataset$n
        
      grad_2 <- g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter + 1) %% dataset$n
        
      Grads <- matrix(0, nrow = length(beta), ncol = 900)
      Grads[,1] <- grad_1
      Grads[,2] <- grad_2
        
      curr_sample_size <- 2 
        
      grad_est <- rowSums(Grads)/curr_sample_size 
      avg_norm_sq <- sum(colSums( (Grads[,1:curr_sample_size] - grad_est)^2 ) )/(curr_sample_size - 1)
        
      var_est <- avg_norm_sq / curr_sample_size
        
      while( (var_est > theta^2*avg_norm_sq^2 || curr_sample_size < grad_n_lower) && curr_sample_size < grad_n_upper) {
        
        curr_sample_size <- curr_sample_size + 1
        new_grad <- g_n_fast(1, beta, dataset, lambda)
        dataset$counter <- (dataset$counter + 1) %% dataset$n
          
        if(curr_sample_size %% 900 == 1) {
          Grads <- cbind(Grads, matrix(0, nrow = length(beta), ncol = 900))
        }
          
        Grads[,curr_sample_size] <- new_grad
        grad_est <- rowSums(Grads)/curr_sample_size 
        avg_norm_sq <- sum(colSums( (Grads[,1:curr_sample_size] - grad_est)^2 ) )/(curr_sample_size - 1)
        var_est <- avg_norm_sq / curr_sample_size
        
      }
      
      grad <- grad_est
      n_hess <- min( round(prop_hess*curr_sample_size), max_n_hess)
      grad_sample_size <- curr_sample_size
        
      } else if(sampling == 'ip') {
        
        grad_1 <- g_n_fast(1, beta, dataset, lambda)
        dataset$counter <- (dataset$counter + 1) %% dataset$n
        
        grad_2 <- g_n_fast(1, beta, dataset, lambda)
        dataset$counter <- (dataset$counter + 1) %% dataset$n
        
        Grads <- matrix(0, nrow = length(beta), ncol = 900)
        Grads[,1] <- grad_1
        Grads[,2] <- grad_2
        
        curr_sample_size <- 2 
        
        grad_est <- rowSums(Grads)/curr_sample_size 
        dot_products <- Grads[,1:curr_sample_size]*grad_est
        norm_sq <- sum(grad_est^2)
        
        var_est <- ((sum( (dot_products - norm_sq)^2 ))/(curr_sample_size - 1))/curr_sample_size
        
        while( (var_est > theta^2*norm_sq^2 || curr_sample_size < grad_n_lower) && curr_sample_size < grad_n_upper) {
    
          # book keeping
          # bookkeeping[dataset$counter,] <- c(dataset$counter, var_est, theta^2*norm_sq^2)
          
          curr_sample_size <- curr_sample_size + 1
          new_grad <- g_n_fast(1, beta, dataset, lambda)
          dataset$counter <- (dataset$counter + 1) %% dataset$n
          
          if(curr_sample_size %% 900 == 1) {
            Grads <- cbind(Grads, matrix(0, nrow = length(beta), ncol = 900))
          }
          
          Grads[,curr_sample_size] <- new_grad
          grad_est <- rowSums(Grads)/curr_sample_size 
          dot_products <- Grads[,1:curr_sample_size]*grad_est
          norm_sq <- sum(grad_est^2)
          var_est <- ((sum( (dot_products - norm_sq)^2 ))/(curr_sample_size - 1))/curr_sample_size
          
        }
        
        grad <- grad_est
        n_hess <- min( round(prop_hess*curr_sample_size), max_n_hess)
        grad_sample_size <- curr_sample_size
        
      } else { # sampling == 'standard'
        grad <- g_n_fast(n_grad, beta, dataset, lambda)
        n_hess <- round(prop_hess*n_grad)
        grad_sample_size <- n_grad
        dataset$counter <- (dataset$counter + grad_sample_size) %% dataset$n
      }
      
      print("grad_sample_size")
      print(grad_sample_size)
  
      nd <- newton_direction(beta, dataset, n_hess, grad, cg_iters, n_grad = grad_sample_size, lambda)
      
      # line search
      beta_new <- beta + alpha*nd
      step <- alpha
      while( (full_likelihood(beta_new,dataset,lambda) > full_likelihood(beta, dataset,lambda) ) && (step > .0005) ){
        step <- step*.7
        beta_new <- beta + step*nd
      }
      beta <- beta_new
      
      #print(full_likelihood(beta,dataset))
      
      if(log_like_mat == TRUE) {
        index_counter <- index_counter + 1
        work_counter <- work_counter + grad_sample_size + n_hess*cg_iters
        log_like[index_counter,1] <- work_counter
        log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
      }
    }
    
    return(list(beta_hat = beta, log_like = log_like )) #bookkeeping = bookkeeping))
}

logistic_SGD_adpt <- function(num_iters, dataset, sampling, theta = 1, start = 1, log_like_mat = FALSE, init_dig = 0, init_beta = NULL, step_size_start = 1, lambda) {
  # num iters is number of iterations of the outer loop of the newton algorithm
  # dataset is a dataset with fields x, y, counter, and n
  # sampling is 'norm' or 'ip'
  # start is the starting index of the dataset from which we are sampling
  # log_like_mat is a boolean that records whether or not to track progress of the algorithm (for plotting)
  # record_every describes how often the likelihood is evaluated (for plotting)

  #bookkeeping <- matrix(0, nrow = 2000, ncol = 3)
  # colnames(bookkeeping) <- c("counter","var_est","norm_sq")
  
  if(!is.null(init_beta)) {
    beta <- init_beta
  }
  else {
    beta <- rep(init_dig,ncol(dataset$x))
  }
  dataset$counter <- start

  # Min and max for adaptive sampling size

  
  # Book keeping
  log_like <- matrix(0, ncol = 2, nrow = num_iters + 1)
  index_counter <- 0
  work_counter <- 0
  log_like[1,1] <- index_counter
  log_like[1,2] <- full_likelihood(beta, dataset, lambda)
  index_counter <- index_counter + 1
  
  # Main iteration
  for(ii in 1:num_iters) {
    cat("\n iteration number:", ii,"\n")
    cat("counter:", dataset$counter, "\n")
    
    grad_n_lower <- 2 + (1/num_iters)*(dataset$n - 2)
    grad_n_upper <- min(2000, dataset$n)
    
    if(sampling == 'norm') {
      grad_1 <- g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter + 1) %% dataset$n
      
      grad_2 <- g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter + 1) %% dataset$n
      
      
      sum_gi <- grad_1 + grad_2
      sum_gisq <- grad_1^2 + grad_2^2
      
      curr_sample_size <- 2 
      
      grad_est <- sum_gi/curr_sample_size 
      samp_variance <-  sum ( (sum_gisq - ( (sum_gi)^2/curr_sample_size ) )/ (curr_sample_size - 1) )
      var_est <- samp_variance / curr_sample_size
      
      while( (var_est > theta^2*(sum(grad_est^2)) || curr_sample_size < grad_n_lower) && curr_sample_size < grad_n_upper) {
        
        curr_sample_size <- curr_sample_size + 1
        new_grad <- g_n_fast(1, beta, dataset, lambda)
        dataset$counter <- (dataset$counter + 1) %% dataset$n
        
        sum_gi <- sum_gi + new_grad
        sum_gisq <- sum_gisq + new_grad^2

        grad_est <- sum_gi/curr_sample_size 
        samp_variance <- sum ( (sum_gisq - ( (sum_gi)^2/curr_sample_size ) )/ (curr_sample_size - 1) )
        var_est <- samp_variance / curr_sample_size
        
      }
      
      grad <- grad_est
      grad_sample_size <- curr_sample_size
      
    } else if(sampling == 'ip') {
      
      grad_1 <- g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter + 1) %% dataset$n
      
      grad_2 <- g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter + 1) %% dataset$n
      
      Grads <- matrix(0, nrow = length(beta), ncol = 900)
      Grads[,1] <- grad_1
      Grads[,2] <- grad_2
      
      curr_sample_size <- 2 
      
      grad_est <- rowSums(Grads)/curr_sample_size 
      dot_products <- Grads[,1:curr_sample_size]*grad_est
      norm_sq <- sum(grad_est^2)
      
      var_est <- ((sum( (dot_products - norm_sq)^2 ))/(curr_sample_size - 1))/curr_sample_size
      
      while( (var_est > theta^2*norm_sq || curr_sample_size < grad_n_lower) && curr_sample_size < grad_n_upper) {
        
        # book keeping
        # bookkeeping[dataset$counter,] <- c(dataset$counter, var_est, theta^2*norm_sq^2)
        
        curr_sample_size <- curr_sample_size + 1
        new_grad <- g_n_fast(1, beta, dataset, lambda)
        dataset$counter <- (dataset$counter + 1) %% dataset$n
        
        if(curr_sample_size %% 900 == 1) {
          Grads <- cbind(Grads, matrix(0, nrow = length(beta), ncol = 900))
        }
        
        Grads[,curr_sample_size] <- new_grad
        grad_est <- rowSums(Grads)/curr_sample_size 
        dot_products <- Grads[,1:curr_sample_size]*grad_est
        norm_sq <- sum(grad_est^2)
        var_est <- ((sum( (dot_products - norm_sq)^2 ))/(curr_sample_size - 1))/curr_sample_size
        
      }
      
      grad_sample_size <- curr_sample_size
      
    } else { # sampling == 'standard'
      print("Select valid sampling option.")
    }
    
    print("grad_sample_size")
    print(grad_sample_size)
  
    
    # line search
    beta <- beta - (step_size_start)*grad_est
    
    #print(full_likelihood(beta,dataset))
    
    if(log_like_mat == TRUE) {
      index_counter <- index_counter + 1
      work_counter <- work_counter + grad_sample_size
      log_like[index_counter,1] <- work_counter
      log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
    }
    
    else {
      work_counter <- work_counter + grad_sample_size
    }
    
    if(work_counter > 20000) {
      break
    }
  }
  
  print("total work")
  print(work_counter)
  
  return(list(beta_hat = beta, log_like = log_like )) #bookkeeping = bookkeeping))
}

ssn_svrg <- function(num_iters, dataset, start = 1, log_like_mat = FALSE, step_size_start, n_sgd, n_grad, prop_hess, cg_iters, init_dig = 0, alpha, lambda) {
  
  beta <- rep(init_dig, ncol(dataset$x))
  dataset$counter <- start
  
  log_like <- matrix(0, ncol = 2, nrow = num_iters*2 + 1)
  work <- 0
  index_counter <- 0
  log_like[1,1] <- index_counter
  log_like[1,2] <- full_likelihood(beta, dataset, lambda)
  index_counter <- index_counter + 1
  
  for(ii in 1:num_iters) {
    
    # sgd step
    for(jj in 1:n_sgd) {
      step_size <- step_size_start/((ii - 1)*n_sgd + jj)
      beta <- beta - step_size*g_n_fast(1, beta, dataset, lambda)
      dataset$counter <- (dataset$counter %% dataset$n) + 1
    }
    
    if( log_like_mat == TRUE)  {
      index_counter <- index_counter + 1
      work <- work + n_sgd
      log_like[index_counter,1] <- work
      log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
    }
    
    # newton_step
    beta <- (newton_CG(1, dataset, sampling = "standard", theta = NA, start = dataset$counter, n_grad = n_grad, prop_hess = prop_hess, cg_iters = cg_iters, init_beta = beta, alpha = alpha, lambda = lambda))$beta_hat
    dataset$counter <- (dataset$counter + n_grad + n_grad*prop_hess*cg_iters) %% dataset$n
    
    if( log_like_mat == TRUE)  {
      index_counter <- index_counter + 1
      work <- work + (n_grad + n_grad*prop_hess*cg_iters)
      log_like[index_counter,1] <- work
      log_like[index_counter,2] <- full_likelihood(beta, dataset, lambda)
    }
  }
  
  return(list( beta_hat = beta, log_like = log_like))
  
}

#### Model evaluation ####

full_likelihood <- function(beta, dataset, lambda) {
  # Returns the negative log likelihood evaluated at beta with L2 parameter lambda
  return( -(1/dataset$n)*sum(dataset$y * (dataset$x %*% beta) - log(1 + exp(dataset$x %*% beta)) ) + .5*lambda*sum(beta^2) )
}

training_error <- function(beta, dataset) {
# Returns proportion of images classified correctly
  return(sum( (log(1/(1 + exp(-dataset$x %*% beta )) ) > log(0.5)) == dataset$y)/dataset$n)
}

plot(log(basic_SGD$log_like[,2])~basic_SGD$log_like[,1])
plot(log(batch$log_like[,2])~batch$log_like[,1])
plot(log(nesterov$log_like[,2])~nesterov$log_like[,1])
plot(log(newton_CG$log_like[,2])~newton_CG$log_like[,1])

#### Plotting with ggplot ####
library(ggplot2)

plot_all <- function(basic, batch, momentum, newton) {
  
  my_list <- list(basic, batch, momentum, newton)
  
  for(ii in my_list) {
    print(str(ii))
    print(ii$log_like[nrow(ii$log_like),1])
  }
  
  x_axis_length <- nrow(basic$log_like)
  
  likelihood <-  c(basic$log_like[1:x_axis_length,2], batch$log_like[1:x_axis_length,2], momentum$log_like[1:x_axis_length,2], newton$log_like[1:x_axis_length,2])
  type <- c(rep("Basic", x_axis_length), rep("Batch", x_axis_length), rep("Momentum", x_axis_length), rep("Newton", x_axis_length))
  my.df <- data.frame(iter = rep(basic$log_like[,1],4), Algorithm = type, logloglikelihood = log(likelihood))
  ggplot(my.df, aes(x = iter, y = logloglikelihood, col = Algorithm)) + 
    geom_line() +
    ggtitle("MNIST") +
    theme(plot.title = element_text(hjust = 0.5)) +
    xlab("Work Complexity") + 
    ylab("log-log-likelihood")
}

plot_all_2 <- function(results, names, title, upper_x_limit) {
# Plots results from running each algorithm.
# results: output from the SGD functions above
# names: names of the algorithms (used for the legend in the plots)
# title: title of the plot
# upper_x_limit: How far the x axis should stretch
#   (this is a parameter because we don't know ahead of time
#    how much sampling the adaptive methods will need)
  full_data <- NULL
  for(ii in 1:length(results)) {
    full_data <- rbind(full_data, cbind(results[[ii]]$log_like, ii) )
  }
  
  colnames(full_data) <- c("work","like","type")
  mydf <- as.data.frame(full_data)
  
  return(ggplot(mydf, aes(x = work, y = log(like), col = factor(type))) +  
         geom_line() +
         ggtitle(title) +
         theme(plot.title = element_text(hjust = 0.5)) +
         xlab("Work Complexity") + 
         ylab("log-log-likelihood") + 
         scale_color_hue(labels = names) +
         scale_x_continuous(limits = c(0, upper_x_limit)))
  
}
