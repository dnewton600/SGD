
queryOracle <- function(dataset,n,oracle_dataset) {
  # Takes n draws from a dataset oracle.
  # Currently allows the possiblity of drawing the same sample twice,
  # which simplifies the code a fair amount
  # and shouldn't be a huge issue for large datasets.
  
  # `dataset` is the working sample used by RA
  # n is the number of requested data points
  # oracle_dataset is the large dataset from which queries are drawn
  
  row_indices <- sample(x = 1:nrow(oracle_dataset$X), size = n)
  dataset$X <- rbind(dataset$X, oracle_dataset$X[row_indices,] )
  dataset$y <- rbind(dataset$y, oracle_dataset$y[row_indices,,drop=FALSE])
  return(dataset)
}

# in a typical OO language the above method
# would be implemented with something like:

#   new_sample <- oracle_dataset.query(n)
#   dataset.update(new_sample)